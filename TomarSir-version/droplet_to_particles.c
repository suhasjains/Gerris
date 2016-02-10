/* droplet_to_particles.c
/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2001-2008 Gaurav Tomar
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include "lparticles.h"
#include "refine.h"
#include "adaptive.h"
#include "config.h"

#ifdef HAVE_MPI
#include "mpi_boundary.h"
#endif /*HAVE_MPI*/

typedef struct {
  FttVector pos,vel;
  gdouble vol,density;
} DropletProp;

typedef struct {
  GfsVariable * tag, * c, *t;
  Particle * drops;
  GfsVariable **u;
  guint * sizes;
  guint n, min;
  gdouble resetval;
  gdouble density;
  GfsFunction *fc;
} DropletsPar;


/*Droplet to Particles: Object*/
static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}


static void reset_small_fraction (FttCell * cell, DropletsPar * p)
{

  gint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = p->resetval;
}

static void compute_droplet_properties (FttCell * cell, DropletsPar * p)
{
  g_return_if_fail (cell != NULL);

  gint i = GFS_VALUE (cell, p->tag);
  gdouble h = ftt_cell_size(cell),vol,avgfact;
  FttVector pos; 
  ftt_cell_pos(cell,&pos);
  GfsVariable **u = p->u;

  if (i > 0){
    avgfact = (gdouble) p->sizes[i - 1];
    p->sizes[i - 1]++;


#if FTT_2D
    vol = h*h*GFS_VALUE(cell,p->t);
#else
    vol = h*h*h*GFS_VALUE(cell,p->t);
#endif

    p->drops[i-1].volume += vol*GFS_VALUE(cell,p->c);


    p->drops[i-1].pos.x = (p->drops[i-1].pos.x*avgfact+ pos.x)/p->sizes[i-1];
    p->drops[i-1].pos.y = (p->drops[i-1].pos.y*avgfact+ pos.y)/p->sizes[i-1];
    p->drops[i-1].vel.x = (p->drops[i-1].vel.x*avgfact+ GFS_VALUE(cell,u[0]))/p->sizes[i-1];
    p->drops[i-1].vel.y = (p->drops[i-1].vel.y*avgfact+ GFS_VALUE(cell,u[1]))/p->sizes[i-1];
#if !FTT_2D
    p->drops[i-1].pos.z = (p->drops[i-1].pos.z*avgfact+pos.z)/p->sizes[i-1];
    p->drops[i-1].vel.z = (p->drops[i-1].vel.z*avgfact+GFS_VALUE(cell,u[2]))/p->sizes[i-1];
#endif
  }

}


static void convert_droplets (GfsDomain * domain, 
			      DropletsPar *p, LParticles * lagrangian)
{

  GfsSimulation *sim = gfs_object_simulation(lagrangian); 

  g_return_if_fail (p->c != NULL);

  guint i;

  lagrangian->idlast = lagrangian->maxid;

  for(i = 0; i < p->n; i++){

    if(p->sizes[i] < p->min && p->sizes[i] >4){
      Particle * drop = g_malloc (sizeof(Particle));

      drop->pos = p->drops[i].pos;
      drop->vel = p->drops[i].vel;
      drop->volume = p->drops[i].volume;
      //drop->move = 1;
      drop->cell = gfs_domain_locate(domain, drop->pos, -1);
  
      if(drop->cell){
	drop->id = ++lagrangian->idlast;
	drop->density = sim->physical_params.alpha ? 1./
	  gfs_function_value(sim->physical_params.alpha,drop->cell) : 1.;
	drop->density = p->density;
	lagrangian->particles = g_slist_append(lagrangian->particles,drop);
	//lagrangian->n++;
      }       
    }
  }

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, p);
}

static void tag_sort_droplets (GfsDomain * domain,
			       DropletsPar *p,
			       gint min)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p->c != NULL);

  p->n = gfs_domain_tag_droplets (domain, p->c, p->tag);

  if (p->n > 0 && -min < (gint) p->n) {

    p->sizes = g_malloc0 (p->n*sizeof (guint));

    p->drops = g_malloc0 (p->n*sizeof (Particle));
  /*Compute Droplet Properties here*/
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_droplet_properties, p);

    if (min >= 0)
      p->min = min;
    else {
      guint * tmp = g_malloc (p->n*sizeof (guint));
      memcpy (tmp, p->sizes, p->n*sizeof (guint));
      qsort (tmp, p->n, sizeof (guint), greater);
      /* fixme: this won't work for parallel jobs */
      p->min = tmp[-1 - min];
      g_free (tmp);
    }

  }
}


static void droplet_to_particles_destroy (GtsObject * o)
{
  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->destroy) (o);

  GfsDropletToParticles * d = DROPLET_TO_PARTICLES(o);

  if (d->fc)
    gts_object_destroy (GTS_OBJECT (d->fc));
  
  if (d->fconvert)
    gts_object_destroy (GTS_OBJECT (d->fconvert));
  

  if (d->shape){
    
    if (d->shape->s)
      gts_object_destroy (GTS_OBJECT (d->shape->s));
    
    if (d->shape->m)
      gts_matrix_destroy (d->shape->m); 
    
    g_free (d->shape);
    
  }
}

static void droplet_to_particles_read (GtsObject ** o, GtsFile * fp)
{

  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->read) (o, fp);

  if (fp->type == GTS_ERROR)
    return;

  while(fp->type == '\n')
    gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsDropletToParticles * r = DROPLET_TO_PARTICLES(*o);

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (min)");
    return;
  }

  r->min = atoi (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type == GTS_FLOAT) {
    r->resetwith = atof (fp->token->str);
    gts_file_next_token (fp);
  }

  if (fp->type == GTS_FLOAT) {
    r->density = atof (fp->token->str);
    gts_file_next_token (fp);
  }

  if (fp->type != '\n' && fp->type != GTS_INT) {
    r->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
  }

  if (fp->type == GTS_INT) {
    r->convert = atoi (fp->token->str);
    gts_file_next_token (fp);

    if (fp->type == GTS_INT){
      r->maxlevel = atoi (fp->token->str);
      gts_file_next_token (fp);
    }

    while(fp->type == '\n')
      gts_file_next_token (fp);

    r->shape = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
    gfs_generic_surface_read (r->shape, gfs_object_simulation (*o), fp);
    
  
    if (fp->type != '\n') {
      r->fconvert = gfs_function_new (gfs_function_class (), 0.);
      gfs_function_read (r->fconvert, gfs_object_simulation (r), fp);
    }
  }

}

static void droplet_to_particles_write (GtsObject * o, FILE * fp)
{

  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->write) (o, fp);

  GfsDropletToParticles * r = DROPLET_TO_PARTICLES(o);

  fprintf (fp, " %s %d %lf %lf ", r->c->name, r->min, r->resetwith, r->density);

  if (r->fc){
    gfs_function_write (r->fc, fp);
  }

  if(r->convert != 0){
    fprintf (fp," %d %d ", r->convert, r->maxlevel);

    if (r->shape) {
      fputs (" ( ", fp);
      gfs_function_write (r->shape->f, fp);
      fputs (" ) ", fp);
    }

    if (r->fconvert){
      gfs_function_write (r->fconvert, fp);
    }
  }

}

static void compute_v (FttCell * cell, GfsDropletToParticles * d)
{
  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
}

void save_p_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  GfsVariable * save_prev = data[3];

  //GFS_VARIABLE (cell, save_prev->i) = GFS_VARIABLE (cell, c->i);
  GFS_VALUE (cell, save_prev) = GFS_VALUE (cell, c);
  
	GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i)) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
}

void restore_p_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  gboolean * not_cut = data[1];
  GfsVariable * status = data[2];

  gdouble * resetwith = data[6];

  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i));

  if (solid) {

    GFS_VARIABLE (cell, c->i) =  solid->a;

    g_free (solid);
    *not_cut = FALSE;
  }
  else if (GFS_VARIABLE (cell, status->i) == 0.) {
    /* fixme: this can fail for non-contiguous domains (e.g. non-connected GfsBoxes) */
    g_assert (*not_cut);
    GFS_VARIABLE (cell, c->i) = 0.;
  }
  else {
    g_assert (GFS_VARIABLE (cell, status->i) == 1. || GFS_VARIABLE (cell, status->i) == 2.);
    GFS_VARIABLE (cell, c->i) =  GFS_VARIABLE (cell, status->i) - 1.;
  }
}

void add_to_prev_void (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  GfsVariable * save_prev = data[3];
  Particle *p = (Particle *)data[4];
  GfsDomain *domain = data[5];
  GfsVariable **u = gfs_domain_velocity(domain);

  if(GFS_VALUE (cell, c) > 0.){
    guint j;
    for(j = 0; j < FTT_DIMENSION; j++)
      GFS_VALUE(cell, u[j]) = (&p->vel.x)[j];
  }

  GFS_VARIABLE (cell, c->i) += GFS_VARIABLE (cell, save_prev->i);

  if(GFS_VALUE (cell, c) > 1.0)
    GFS_VARIABLE (cell, c->i) = 1.0;
  else if (GFS_VALUE (cell, c) < 0.0)
      GFS_VARIABLE (cell, c->i) = 0.0;
}

void gfs_domain_assign_fraction (GfsDomain * domain,
			       GfsGenericSurface * s,
			       GfsVariable * c, gdouble resetwith, Particle *p)
{
  gboolean not_cut = TRUE;
  gpointer data[7];
  GfsVariable * status, *save_prev;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (c != NULL);

  status = gfs_temporary_variable (domain);

  save_prev = gfs_temporary_variable (domain);
  data[0] = c;
  data[1] = &not_cut;
  data[2] = status;
  data[3] = save_prev;
  data[4] = p;
  data[5] = domain;
  data[6] = &resetwith;

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) save_p_solid, data);
  GSList * l = g_slist_prepend (NULL, s);

  gfs_domain_init_solid_fractions (domain, l, FALSE, NULL, NULL, status);
  g_slist_free (l);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) restore_p_solid, data);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) add_to_prev_void, data);

  gts_object_destroy (GTS_OBJECT (status));
  gts_object_destroy (GTS_OBJECT (save_prev));

}


static void generate_surface(GfsSurface *surface, GfsSurface *shape, gdouble rad, Particle *p)
{

  memcpy(surface, shape, sizeof(GfsSurface));

  surface->translate[0] = p->pos.x;
  surface->translate[1] = p->pos.y;
  surface->translate[2] = p->pos.z;

  GtsMatrix * m = gts_matrix_translate (NULL, surface->translate);

  surface->scale[0] = rad;
  surface->scale[1] = rad;
  surface->scale[2] = rad;

  GtsMatrix * ms = gts_matrix_scale (NULL, surface->scale);
  if (surface->m)
    gts_matrix_destroy (surface->m);
  surface->m = gts_matrix_product (m, ms);
  gts_matrix_destroy (m);
  gts_matrix_destroy (ms);
      
  if (surface->s) {
    gts_surface_foreach_vertex (surface->s, (GtsFunc) gts_point_transform, surface->m);
    gts_matrix_destroy (surface->m);
    surface->m = NULL;
    if (surface->flip)
      gts_surface_foreach_face (surface->s, (GtsFunc) gts_triangle_revert, NULL);
  }
  else {
    GtsMatrix * i = gts_matrix_inverse (surface->m);
    gts_matrix_destroy (surface->m);
    surface->m = i;
  }
}

static double compute_xyz (FttCell * cell, GfsDropletToParticles * d)
{
  if(d->fconvert)
    return gfs_function_value (d->fconvert, cell);
  else
    return -1.0;
}

#define PROXIMITY_FRAC 0.25

static void check_cut_cells(FttCell *cell, GfsGenericSurface *s, gpointer *data)
{
  GfsVariable *v = (GfsVariable *)data[0];
  GfsVariable **u = (GfsVariable **)data[1];
  gboolean * check = (gboolean *)data[2];
  Particle * p = (Particle *)data[3];  
  
  if(!(*check)){
    if(GFS_VALUE(cell, v) > PROXIMITY_FRAC){
      *check = TRUE;
    }
  } 
}

static gboolean check_proximity(GfsVariable *v, Particle *p, GfsSurface *surface)
{ 

  GfsDomain *domain = GFS_DOMAIN(gfs_object_simulation(v));

  gboolean check = FALSE;
  gpointer data[4];
  data[0] = v;
  data[1] = gfs_domain_velocity(domain);
  data[2] = &check;
  data[3] = p;

  gfs_domain_traverse_cut (domain,
			   surface,
			   FTT_PRE_ORDER,
			   FTT_TRAVERSE_LEAFS,
			   (FttCellTraverseCutFunc) check_cut_cells,
			   data);
  return check;
}

void refine_implicit_p_cell (FttCell * cell, RefineCut * p)
{
  guint maxlevel = p->maxlevel;

  if (ftt_cell_level (cell) < maxlevel && (gfs_cell_is_cut (cell, p->surface, FALSE, maxlevel))){
    ftt_cell_refine_single (cell, p->domain->cell_init, p->domain->cell_init_data);
    p->check = TRUE;
  }
}


static void convert_particles( GfsDomain *domain, GfsDropletToParticles *d)
{

  GfsLagrangianParticles *lagrangian = LAGRANGIAN_PARTICLES(d);
  
  GSList *i = lagrangian->particles;
  Particle *p;

  while(i){
    p = (Particle *) (i->data);
    i = i->next;

    gdouble rad = pow(p->volume/M_PI,1./2.);
#if !FTT_2D
    rad = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
#endif
    GfsSurface * surface = GFS_GENERIC_SURFACE (gts_object_new 
						  (GTS_OBJECT_CLASS (gfs_surface_class ())));
    generate_surface(surface, d->shape, rad, p);

    gboolean check = FALSE;


    if(rad > ftt_cell_size(p->cell)){
      check = check_proximity(d->c, p, surface);
      if(check)
	g_warning("Transforming Particle into Droplet-proximity based\n");
    }
    else{
      if(GFS_VALUE(p->cell, d->c) > PROXIMITY_FRAC)
	check = TRUE;
      
      if(check)
	g_warning("Transforming Particle into Droplet-proximity based\n");
    }
 
    if(compute_xyz(p->cell, d) > 0 || check){

      RefineCut prefine;
      prefine.domain = domain;
      prefine.maxlevel = d->maxlevel;
      prefine.surface = surface;
      prefine.check = TRUE;

      while(prefine.check){
      	prefine.check = FALSE;
	gfs_domain_cell_traverse (domain,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) refine_implicit_p_cell, &prefine);
      }

      gfs_domain_assign_fraction (domain,
				  surface,
				  d->c, d->resetwith, p);

      lagrangian->particles = g_slist_remove(lagrangian->particles, p);
      g_free(p);
    }

    if (surface->s)
      gts_object_destroy (GTS_OBJECT (surface->s));

    if (surface->m)
      gts_matrix_destroy (surface->m);
  }  

  gfs_domain_reshape (domain, gfs_domain_depth(domain));

}

static gboolean droplet_to_particles_event (GfsEvent * event, GfsSimulation * sim)
{


  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class)->event) 
      (event, sim)) {

    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsLagrangianParticles * lagrangian = LAGRANGIAN_PARTICLES(event);
    GfsDropletToParticles *d = DROPLET_TO_PARTICLES(event);

    d->v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
 
    DropletsPar p ;

    p.resetval = d->resetwith;
    p.tag = gfs_temporary_variable (domain);
    p.u = gfs_domain_velocity(domain);
    p.density = d->density;
    p.t = d->c;

    if (d->v){
      p.c = d->v;
      tag_sort_droplets (domain, &p, d->min);

      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;
	convert_droplets (domain, &p, lagrangian);
      }
    }
    else {

      d->v = gfs_temporary_variable (domain);

      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_v, d);

      p.c = d->v;
 
      tag_sort_droplets (domain, &p, d->min);

       if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;

	convert_droplets (domain, &p, lagrangian);

      }     
 
       gts_object_destroy (GTS_OBJECT (d->v));

    }

    if(d->convert != 0 )
      convert_particles(domain, d);
    
#ifdef HAVE_MPI
    mpi_particle_numbering(domain, lagrangian);
#endif  /*HAVE_MPI*/
    gts_object_destroy(GTS_OBJECT(p.tag));

    if(p.n > 0){
      g_free(p.drops);
      g_free(p.sizes);	
    }

    return TRUE;
  }
  return FALSE;
}

static void droplet_to_particles_class_init (GfsDropletToParticlesClass * klass)
{
  /* define new methods and overload inherited methods here */
  
  GFS_EVENT_CLASS (klass)->event = droplet_to_particles_event;
  GTS_OBJECT_CLASS (klass)->read = droplet_to_particles_read;
  GTS_OBJECT_CLASS (klass)->write = droplet_to_particles_write;
  GTS_OBJECT_CLASS (klass)->destroy = droplet_to_particles_destroy;
}

static void droplet_to_particles_init (GfsDropletToParticles * r)
{
  r->resetwith = 0.;
  r->convert = 0;
  r->maxlevel = 8;

  r->min = 0;
}

GfsDropletToParticlesClass * droplet_to_particles_class (void)
{
  static GfsDropletToParticlesClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo droplet_to_particles_info = {
      "GfsDropletToParticles",
      sizeof (GfsDropletToParticles),
      sizeof (GfsDropletToParticlesClass),
      (GtsObjectClassInitFunc) droplet_to_particles_class_init,
      (GtsObjectInitFunc) droplet_to_particles_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (lagrangian_particles_class ()),
				  &droplet_to_particles_info);
  }

  return klass;
}

