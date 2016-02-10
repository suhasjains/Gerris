/* convert_particles.c
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

#include "lagrangian.h"
#include "refine.h"
#include "adaptive.h"

static void convert_particles_destroy (GtsObject * o)
{
  (* GTS_OBJECT_CLASS (convert_particles_class ())->parent_class->destroy) (o);
  GfsConvertParticles * d = CONVERT_PARTICLES(o);
  
  if (d->fconvert)
    gts_object_destroy (GTS_OBJECT (d->fconvert));

  if (d->shape->s)
    gts_object_destroy (GTS_OBJECT (d->shape->s));

  if (d->shape->m)
    gts_matrix_destroy (d->shape->m);

  if (d->shape)
  g_free (d->shape);

}

static void convert_particles_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (convert_particles_class ())->parent_class->read) (o, fp);

  if (fp->type == GTS_ERROR)
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsConvertParticles * r = CONVERT_PARTICLES(*o);

   while(fp->type == '\n')
     gts_file_next_token (fp);


  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }

  gts_file_next_token (fp);

  if (fp->type == GTS_FLOAT) {
    r->resetwith = atof (fp->token->str);
    gts_file_next_token (fp);
  }

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


static double compute_xyz (FttCell * cell, GfsConvertParticles * d)
{
  if(cell)
    return gfs_function_value (d->fconvert, cell);
  else
    return 1.;
}

static void convert_particles_write (GtsObject * o, FILE * fp)
{

  (* GTS_OBJECT_CLASS (convert_particles_class ())->parent_class->write) (o, fp);

  GfsConvertParticles * r = CONVERT_PARTICLES(o);

  fprintf (fp, " %s %lf %d ", r->c->name, r->resetwith, r->maxlevel);

  if (r->shape) {
    fputs (" ( ", fp);
    gfs_function_write (r->shape->f, fp);
    fputs (" ) ", fp);
  }

  if (r->fconvert){
    gfs_function_write (r->fconvert, fp);
  }

}

/*Reduce the force in the domain due to the particle conversion*/
/*Smoothening the force using Gaussian smoothing Exp(-x^2/sigma^2)*/
static void compute_smooth_force(FttCell *cell, gpointer * data)
{
  if(!cell) return;

  Particle *p = (Particle *)(data[0]);
  GfsVariable **f = data[1];
  gdouble sigma = *(gdouble *)data[2];
  FttVector pos;
  ftt_cell_pos(cell, &pos);


  gdouble dist = (p->pos.x - pos.x)*(p->pos.x - pos.x)
    + (p->pos.y - pos.y)*(p->pos.y - pos.y);
#if !FTT_2D
  dist += (p->pos.z - pos.z)*(p->pos.z - pos.z);
#endif
  dist = exp(-dist/(sigma*sigma))/(2.*M_PI*sigma*sigma);
#if !FTT_2D
  dist /= (pow(2.*M_PI,0.5)*sigma);
#endif

  GFS_VARIABLE(cell, f[0]->i) -= p->phiforce.x*dist;
  GFS_VARIABLE(cell, f[1]->i) -= p->phiforce.y*dist;
#if !FTT_2D
  GFS_VARIABLE(cell, f[2]->i) -= p->phiforce.z*dist;
#endif
}

static gboolean check_stencil(FttCell * cell, FttVector pos0, gdouble sigma, FttDirection *d)
{
  if(!cell) return FALSE;

  gdouble size = ftt_cell_size(cell);
  FttVector pos1;
  ftt_cell_pos(cell, &pos1);

  FttDirection d0;
  gdouble dist, dist1;
  gboolean check = TRUE;
  for(d0 = 0; d0 < FTT_DIMENSION; d0++){
    dist = abs(((&pos1.x)[d0] - size/2.0 - (&pos0.x)[d0]));
    dist1 = abs(((&pos1.x)[d0] + size/2.0 - (&pos0.x)[d0]));

    if(dist < dist1)
      d[d0] = 2.0*d0;
    else
      d[d0] = 2.0*d0 + 1;

    if(dist > sigma || dist1 > sigma)
      check = FALSE;
  }

  return check;
}

/*Gaussian Smoothed Two-way Coupling Force: Applied only on 3Sigma surroundings*/
static void compute_coupling_force (Particle *p, GfsVariable **f)
{
  if(!p->cell) return;

  gpointer data[3];
  data[0] = p;
  data[1] = f;

  gdouble size = ftt_cell_size(p->cell);
  gdouble radius = pow(p->volume/M_PI,1./2.);
#if !FTT_2D
  radius = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
#endif

  gdouble sigma = MAX(radius, size);
  data[2] = &sigma;

  FttCell *cell = p->cell, *neighbor1, *neighbor2;
  FttDirection d[FTT_DIMENSION];
  while(!FTT_CELL_IS_ROOT(cell) && ftt_cell_parent(cell) 
	&& check_stencil(cell, p->pos, 3.*sigma, d)){
    cell = ftt_cell_parent(cell);
  }

/*   if(FTT_CELL_IS_ROOT(cell)) */
/*     return; */

  g_assert(cell!=NULL);


  ftt_cell_traverse (cell, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) compute_smooth_force, data);

  FttDirection d0, d1;
  for(d0 = 0; d0 < FTT_DIMENSION; d0++){
    neighbor1 = ftt_cell_neighbor(cell, d[d0]);
    if(neighbor1)
      ftt_cell_traverse (neighbor1, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			 (FttCellTraverseFunc) compute_smooth_force, data);
    else
      continue;
  /*   Do for neighbor if it exists else continue */
    for(d1 = d0 + 1; d1 < FTT_DIMENSION; d1++){
      neighbor2 = ftt_cell_neighbor(neighbor1, d[d1]);
      if(neighbor2){
	ftt_cell_traverse (neighbor2, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			   (FttCellTraverseFunc) compute_smooth_force, data);
      }
   /*    Do for neighbor if it exists else continue */
    }
  }

#if !FTT_2D
  /*Corner Cell*/
  neighbor1 = ftt_cell_neighbor(cell, d[0]);
  if(neighbor1)
    neighbor1 = ftt_cell_neighbor(neighbor1, d[1]);
  if(neighbor1)
    neighbor1 = ftt_cell_neighbor(neighbor1, d[2]);
  if(neighbor1)
    ftt_cell_traverse (neighbor1, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		       (FttCellTraverseFunc) compute_smooth_force, data);
#endif

}
/******************************************/
static void convert_particles( GfsDomain *domain, GfsConvertParticles *d)
{

  GfsLagrangianParticles *lagrangian = LAGRANGIAN_PARTICLES(d);
  
  GSList *i = lagrangian->particles;
  Particle *p;

  while(i){
    p = (Particle *) (i->data);
    i = i->next;

    if(compute_xyz(p->cell, d) != 0){
 
      gdouble rad = pow(p->volume/M_PI,1./2.);
#if !FTT_2D
      rad = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
#endif
      GfsSurface * surface = GFS_GENERIC_SURFACE (gts_object_new 
						  (GTS_OBJECT_CLASS (gfs_surface_class ())));

      memcpy(surface, d->shape, sizeof(GfsSurface));


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

      RefineCut prefine;
      prefine.domain = domain;
      prefine.maxlevel = d->maxlevel;
      prefine.surface = surface;
      prefine.check = TRUE;

      p->cell = gfs_domain_locate(domain, p->pos, -1);
      compute_coupling_force (p, lagrangian->couplingforce);

      while(prefine.check){
      	prefine.check = FALSE;
	gfs_domain_cell_traverse (domain,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) refine_implicit_p_cell, &prefine);
      }

      gfs_domain_assign_fraction (domain,
				  surface,
				  d->c, (1.- d->resetwith), p); 


      lagrangian->particles = g_slist_remove(lagrangian->particles, p);
      g_free(p);

      if (surface->s)
	gts_object_destroy (GTS_OBJECT (surface->s));

      if (surface->m)
	gts_matrix_destroy (surface->m);
    }
  }  
 /*  gfs_domain_reshape (domain, d->maxlevel); */
  gfs_domain_reshape (domain, gfs_domain_depth(domain));
}

static gboolean convert_particles_event (GfsEvent * event, GfsSimulation * sim)
{

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (convert_particles_class ())->parent_class)->event) 
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsConvertParticles *d = CONVERT_PARTICLES(event);

    convert_particles(domain, d);
    return TRUE;
  }
  return FALSE;
}

static void convert_particles_class_init (GfsConvertParticlesClass * klass)
{
  /* define new methods and overload inherited methods here */

  GFS_EVENT_CLASS (klass)->event = convert_particles_event;
  GTS_OBJECT_CLASS (klass)->read = convert_particles_read;
  GTS_OBJECT_CLASS (klass)->write = convert_particles_write;
  GTS_OBJECT_CLASS (klass)->destroy = convert_particles_destroy;
}

static void convert_particles_init (GfsConvertParticles * r)
{
  r->resetwith = 0.;
  r->maxlevel = 8;
}

GfsConvertParticlesClass * convert_particles_class (void)
{
  static GfsConvertParticlesClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo convert_particles_info = {
      "GfsConvertParticles",
      sizeof (GfsConvertParticles),
      sizeof (GfsConvertParticlesClass),
      (GtsObjectClassInitFunc) convert_particles_class_init,
      (GtsObjectInitFunc) convert_particles_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (feed_particles_class ()),
				  &convert_particles_info);
  }

  return klass;
}

/* Initialize modules */

/* const gchar * g_module_check_init (void); */

/* const gchar * g_module_check_init (void) */
/* { */
/*   convert_particles_class (); */
/*   return NULL; */
/* }  */
