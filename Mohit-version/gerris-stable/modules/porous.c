#include <math.h>
#include <stdlib.h>

#include "simulation.h"
#include "timestep.c"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"
#include "mpi_boundary.h"
#include "event.h"
#include "ftt.h"
#include "domain.h"
#include "gfs.h"
#include "refine.h"




/* GfsPorous: Header */

typedef struct _GfsPorous         GfsPorous;
typedef struct _GfsPorousClass    GfsPorousClass;

struct _GfsPorous {
  /*< private >*/
  GfsSimulation parent;

  /*< public >*/
  /* add extra data here (if public) */
  GfsFunction * porosity ;                          /* porosity given by epsilon*/
  GfsFunction * K ;                                 /* permeability of medium */
};

struct _GfsPorousClass {
  /*< private >*/
  GfsSimulationClass parent_class;

  /*< public >*/
  /* add extra methods here */
};
typedef struct {
  gdouble lambda2[FTT_DIMENSION];
  GfsFunction * alpha, * phi;
  GfsDomain * domain;
  gboolean positive;
} PoissonCoeff_por;

#define GFS_POROUS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsPorous,\
					         gfs_porous_class ())
#define GFS_POROUS_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsPorousClass,\
						 gfs_porous_class())
#define IS_GFS_POROUS(obj)         (gts_object_is_from_class (obj,\
						 gfs_porous_class ()))

static void gfs_porous_class_init (GfsSimulationClass * klass);
static void gfs_porous_read (GtsObject ** o, GtsFile * fp);
static void gfs_porous_write (GtsObject * o, FILE * fp);
static void gfs_porous_run (GfsSimulation * sim);
GfsPorousClass * gfs_porous_class  (void);
GfsPorous * gfs_porous_new    (GfsPorousClass * klass);


void gfs_pre_projection (GfsDomain *domain, GfsPorous *por, guint dimension);
void gfs_post_projection (GfsDomain *domain, GfsPorous *por, guint dimension);

static void face_coeff_from_below_por (FttCell * cell);
static void correct_por (FttCell *cell, gpointer *data);
static void correct_por_undo (FttCell *cell, gpointer *data);



void gfs_approximate_projection_por (GfsDomain * domain,
                                     GfsPorous * por,
                                     GfsMultilevelParams * par,
                                     gdouble dt,
                                     GfsVariable * p,
                                     GfsFunction * alpha,
                                     GfsVariable * res,
                                     GfsVariable ** g,
                                     void (* divergence_hook) (GfsDomain * domain,
                                                           gdouble dt,
                                                           GfsVariable * div)
                                     );


static void mac_projection_por (GfsDomain * domain,
                                GfsPorous * por,
                                GfsMultilevelParams * par,
                                gdouble dt,
                                GfsVariable * p,
                                GfsFunction * alpha,
                                GfsVariable * res,
                                GfsVariable ** g,
                                void (* divergence_hook) (GfsDomain * domain,
                                                      gdouble dt,
                                                      GfsVariable * div)
                            );


void gfs_poisson_coefficients_por (GfsDomain * domain,
                                   GfsPorous *por,
                                   GfsFunction * alpha,
                                   gboolean positive,
                                   gboolean centered,
                                   gboolean reset);

static void poisson_mixed_coeff_por (FttCell * cell, PoissonCoeff_por * p);

static void poisson_coeff_por (FttCellFace * face,
                           PoissonCoeff_por * p);


static void reset_coeff_por (FttCell * cell, PoissonCoeff_por * p);


void gfs_mac_projection_por (GfsDomain * domain,
                             GfsPorous * por,
                             GfsMultilevelParams * par,
                             gdouble dt,
                             GfsVariable * p,
                             GfsFunction * alpha,
                             GfsVariable ** g,
                             void (* divergence_hook) (GfsDomain * domain,
                                                   gdouble dt,
                                                   GfsVariable * div)
                             );


void gfs_update_gradients_por (GfsDomain * domain,
                               GfsPorous * por,
                               GfsVariable * p, 
                               GfsFunction * alpha,
                               GfsVariable ** g);


/*               GfsSourceDarcy Header               */

typedef struct _GfsSourceDarcy         GfsSourceDarcy;

struct _GfsSourceDarcy {
	/*< private >*/
	GfsSourceVelocity parent;
	GfsVariable * u[FTT_DIMENSION];
	/*< public >*/
	GfsFunction * darcycoeff; /* Linear Drag coefficient */
	GfsFunction * forchhicoeff; /* Non-linear Forchheimer coefficient */
	gdouble         beta; /* "implicitness coefficient" 0.5 CN, 1 backward Euler */
	/* add extra data here (if public) */
};
static gboolean gfs_source_darcy_event (GfsEvent * event, GfsSimulation * sim);
static void gfs_source_darcy_class_init (GfsSourceGenericClass * klass);
static void gfs_source_darcy_init (GfsSourceGeneric * s);
static gdouble gfs_source_darcy_centered_value (GfsSourceGeneric * s,
						   FttCell * cell,
						   GfsVariable * v);

static void implicit_darcy_2D (FttCell * cell, GfsSourceDarcy * s);



/*Definitions of functions modified for Porous Media */

static void face_coeff_from_below_por (FttCell * cell)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  guint neighbors = 0;

  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellChildren child;
    guint i, n;

    f[d].v = 0.;
    n = ftt_cell_children_direction (cell, d, &child);
    for (i = 0; i < n; i++)
      if (child.c[i])
        f[d].v += GFS_STATE (child.c[i])->f[d].v;
    f[d].v /= n;
    /* fixme: this stuff may not be necessary anymore? The 'dumbell'
       test case seems to work fine without this */
    FttCell * neighbor;
    if (f[d].v != 0. &&
        (neighbor = ftt_cell_neighbor (cell, d)) && !GFS_CELL_IS_BOUNDARY (neighbor))
      neighbors++;
  }

  if (neighbors == 1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      f[d].v = 0.;
}


static void correct_por (FttCell * cell, gpointer * data)
{
  FttComponent c;
  GfsVariable **v = data[0];
  guint * dimension = data[1];
  GfsPorous *por = data[2];
  
  for (c = 0; c < *dimension; c++)
  GFS_VALUE (cell, v[c]) *= gfs_function_value (por->porosity, cell); 
}

static void correct_por_undo (FttCell * cell, gpointer * data)
{
  FttComponent c;
  GfsVariable **v = data[0];
  guint * dimension = data[1];
  GfsPorous *por = data[2];

  for (c = 0; c < *dimension; c++)
  GFS_VALUE (cell, v[c]) *= (1/gfs_function_value (por->porosity, cell)); 
}

void gfs_pre_projection(GfsDomain *domain, GfsPorous *por, guint dimension)
{
GfsVariable **v;
FttComponent c;
gpointer data[3];

g_return_if_fail (domain != NULL);
data[0] = v = gfs_domain_velocity (domain);
data[1] = &dimension;
data[2] = por;
    
gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, (FttCellTraverseFunc) correct_por, data);

for (c = 0; c < dimension; c++)
gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v[c]);
printf("\nvelocity correction working!!!\n");

}

void gfs_post_projection (GfsDomain *domain, GfsPorous *por, guint dimension)
{
GfsVariable **v;
FttComponent c;
gpointer data[3];

g_return_if_fail (domain != NULL);
data[0] = v = gfs_domain_velocity (domain);
data[1] = &dimension;
data[2] = por;
    
gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, (FttCellTraverseFunc) correct_por_undo, data);

for (c = 0; c < dimension; c++)
gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v[c]);
printf("\nvelocity correction undo !!!\n");

}


static void reset_coeff_por (FttCell * cell, PoissonCoeff_por * p)
{
  FttDirection d;
  GfsFaceStateVector * f = GFS_STATE (cell)->f;
  if (GFS_IS_MIXED (cell)) {
    FttVector v = {0.,0.,0.};
    GFS_STATE (cell)->solid->v = v;
  }
  for (d = 0; d < FTT_NEIGHBORS; d++)
    f[d].v = 0.;
}

static void poisson_coeff_por (FttCellFace * face,
                           PoissonCoeff_por * p)
{
  gdouble alpha = p->alpha ? (gfs_function_face_value (p->alpha, face)*gfs_function_face_value (p->phi, face)) : gfs_function_face_value (p->phi, face);
  gdouble v = p->lambda2[face->d/2]*alpha*gfs_domain_face_fraction (p->domain, face)/
    gfs_domain_face_scale_metric (p->domain, face, face->d/2);

  if (alpha <= 0. && p->positive) {
    FttVector p;
    ftt_face_pos (face, &p);
    g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
           "alpha is negative (%g) at face (%g,%g,%g).\n"
           "Please check your definition.",
           alpha, p.x, p.y, p.z);
  }
  GFS_STATE (face->cell)->f[face->d].v += v;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v += v;
    break;
  case FTT_FINE_COARSE:
    GFS_STATE (face->neighbor)->f[FTT_OPPOSITE_DIRECTION (face->d)].v +=
      v/FTT_CELLS_DIRECTION (face->d);
    break;
  default:
    g_assert_not_reached ();
  }
}

static void poisson_mixed_coeff_por (FttCell * cell, PoissonCoeff_por * p)
{
  if (GFS_IS_MIXED (cell)) {
    gdouble alpha = p->alpha ? (gfs_function_value (p->alpha, cell)*gfs_function_value (p->phi , cell)) : gfs_function_value (p->phi, cell);
    if (((cell)->flags & GFS_FLAG_DIRICHLET) == 0)
      /* Neumann condition (prescribed flux) */
      GFS_STATE (cell)->solid->v.x += alpha;
    else {
      /* Dirichlet */
      GfsSolidVector * s = GFS_STATE (cell)->solid;
      FttVector m = {1.,1.,1.};
      gfs_domain_solid_metric (p->domain, cell, &m);
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
        (&s->v.x)[c] += alpha*(&m.x)[c]*(s->s[2*c + 1] - s->s[2*c]);
    }

    if (alpha <= 0. && p->positive) {
      FttVector p;
      ftt_cell_pos (cell, &p);
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
             "alpha is negative (%g) at cell (%g,%g,%g).\n"
             "Please check your definition.",
             alpha, p.x, p.y, p.z);
    }
  }
}

void gfs_poisson_coefficients_por (GfsDomain * domain,
                                   GfsPorous *por,
                                   GfsFunction * alpha,
                                   gboolean positive,
                                   gboolean centered,
                                   gboolean reset)
{
  PoissonCoeff_por p;
  FttComponent i;

  g_return_if_fail (domain != NULL);

  for (i = 0; i < FTT_DIMENSION; i++) {
    gdouble lambda = (&domain->lambda.x)[i];

    p.lambda2[i] = lambda*lambda;
  }
  p.alpha = alpha;
  p.domain = domain;
  p.positive = positive;
  p.phi = por->porosity;

  if (reset)
    gfs_domain_cell_traverse (domain,
                              FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                              (FttCellTraverseFunc) reset_coeff_por, &p);
  if (!centered)
    gfs_domain_cell_traverse (domain,
                              FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                              (FttCellTraverseFunc) poisson_mixed_coeff_por, &p);
  gfs_domain_face_traverse (domain, FTT_XYZ,
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttFaceTraverseFunc) poisson_coeff_por, &p);
  gfs_domain_cell_traverse (domain,
                            FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
                            (FttCellTraverseFunc) face_coeff_from_below_por, NULL);
}

static void mac_projection_por (GfsDomain * domain,
                                GfsPorous * por,
                                GfsMultilevelParams * par,
                                gdouble dt,
                                GfsVariable * p,
                                GfsFunction * alpha,
                                GfsVariable * res,
                                GfsVariable ** g,
                                void (* divergence_hook) (GfsDomain * domain,
                                                      gdouble dt,
                                                      GfsVariable * div)
                            )
{
  /* Add face sources */
  gfs_reset_gradients (domain, FTT_DIMENSION, g);
  gfs_velocity_face_sources (domain, gfs_domain_velocity (domain), dt, alpha, g);
    
  GfsVariable * dia = gfs_temporary_variable (domain);
  GfsVariable * div = gfs_temporary_variable (domain);
  GfsVariable * res1 = res ? res : gfs_temporary_variable (domain);
  
  /* Initialize face coefficients */
  gfs_poisson_coefficients_por (domain, por, alpha, TRUE, TRUE, TRUE);
  
  /* hydrostatic pressure */ 
  GSList * i = domain->variables;
  while (i) {
    if (GFS_IS_HYDROSTATIC_PRESSURE (i->data))
      gfs_correct_normal_velocities (domain, FTT_DIMENSION, i->data, g, dt);
    i = i->next;
  }
  
  /* Initialize diagonal coefficient */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                            (FttCellTraverseFunc) gfs_cell_reset, dia);

  /* compute MAC divergence */
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) gfs_normal_divergence, div);

  /* Divergence hook */
  if (divergence_hook)
    (* divergence_hook) (domain, dt, div);

  /* add volume sources (if any) */
  if (p->sources)
    volume_sources (domain, p, div);

  /* Scale divergence */
  gpointer data[2];
  data[0] = div;
  data[1] = &dt;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) scale_divergence, data);

#if 0
  {
    FILE * fp = fopen ("/tmp/mac", "wt");
    GfsNorm norm;

    gfs_write_mac_velocity (domain, 0.9, FTT_TRAVERSE_LEAFS, -1, NULL, fp);
    fclose (fp);
    norm = gfs_domain_norm_variable (domain, div, FTT_TRAVERSE_LEAFS, -1);
    fprintf (stderr, "mac div before: %g %g %g\n",
             norm.first, norm.second, norm.infty);
  }
#endif

par->poisson_solve (domain, par, p, div, res1, dia, dt);

  gts_object_destroy (GTS_OBJECT (dia));
  gts_object_destroy (GTS_OBJECT (div));
  if (!res)
    gts_object_destroy (GTS_OBJECT (res1));

  gfs_correct_normal_velocities (domain, FTT_DIMENSION, p, g, dt);
  gfs_scale_gradients (domain, FTT_DIMENSION, g);
}

void gfs_approximate_projection_por (GfsDomain * domain,
                                     GfsPorous * por,
                                     GfsMultilevelParams * par,
                                     gdouble dt,
                                     GfsVariable * p,
                                     GfsFunction * alpha,
                                     GfsVariable * res,
                                     GfsVariable ** g,
                                     void (* divergence_hook) (GfsDomain * domain, 
                                                           gdouble dt,
                                                           GfsVariable * div)
                                     )
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);

  gfs_domain_timer_start (domain, "approximate_projection");
  
  gfs_pre_projection (domain, por, FTT_DIMENSION);

  /* compute MAC velocities from centered velocities */
  gfs_domain_face_traverse (domain, FTT_XYZ,
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
  gfs_domain_face_traverse (domain, FTT_XYZ,
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity, 
                            gfs_domain_velocity (domain));
  
  mac_projection_por (domain, por, par, dt, p, alpha, res, g, divergence_hook);
  
  gfs_correct_centered_velocities (domain, FTT_DIMENSION, g, dt);

  gfs_post_projection (domain, por, FTT_DIMENSION);

  gfs_domain_timer_stop (domain, "approximate_projection");

  if (par->residual.infty > par->tolerance)
    g_warning ("approx projection: max residual %g > %g", par->residual.infty, par->tolerance);
}

void gfs_mac_projection_por (GfsDomain * domain,
                             GfsPorous * por,
                             GfsMultilevelParams * par,
                             gdouble dt,
                             GfsVariable * p,
                             GfsFunction * alpha,
                             GfsVariable ** g,
                             void (* divergence_hook) (GfsDomain * domain,
                                                   gdouble dt,
                                                   GfsVariable * div)
                             ) 
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);

  gfs_domain_timer_start (domain, "mac_projection");

  mac_projection_por (domain, por, par, dt, p, alpha, NULL, g, divergence_hook); 

  gfs_domain_timer_stop (domain, "mac_projection");

  if (par->residual.infty > par->tolerance)
    g_warning ("MAC projection: max residual %g > %g", par->residual.infty, par->tolerance);
}


void gfs_update_gradients_por (GfsDomain * domain,
                               GfsPorous * por,
                               GfsVariable * p,  
                               GfsFunction * alpha,
                               GfsVariable ** g)
{        
  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);
  g_return_if_fail (g != NULL);
  
  /* Add face sources */
  gfs_reset_gradients (domain, FTT_DIMENSION, g);
  gfs_velocity_face_sources (domain, gfs_domain_velocity (domain), 0., alpha, g);
  /* Initialize face coefficients */
  gfs_poisson_coefficients_por (domain, por, alpha, TRUE, TRUE, TRUE);
  /* Add pressure gradient */
  gfs_correct_normal_velocities (domain, FTT_DIMENSION, p, g, 0.); 
  gfs_scale_gradients (domain, FTT_DIMENSION, g);
}                               



/* GfsSourceDarcy: Header */


GfsSourceGenericClass * gfs_source_darcy_class(void);

#define GFS_SOURCE_DARCY(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceDarcy,\
					         gfs_source_darcy_class ())
#define GFS_IS_SOURCE_DARCY(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_darcy_class ()))
/*Header for Source Darcy Funcs*/
void gfs_source_darcy_implicit (GfsDomain * domain, gdouble dt);

GfsSourceDarcy * gfs_has_source_darcy (GfsDomain * domain);





/* GfsPorous: Object */
static void gfs_porous_class_init (GfsSimulationClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = gfs_porous_read;
  printf("\npass-1\n");
  GTS_OBJECT_CLASS (klass)->write = gfs_porous_write;
  printf("\npass-2\n");
  klass->run = gfs_porous_run ;
}

static void gfs_porous_init (GfsPorous * object)
{
  /* initialize object here */
  object->porosity = gfs_function_new (gfs_function_class (), 0.);
  object->K = gfs_function_new (gfs_function_class (), 0.);
}

GfsPorousClass * gfs_porous_class (void)
{
  static GfsPorousClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_porous_info = {
      "GfsPorous",
      sizeof (GfsPorous),
      sizeof (GfsPorousClass),
      (GtsObjectClassInitFunc) gfs_porous_class_init,
      (GtsObjectInitFunc) gfs_porous_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()),
				  &gfs_porous_info);
  }

  return klass;
}

static void gfs_porous_read (GtsObject ** o, GtsFile * fp)
{


	/* call read method of parent */
	(* GTS_OBJECT_CLASS (gfs_porous_class ())->parent_class->read) (o, fp);
	if (fp->type == GTS_ERROR)
		return;

	GfsPorous * por = GFS_POROUS (*o);
	GfsSimulation * sim = GFS_SIMULATION (por);
	/* do object-specific read here */

	if (fp->type != '{') {
		gts_file_error (fp, "expecting an opening brace");
		return;
	}

	fp->scope_max++;
	gts_file_next_token (fp);

	while (fp->type != GTS_ERROR && fp->type != '}') {
		if (fp->type == '\n') {
			gts_file_next_token (fp);
			continue;
		}

		if (fp->type != GTS_STRING) {
			gts_file_error (fp, "expecting a keyword");
			return;
		}

		if (!strcmp (fp->token->str, "porosity")) {
			gts_file_next_token (fp);
			if (fp->type != '=')
				gts_file_error (fp, "expecting `='");
			else{
				gts_file_next_token (fp);
				gfs_function_read (por->porosity, sim, fp);
			}
		}


		if (!strcmp (fp->token->str, "K")) {
			gts_file_next_token (fp);
			if (fp->type != '=')
				gts_file_error (fp, "expecting `='");
			else{
				gts_file_next_token (fp);
				gfs_function_read (por-> K, sim , fp);
			}
		}
	}
	if (fp->type == GTS_ERROR)
		return;
	if (fp->type != '}') {
		gts_file_error (fp, "expecting a closing brace");
		return;

		fp->scope_max--;
		gts_file_next_token (fp);
		/* do not forget to prepare for next read */
	}
}

static void gfs_porous_write (GtsObject * o, FILE * fp)
{
  /* call write method of parent */
    (* GTS_OBJECT_CLASS (gfs_porous_class ())->parent_class->write)(o, fp);
  /* do object specific write here */
  GfsPorous * por = GFS_POROUS (o);
  
  fputs (" {\n"
	 "  porosity =", fp);
  gfs_function_write (por->porosity, fp);
  fputs ("\n  K =", fp);
  gfs_function_write (por->K, fp);
  fputs("}\n ", fp);
}


static void gfs_porous_run (GfsSimulation * sim)
{
  GfsVariable * p, * pmac, * res = NULL, * g[FTT_DIMENSION], * gmac[FTT_DIMENSION];
  GfsVariable ** gc = sim->advection_params.gc ? g : NULL;
  GfsDomain * domain;
  GfsPorous *por;
  GSList * i;
    
  domain = GFS_DOMAIN (sim);
  por = GFS_POROUS (sim);

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gmac[c] = gfs_temporary_variable (domain);
    if (sim->advection_params.gc)
      g[c] = gfs_temporary_variable (domain);
    else
      g[c] = gmac[c];
  }
  gfs_variable_set_vector (gmac, FTT_DIMENSION);
  gfs_variable_set_vector (g, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  gfs_simulation_set_timestep (sim);
  if (sim->time.i == 0) {

    /*inserted changes inside this function*/
    gfs_approximate_projection_por (domain, por,
				&sim->approx_projection_params,
				sim->advection_params.dt,
				p, sim->physical_params.alpha, res, g, NULL);


    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt/2.);
  }
  else if (sim->advection_params.gc)
    gfs_update_gradients_por (domain, por, p, sim->physical_params.alpha, g);


  while (sim->time.t < sim->time.end &&
	 sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);
    
    /*inserted changes */
      gfs_pre_projection (domain, por, FTT_DIMENSION);
      
    if (sim->advection_params.linear) {
      /* linearised advection */

      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
      gfs_domain_face_traverse (domain, FTT_XYZ,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity,
				sim->u0);
    }

    else
      gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);
      
    gfs_variables_swap (p, pmac);


    gfs_mac_projection_por (domain, por,
    			&sim->projection_params, 
    			sim->advection_params.dt/2.,
			p, sim->physical_params.alpha, gmac, NULL);


    gfs_variables_swap (p, pmac);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_centered_velocity_advection_diffusion (domain,
					       FTT_DIMENSION,
					       &sim->advection_params,
					       gmac,
					       sim->time.i > 0 || !gc ? gc : gmac,
					       sim->physical_params.alpha);
    if (gc) {
      gfs_source_darcy_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, sim->time.i > 0 ? gc : gmac, 
				       -sim->advection_params.dt);
      /*inserted changes*/
      gfs_post_projection (domain, por, FTT_DIMENSION);
}
    else if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, sim->advection_params.dt);
      gfs_source_darcy_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, -sim->advection_params.dt);
      /*inserted changes*/
      gfs_post_projection (domain, por, FTT_DIMENSION);
   
 }

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    /*inserted changes */

    gfs_approximate_projection_por (domain, por,
   				&sim->approx_projection_params, 
    				sim->advection_params.dt, 
				p, sim->physical_params.alpha, res, g, NULL);

    /*inserted changes */

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) {
    gts_object_destroy (GTS_OBJECT (gmac[c]));
    if (sim->advection_params.gc)
      gts_object_destroy (GTS_OBJECT (g[c]));
  }
}


/* GfsSourceDarcy: Object */

static void gfs_source_darcy_destroy (GtsObject * o)
{
  /* do object-specific cleanup here */
  GfsSourceDarcy * s = GFS_SOURCE_DARCY (o);
  FttComponent c;

  if (s->darcycoeff)
    gts_object_destroy (GTS_OBJECT (s->darcycoeff));
  if (s->forchhicoeff)
    gts_object_destroy (GTS_OBJECT (s->forchhicoeff));

  for (c = 0; c <  FTT_DIMENSION; c++)
    if (s->u[c])
      gts_object_destroy (GTS_OBJECT (s->u[c]));
  /* do not forget to call destroy method of the parent */
  (* GTS_OBJECT_CLASS (gfs_source_darcy_class ())->parent_class->destroy) 
    (o);
}
	
static void gfs_source_darcy_read (GtsObject ** o, GtsFile * fp)
{
    (* GTS_OBJECT_CLASS (gfs_source_darcy_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  printf("\ntesting1\n");
  /*Check if source darcy has already been added or not*/
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = GFS_SOURCE_VELOCITY (*o)->v[c];

    if (v->sources) {
      GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
      
      while (i) {
	if (i->data != *o && GFS_IS_SOURCE_DARCY (i->data)) {
	  gts_file_error (fp, "variable '%s' cannot have multiple Darcy source terms", v->name);
	  return;
	}
	i = i->next;
      }
    }
  }

  printf("\ntesting2\n");
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSourceDarcy * s = GFS_SOURCE_DARCY (*o);
  printf("\ntesting3\n");
  s->darcycoeff = gfs_function_new (gfs_function_class (), 0.);
  gfs_function_read (s->darcycoeff, gfs_object_simulation (s), fp);
  printf("\ntesting4\n");
  
  if (fp->type != '\n') {
    s->forchhicoeff = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (s->forchhicoeff, gfs_object_simulation (s), fp);
  }

  if (s->beta < 1.)
    for (c = 0; c <  FTT_DIMENSION; c++)
      s->u[c] = gfs_temporary_variable (domain);
  else {
    GFS_SOURCE_GENERIC (s)->centered_value = NULL;
    GFS_SOURCE_GENERIC (s)->mac_value = NULL;
  }
  printf("\ntesting5\n");
}


static void gfs_source_darcy_write (GtsObject * o, FILE * fp)
{
  GfsSourceDarcy * s = GFS_SOURCE_DARCY (o);

  (* GTS_OBJECT_CLASS (gfs_source_darcy_class ())->parent_class->write) (o, fp);
  gfs_function_write (s->darcycoeff, fp);
  if (s->forchhicoeff)
    gfs_function_write (s->forchhicoeff, fp);
}

static GfsSourceDiffusion * source_diffusion_viscosity (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;
    
    while (i) {
      GtsObject * o = i->data;
      
      if (GFS_IS_SOURCE_DIFFUSION (o))
	return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static void darcy_coefficients (GfsSourceDarcy * sd, FttCell * cell, GfsVariable ** u, FttComponent c1, gdouble f[2])
{
    GfsPorous * por = GFS_POROUS (gfs_object_simulation (sd));
    GfsDomain * domain = GFS_DOMAIN(por);
    GfsFunction * porosity = por->porosity;
    GfsFunction * K = por->K;
    gdouble viscosity = 0.001;
    GfsVariable ** U = gfs_domain_velocity (domain);
    GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (domain);
    GfsSourceDiffusion * d = source_diffusion_viscosity (U[0]); 

    if (d)
    viscosity = gfs_diffusion_cell (d->D, cell);
       
    f[0] = (viscosity * gfs_function_value (porosity,cell))/gfs_function_value (K,cell);
    
    gdouble modu;
    modu = fabs(sqrt(GFS_VALUE (cell, sv->v[0])*GFS_VALUE (cell, sv->v[0]) + 
		   GFS_VALUE (cell, sv->v[1])*GFS_VALUE (cell, sv->v[1])));

    f[1] = (gfs_function_value (porosity,cell))/sqrt(gfs_function_value (K,cell))*modu;
    
  
}

static gdouble gfs_source_darcy_mac_value (GfsSourceGeneric * s,
					      FttCell * cell,
					      GfsVariable * v)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
  GfsSourceDarcy * sd = GFS_SOURCE_DARCY (s);

#if FTT_2D
  gdouble f[2];
  darcy_coefficients (sd, cell, sv->v, v->component, f);
  switch (v->component) {
  case FTT_X: return    -(f[0]+f[1])*GFS_VALUE (cell, sv->v[0]);
  case FTT_Y: return    -(f[0]+f[1])*GFS_VALUE (cell, sv->v[1]);
  default: g_assert_not_reached ();
  }
#else  /* 3D */
  gdouble f = gfs_function_value (sd->darcycoeff, cell);

//gdouble e = sd->forchhi ? gfs_function_value (sd->forchhi, cell) : 0.;

  switch (v->component) {
  case FTT_X: return - f*GFS_VALUE (cell, sv->v[0]);
  case FTT_Y: return - f*GFS_VALUE (cell, sv->v[1]);
  case FTT_Z: return - f*GFS_VALUE (cell, sv->v[2]);
  default: g_assert_not_reached ();
  }
#endif /* 3D */
  return 0.;
}

static void save_darcy (FttCell * cell, GfsSourceDarcy * s)
{
  GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);

  gdouble f[2];
  darcy_coefficients (s, cell, sv->v, FTT_X, f);
  darcy_coefficients (s, cell, sv->v, FTT_Y, f);

  gdouble d = s->darcycoeff ? gfs_function_value (s->darcycoeff, cell)*(1. - s->beta):0.;
  gdouble e = s->forchhicoeff ? gfs_function_value (s->forchhicoeff, cell)*(1. - s->beta) : 0.;

  d *= f[0];
  e *= f[1];

  GFS_VALUE (cell, s->u[0]) = - (d+e)*GFS_VALUE (cell, sv->v[0]);
  GFS_VALUE (cell, s->u[1]) = - (d+e)*GFS_VALUE (cell, sv->v[1]);
#if !FTT_2D
  GFS_VALUE (cell, s->u[2]) = - (d+e)*GFS_VALUE (cell, sv->v[2]);
#endif /* 3D */
}


static gboolean gfs_source_darcy_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_darcy_class ())->parent_class)->event) (event, sim)) {
    /* do object-specific event here */
    if (GFS_SOURCE_DARCY (event)->beta < 1.) {
      gfs_catch_floating_point_exceptions ();
      gfs_domain_traverse_layers (GFS_DOMAIN (sim), (FttCellTraverseFunc) save_darcy, event);
      if (gfs_restore_floating_point_exceptions ()) {
	GfsSourceDarcy * c = GFS_SOURCE_DARCY (event);
	gchar * s = g_strconcat ("\n", gfs_function_description (c->darcycoeff, FALSE), NULL);
	if (c->forchhicoeff)
	  s = g_strconcat (s, "\n", gfs_function_description (c->forchhicoeff, FALSE), NULL);
	/* fixme: memory leaks */
	g_message ("floating-point exception in user-defined function(s):%s", s);
	exit (1);
      }
    }			  
   return TRUE;			      
  }
  return FALSE;
}

static gdouble gfs_source_darcy_centered_value (GfsSourceGeneric * s,
						   FttCell * cell,
						   GfsVariable * v)
{
  return GFS_VALUE (cell, GFS_SOURCE_DARCY (s)->u[v->component]);
}

static void gfs_source_darcy_class_init (GfsSourceGenericClass * klass)
{
  /* define new methods and overload inherited methods here */

  GFS_EVENT_CLASS (klass)->event = gfs_source_darcy_event;
  printf("\ndarcy init\n");
  GTS_OBJECT_CLASS (klass)->read = gfs_source_darcy_read;
  printf("\ndarcy read\n");
  GTS_OBJECT_CLASS (klass)->write = gfs_source_darcy_write;
  printf("\ndarcy written\n");
  GTS_OBJECT_CLASS (klass)->destroy = gfs_source_darcy_destroy;
}

static void gfs_source_darcy_init (GfsSourceGeneric * s)
{
  /* initialize object here */
  s->mac_value = gfs_source_darcy_mac_value;
  GFS_SOURCE_DARCY (s)->beta = 0.5; /* Crank-Nicholson */
  s->centered_value = gfs_source_darcy_centered_value;
}

GfsSourceGenericClass * gfs_source_darcy_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_darcy_info = {
      "GfsSourceDarcy",
      sizeof (GfsSourceDarcy),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) gfs_source_darcy_class_init,
      (GtsObjectInitFunc) gfs_source_darcy_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
				  &gfs_source_darcy_info);
  }
  return klass;
}

/**
 * gfs_has_source_darcy:
 * @domain: a #GfsDomain.
 *
 * Returns: the #GfsSourceDarcy associated with @domain or %NULL.
 */

GfsSourceDarcy * gfs_has_source_darcy (GfsDomain * domain)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  v = gfs_variable_from_name (domain->variables, "U");
  g_return_val_if_fail (v != NULL, NULL);

  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    while (i) {
      if (GFS_IS_SOURCE_DARCY (i->data))
	return i->data;
      i = i->next;
    }
  }
  return NULL;
}

static void implicit_darcy_2D (FttCell * cell, GfsSourceDarcy * s)
{
    
    GfsSourceVelocity * sv = GFS_SOURCE_VELOCITY (s);
    gdouble m[2][2];
    gdouble f[2];
    GfsSimulation * sim = gfs_object_simulation (s);
    gdouble dt = sim->advection_params.dt;

    darcy_coefficients (s, cell, sv->v, FTT_X, f);
    gdouble d = s->darcycoeff ? gfs_function_value (s->darcycoeff, cell):0.;
    gdouble e = s->forchhicoeff ? gfs_function_value (s->forchhicoeff, cell):0.;

    f[0] *= d;
    f[1] *= e;

    m[0][0] = 1. + (f[0]+f[1])*dt;
    m[0][1] = 0;
    darcy_coefficients (s, cell, sv->v, FTT_Y, f);
    m[1][0] = 0;
    m[1][1] = 1. + (f[0]+f[1])*dt;

    gdouble det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
    gdouble u = GFS_VALUE (cell, sv->v[0]);
    gdouble v = GFS_VALUE (cell, sv->v[1]);

    GFS_VALUE (cell, sv->v[0]) = (  m[1][1]*u - m[0][1]*v)/det;
    GFS_VALUE (cell, sv->v[1]) = (- m[1][0]*u + m[0][0]*v)/det;
}

static void implicit_darcy_3D (FttCell * cell, GfsSourceDarcy * s)
{

}

void gfs_source_darcy_implicit (GfsDomain * domain,
				   gdouble dt)
{
  GfsSourceDarcy * s;

  g_return_if_fail (domain != NULL);

  if ((s = gfs_has_source_darcy (domain))) {
    GfsSimulation * sim = GFS_SIMULATION (domain);
    gdouble olddt = sim->advection_params.dt;
    sim->advection_params.dt = dt;
    gfs_catch_floating_point_exceptions ();
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) (FTT_DIMENSION==3? 
						     implicit_darcy_3D : implicit_darcy_2D), s);
    if (gfs_restore_floating_point_exceptions ()) {
      gchar * c = g_strconcat ("\n", gfs_function_description (s->darcycoeff, FALSE), NULL);
      if (s->forchhicoeff)
	c = g_strconcat (c, "\n", gfs_function_description (s->forchhicoeff, FALSE), NULL);
      /* fixme: memory leaks */
      g_message ("floating-point exception in user-defined function(s):%s", c);
      exit (1);
    }
    sim->advection_params.dt = olddt;
  }
}

/* Initialize module */
const gchar gfs_module_name[] = "porous";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
 gfs_porous_class ();
 gfs_source_darcy_class ();
  return NULL;
}
