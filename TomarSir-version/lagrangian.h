#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ftt.h"
#include "domain.h"
#include "output.h"
#include "fluid.h"
#include "solid.h"
#include "simulation.h"
#include "source.h"
#include "boundary.h"
#include "vof.h"
#include "spatial.h"

/* Particle */
typedef struct _Particle Particle;
struct _Particle{

  FttVector pos,acc,vel,phiforce;
  guint id;
  gdouble density,volume;  
  FttCell *cell;
 
  /*IBM params & grid*/
  gint move;
  FttVector q;
};

/* LagrangianParticles: Header */

typedef struct _GfsLagrangianParticles         GfsLagrangianParticles;

typedef struct _ForceCoefficients ForceCoefficients;

struct _ForceCoefficients {

  gdouble cl, cd, cm, cf; 
  guint lift, drag, buoy, inertial, faxen, amf, basset, init, fluidadv;
  guint collision, merging;
  guint imsolid;
  gdouble bending, tension;

  GfsFunction * cdrag, * clift, * camf , * cfaxen;
};

struct _GfsLagrangianParticles {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  guint n;
  gdouble time;
  ForceCoefficients   fcoeff;
  GSList * particles;
  GSList * forces;
  GfsVariable **un;
  GfsVariable **couplingforce;
  GfsVariable * density;

  guint maxid, idlast;
  GfsVariable * reynolds, * urel, * vrel, * wrel, * pdia;
  GString *name;

  gboolean first_call;
  GHashTable * cell_particles;
  /* add extra data here (if public) */
};


typedef struct _GfsLagrangianParticlesClass    GfsLagrangianParticlesClass;

struct _GfsLagrangianParticlesClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
  /* add extra methods here */
};


#define LAGRANGIAN_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsLagrangianParticles,\
					         lagrangian_particles_class ())
#define LAGRANGIAN_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsLagrangianParticlesClass,\
						 lagrangian_particles_class())
#define IS_LAGRANGIAN_PARTICLES(obj)         (gts_object_is_from_class (obj,\
						 lagrangian_particles_class ()))

GfsLagrangianParticlesClass * lagrangian_particles_class  (void);


/*GfsFeedParticles: Header*/
typedef struct _GfsFeedParticles         GfsFeedParticles;

struct _GfsFeedParticles {
  /*< private >*/
  GfsLagrangianParticles parent;

  gdouble feed, time;
  GSList *particles;
  /*< public >*/
  
};

typedef struct _GfsFeedParticlesClass    GfsFeedParticlesClass;
struct _GfsFeedParticlesClass {
  /*< private >*/
  GfsLagrangianParticlesClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define FEED_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsFeedParticles,\
					         feed_particles_class ())
#define FEED_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsFeedParticlesClass,\
						 feed_particles_class())
#define IS_FEED_PARTICLES(obj)         (gts_object_is_from_class (obj,\
						 feed_particles_class ()))

GfsFeedParticlesClass * feed_particles_class (void);

/*GfsSourceLagrangian : Header*/
typedef struct _GfsSourceLagrangian         GfsSourceLagrangian;

struct _GfsSourceLagrangian {
  /*< private >*/
  GfsSourceVelocity parent;

  GString * name;
  /*< public >*/
  
};

#define GFS_SOURCE_LAGRANGIAN(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceLagrangian,\
					         gfs_source_lagrangian_class ())

#define GFS_IS_SOURCE_LAGRANGIAN(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_lagrangian_class ()))

GfsSourceGenericClass * gfs_source_lagrangian_class (void);

/*Droplet to Particles: Header*/
typedef struct _GfsDropletToParticles         GfsDropletToParticles;
typedef struct _GfsDropletToParticlesClass    GfsDropletToParticlesClass;
 
struct _GfsDropletToParticles {
  /*< private >*/
  GfsLagrangianParticles parent;
  GfsVariable * v;

  /*< public >*/
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble resetwith;
  gdouble density;
  /*Convert Particles*/
  GfsFunction * fconvert;
  GfsSurface *shape;
  gint convert, maxlevel;
};

struct _GfsDropletToParticlesClass {
  /*< private >*/
  GfsLagrangianParticlesClass parent_class;
  /*< public >*/
  /* add extra methods here */
};


#define DROPLET_TO_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDropletToParticles,\
					         droplet_to_particles_class ())
#define DROPLET_TO_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsDropletToParticlesClass,\
						 droplet_to_particles_class())
#define IS_DROPLET_TO_PARTICLES(obj)         (gts_object_is_from_class (obj,\
						 droplet_to_particles_class ()))


GfsDropletToParticlesClass * droplet_to_particles_class (void);



/*Convert Particles: Header*/

typedef struct _GfsConvertParticles        GfsConvertParticles;
typedef struct _GfsConvertParticlesClass    GfsConvertParticlesClass;
 
struct _GfsConvertParticles {
  /*< private >*/
  GfsFeedParticles parent;

  /*< public >*/
  /*Convert Particles*/
  GfsFunction * fconvert;
  GfsSurface *shape;
  gint maxlevel;
  GfsVariable * c;
  gdouble resetwith;
};

struct _GfsConvertParticlesClass {
  /*< private >*/
  GfsFeedParticlesClass parent_class;
  /*< public >*/
  /* add extra methods here */
};


#define CONVERT_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
					         GfsConvertParticles,\
					         convert_particles_class ())
#define CONVERT_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsConvertParticlesClass,\
						 convert_particles_class())
#define IS_CONVERT_PARTICLES(obj)         (gts_object_is_from_class (obj,\
						 convert_particles_class ()))


GfsConvertParticlesClass * convert_particles_class (void);


/********************/

/*Output Particles: Header*/

typedef struct _GfsOutputLParticle     GfsOutputLParticle;
typedef struct _GfsOutputLParticleClass    GfsOutputLParticleClass;
 
struct _GfsOutputLParticle {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  GString *name;
  GfsLagrangianParticles *lagrangian; 

};



#define GFS_OUTPUT_LPARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputLParticle,\
					         gfs_output_lparticle_class ())
#define GFS_OUTPUT_LPARTICLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsOutputClass,\
						 gfs_output_lparticle_class())
#define GFS_IS_OUTPUT_LPARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_lparticle_class ()))


GfsOutputClass * gfs_output_lparticle_class (void);


/********************/
typedef struct {
  guint maxlevel;
  GfsDomain * domain;
  GfsGenericSurface * surface;
  gboolean check;
} RefineCut;

/*Shared Methods*/

void mpi_particle_numbering(GfsDomain *domain, GfsLagrangianParticles *lagrangian);

void save_p_solid (FttCell * cell, gpointer * data);

void restore_p_solid (FttCell * cell, gpointer * data);

void gfs_domain_assign_fraction (GfsDomain * domain,
			       GfsGenericSurface * s,
				 GfsVariable * c, gdouble resetwith,Particle *p);

void refine_implicit_p_cell (FttCell * cell, RefineCut * p);

void add_to_prev_void (FttCell * cell, gpointer * data);
