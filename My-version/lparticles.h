#ifndef _LPARTICLES_H_
#define _LPARTICLES_H_
//#include "gfs.h"
#include "event.h"
#include "source.h"
#include "vof.c"

/* LParticles: Header */
typedef struct _Particle Particle;
typedef struct _ForceCoefficients ForceCoefficients;
typedef struct _LParticles         LParticles;
typedef struct _LParticlesClass    LParticlesClass;

/*Various forces acting on the particle due to its motion in the fluid*/
typedef struct {
        Particle *p;
        GfsVariable **u, **un;
        ForceCoefficients *fcoeffs;
        FttVector *force;
	gdouble dt, dtn;
	gdouble fluid_rho, viscosity;
	LParticles *lagrangian;
} ForceParams;

typedef void (*force_pointer) (ForceParams *pars);


struct _Particle {

        FttVector pos,vel,acc,phiforce;
        guint id;
        gdouble density,volume;
        FttCell *cell;

};

struct _ForceCoefficients {

        guint init, fluidadv, RK4;
	gdouble cl, cd, cm;
	guint lift, drag, inertial, amf, buoy;


	GfsFunction *cdrag, *clift, *camf;


};

struct _LParticles {
  	
	/*< private >*/
  	GfsEvent parent;

  	/*< public >*/
  	/* add extra data here (if public) */
	GString *name;
        GfsVariable *density;
        GfsVariable *reynolds;
	GfsVariable **un;
	GfsVariable **couplingforce;
	GSList *forces;
        GSList *particles;
        guint maxid, idlast;
        gboolean first_call;
        ForceCoefficients fcoeff;

};



struct _LParticlesClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
  /* add extra methods here */
};


#define L_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
                                                 LParticles,\
                                                 l_particles_class ())
#define L_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
                                                 LParticlesClass,\
                                                 l_particles_class())
#define IS_L_PARTICLES(obj)         (gts_object_is_from_class (obj,\
                                                 l_particles_class ()))

LParticlesClass * l_particles_class  (void);
//LParticles * l_particles_new    (LParticlesClass * klass);




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




///*Droplet to Particles: Header*/
//typedef struct _GfsDropletToParticles         GfsDropletToParticles;
//typedef struct _GfsDropletToParticlesClass    GfsDropletToParticlesClass;
//
//struct _GfsDropletToParticles {
//  /*< private >*/
//  LParticles parent;
//  GfsVariable * v;
//
//  /*< public >*/
//  GfsFunction * fc;
//  GfsVariable * c;
//  gint min;
//  gdouble resetwith;
//  gdouble density;
//  /*Convert Particles*/
//  GfsFunction * fconvert;
//  GfsSurface *shape;
//  gint convert, maxlevel;
//};
//
//struct _GfsDropletToParticlesClass {
//  /*< private >*/
//  LParticlesClass parent_class;
//  /*< public >*/
//  /* add extra methods here */
//};
//
//
//#define DROPLET_TO_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
//                                                 GfsDropletToParticles,\
//                                                 droplet_to_particles_class ())
//#define DROPLET_TO_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
//                                                 GfsDropletToParticlesClass,\
//                                                 droplet_to_particles_class())
//#define IS_DROPLET_TO_PARTICLES(obj)         (gts_object_is_from_class (obj,\
//                                                 droplet_to_particles_class ()))
//
//
//GfsDropletToParticlesClass * droplet_to_particles_class (void);
/********************/

/* DropletToParticles: Header */

typedef struct _DropletToParticles         DropletToParticles;

struct _DropletToParticles {
  /*< private >*/
  LParticles parent;

  /*< public >*/
  /* add extra data here (if public) */
};

typedef struct _DropletToParticlesClass    DropletToParticlesClass;

struct _DropletToParticlesClass {
  /*< private >*/
  LParticlesClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define DROPLET_TO_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
					         DropletToParticles,\
					         droplet_to_particles_class ())
#define DROPLET_TO_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 DropletToParticlesClass,\
						 droplet_to_particles_class())
#define IS_DROPLET_TO_PARTICLES(obj)         (gts_object_is_from_class (obj,\
						 droplet_to_particles_class ()))

DropletToParticlesClass * droplet_to_particles_class  (void);
DropletToParticles * droplet_to_particles_new    (DropletToParticlesClass * klass);




/*Shared Methods*/

void mpi_particle_numbering(GfsDomain *domain, LParticles *lagrangian);

typedef struct {
  guint maxlevel;
  GfsDomain * domain;
  GfsGenericSurface * surface;
  gboolean check;
} RefineCut;

void save_p_solid (FttCell * cell, gpointer * data);

void restore_p_solid (FttCell * cell, gpointer * data);

void gfs_domain_assign_fraction (GfsDomain * domain,
                               GfsGenericSurface * s,
                                 GfsVariable * c, gdouble resetwith,Particle *p);

void refine_implicit_p_cell (FttCell * cell, RefineCut * p);

void add_to_prev_void (FttCell * cell, gpointer * data);

/********************/




#endif
