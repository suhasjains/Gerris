#ifndef _LPARTICLES_H_
#define _LPARTICLES_H_
#include "gfs.h"

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
        gdouble dt, dtn;
} ForceParams;

struct _Particle {

        FttVector pos,vel;
        guint id;
        gdouble density,volume;
        FttCell *cell;

};

struct _ForceCoefficients {

        guint init, fluidadv, RK4;

};

struct _LParticles {
  	
	/*< private >*/
  	GfsEvent parent;

  	/*< public >*/
  	/* add extra data here (if public) */
	GString *name;
        GfsVariable *density;
        GfsVariable *reynolds;

        GSList *particles;
        guint maxid;
        gboolean first_call;
        ForceCoefficients fcoeff;
	ForceParams pars;

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
LParticles * l_particles_new    (LParticlesClass * klass);



#endif
