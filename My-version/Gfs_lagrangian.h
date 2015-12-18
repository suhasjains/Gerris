#include "gts.h"
#include "event.h"

typedef struct _Particle Particle;
typedef struct _ForceCoefficients ForceCoefficients;
typedef struct _GfsLagrangianParticles GfsLagrangianParticles;
typedef struct _GfsLagrangianParticlesClass    GfsLagrangianParticlesClass;
GfsLagrangianParticlesClass * lagrangian_particles_class  (void);


/*Various forces acting on the particle due to its motion in the fluid*/
typedef struct {
        Particle *p;
        GfsVariable **u;
        ForceCoefficients *fcoeffs;
        gdouble dt;
        GfsLagrangianParticles *lagrangian;
} ForceParams;

struct _Particle {
        
	FttVector pos,vel;
        guint id;
        gdouble density,volume;  
	FttCell *cell;

};   

struct _ForceCoefficients {
	
	guint init, fluidadv;

};


struct _GfsLagrangianParticles {

        /*< private >*/
        GfsEvent parent;

        /*< public >*/
	GString *name;
	GfsVariable *density;
	GfsVariable *reynolds;

	GSList *particles;
	guint maxid;
	gboolean first_call;
	ForceCoefficients fcoeff;
};


struct _GfsLagrangianParticlesClass {
  /*< private >*/
  GfsEventClass parent_class;

  /*< public >*/
};

#define LAGRANGIAN_PARTICLES(obj)            GTS_OBJECT_CAST (obj,\
                                                 GfsLagrangianParticles,\
                                                 lagrangian_particles_class ())
#define LAGRANGIAN_PARTICLES_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
                                                 GfsLagrangianParticlesClass,\
                                                 lagrangian_particles_class())
#define IS_LAGRANGIAN_PARTICLES(obj)         (gts_object_is_from_class (obj,\
                                                 lagrangian_particles_class ()))

