#include "gts.h"
#include "event.h"



typedef struct _GfsLagrangianParticles GfsLagrangianParticles;
typedef struct _GfsLagrangianParticlesClass    GfsLagrangianParticlesClass;
GfsLagrangianParticlesClass * lagrangian_particles_class  (void);


struct _GfsLagrangianParticles {

        /*< private >*/
        GfsEvent parent;

        /*< public >*/
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

