#ifndef _LPARTICLES_H_
#define _LPARTICLES_H_


/* LParticles: Header */

typedef struct _LParticles         LParticles;


struct _LParticles {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  /* add extra data here (if public) */
};


typedef struct _LParticlesClass    LParticlesClass;


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
