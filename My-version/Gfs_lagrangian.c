#include "Gfs_lagrangian.h"


/* Initialize the lagrangian module*/
const gchar gfs_module_name[] = "lagrangian";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
   lagrangian_particles_class ();

        return NULL;
}

/*Reading lagrangian particles */
static void lagrangian_particles_read (GtsObject ** o, GtsFile * fp) {



}

/*Initializing lagrangian particles class*/
static void lagrangian_particles_class_init (GfsLagrangianParticlesClass * klass) {
  /* define new methods and overload inherited methods here */

	GTS_OBJECT_CLASS (klass)->read = lagrangian_particles_read;


}

/*Initializing lagrangian particles object*/
static void lagrangian_particles_init (GfsLagrangianParticles * object) {

  /* initialize object here */
}

/*Lagrangian particles class*/
GfsLagrangianParticlesClass * lagrangian_particles_class (void) {
static GfsLagrangianParticlesClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo lagrangian_particles_info = {
      "GfsLagrangianParticles",
      sizeof (GfsLagrangianParticles),
      sizeof (GfsLagrangianParticlesClass),
      (GtsObjectClassInitFunc) lagrangian_particles_class_init,
      (GtsObjectInitFunc) lagrangian_particles_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
                                  &lagrangian_particles_info);
  }

  return klass;
}

