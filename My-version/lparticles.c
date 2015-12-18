#include <gfs.h>
#include "lparticles.h"
/* LParticles: Object */


/*Particle data read method*/
static gboolean particle_read (GtsFile * fp, guint * id, FttVector * p, FttVector * v, gdouble *density, gdouble *volume) {

  	if (fp->type != GTS_INT) {
    		gts_file_error (fp, "expecting an integer (Id)");
    		return FALSE;
  	}

  	*id = atoi (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (p.x)");
    		return FALSE;
  	}
  	p->x = atof (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (p.y)");
    		return FALSE;
  	}
  	p->y = atof (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (p.z)");
    		return FALSE;
  	}
  	p->z = atof (fp->token->str);
  	gts_file_next_token (fp);

	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (v.x)");
    		return FALSE;
  	}
  	v->x = atof (fp->token->str);
  	gts_file_next_token (fp);
  
	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (v.y)");
    		return FALSE;
  	}
  	v->y = atof (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (v.z)");
    		return FALSE;
  	}
  	v->z = atof (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (density)");
    		return FALSE;
  	}
  	*density = atof (fp->token->str);
  	gts_file_next_token (fp);

  	if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    		gts_file_error (fp, "expecting a number (volume)");
    		return FALSE;
  	}
  	*volume = atof (fp->token->str);
  	gts_file_next_token (fp);

	return TRUE;
}


//Assigns values to variables after being read
static void assign_val_vars (guint * f, GtsFile *fp, GtsObject *o) {

  	gts_file_next_token (fp);
  	if (fp->type != '=') {
    		gts_file_error (fp, "expecting `='");
    		return;
  	}
  	gts_file_next_token (fp);

  	if (fp->type == GTS_INT)
    		*f = atoi(fp->token->str);

  	gts_file_next_token (fp);

}



static void l_particles_read (GtsObject ** o, GtsFile * fp)
{
  	/* call read method of parent */
  	if (GTS_OBJECT_CLASS (l_particles_class ())->parent_class->read) (* GTS_OBJECT_CLASS (l_particles_class ())->parent_class->read) (o, fp);
  	if (fp->type == GTS_ERROR)
    		return;


	/* do object-specific read here */
  	GfsLagrangianParticles * lagrangian = LAGRANGIAN_PARTICLES(*o);
	GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (lagrangian));


	//Name
  	if (fp->type == GTS_STRING) {
    		lagrangian->name = g_string_new(fp->token->str);
		//printf("Name of the new lagrangian variable is = %s\n", (lagrangian->name)->str);	
    		gts_file_next_token (fp);
  	}
  	else {
    		gts_file_error (fp, "expecting a string");
    		return;
  	}




	//density
	if (fp->type == GTS_STRING) {
		lagrangian->density = gfs_domain_get_or_add_variable (domain, fp->token->str, "Lagrangian Density variable");
    		if(lagrangian->density == NULL) {
      			gts_file_error (fp, "Expecting a defined variable name for fluid density evaluation.");
      			return;
    		}

		//printf("New domain variable added is ' %s '\n", lagrangian->density->name);
    		gts_file_next_token (fp);
  	}

	
	//Opening brackets
	if (fp->type != '{') {
    		gts_file_error (fp, "expecting an opening brace");
    		return;
  	}
	fp->scope_max++;
	gts_file_next_token (fp);	
	
	//Lagrangian particles initialize 
	lagrangian->particles = NULL;
	lagrangian->maxid = 0;


	//Particle data within brackets in the script
	while (fp->type != GTS_NONE && fp->type != '}') {					
		
		guint id;
		FttVector p,v;
		gdouble density, volume;
		if (!particle_read (fp, &id, &p, &v, &density, &volume)) 
        		return;
		
		Particle * particles = g_malloc(sizeof(Particle));
   		particles->id = id;
   		particles->pos = p;
   		particles->vel = v;
   		particles->density = density;
   		particles->volume = volume;

		//printf("Particle ID is %d\n",particles->id);

		//assigning the maximum of ids to maxid
		if(id > lagrangian->maxid)
		       	lagrangian->maxid = id;	

			
		lagrangian->particles = g_slist_append(lagrangian->particles, particles);		
		
		do
        		gts_file_next_token (fp);
      		while (fp->type == '\n');

	}

	if (fp->type != '}') {
    		gts_file_error (fp, "expecting a closing brace");
    		return;
  	}
	fp->scope_max--;
	gts_file_next_token (fp);	

	//lagrangian->reynolds = gfs_domain_get_or_add_variable (domain, "Rep", "Particle Reynolds number");
	//printf("%s\n", lagrangian->reynolds->name);

	
	//Opening a new file to input particles data
	if (fp->type != '{' && fp->type == GTS_STRING) {

		//printf("Particles input file is '%s'\n",fp->token->str);
    
	
		FILE *fptr = fopen(fp->token->str, "r");
    		GtsFile *fp1;
    		if(fptr == NULL) {
      			gts_file_error (fp, "cannot open file `%s'", fp->token->str);
      			return;
    		}
    		fp1 = gts_file_new (fptr);


		while (fp1->type != GTS_NONE) {
					
			guint id;
			FttVector p,v;
			gdouble density, volume;
			if (!particle_read (fp1, &id, &p, &v, &density, &volume)) {
        			gts_file_error (fp, "%s:%d:%d: %s", fp->token->str, fp1->line, fp1->pos, fp1->error);
        			return;
      			}
		
			Particle * particles = g_malloc(sizeof(Particle));
   			particles->id = id;
   			particles->pos = p;
   			particles->vel = v;
   			particles->density = density;
   			particles->volume = volume;

		//	printf("Particle ID is %d\n",particles->id);

			//assigning the maximum of ids to maxid
			if(id > lagrangian->maxid)
		        	lagrangian->maxid = id;	

			
			lagrangian->particles = g_slist_append(lagrangian->particles, particles);		
			
			do
        			gts_file_next_token (fp1);
      			while (fp1->type == '\n');

		}	
		gts_file_destroy (fp1);
 		fclose (fptr);
 
	}

 	gts_file_next_token (fp);



	//Coefficients	
	if (fp->type != '{') {
    		gts_file_error (fp, "expecting an opening brace");
    		return;
  	}
  
	fp->scope_max++;
  	gts_file_next_token (fp);

  	while (fp->type != GTS_ERROR && fp->type != '}') {

    		if (fp->type != GTS_STRING) {
      			gts_file_error (fp, "expecting a variable name");
      			return;
    		}
    		else {

			if(g_ascii_strcasecmp(fp->token->str, "init") == 0)
        			assign_val_vars (&lagrangian->fcoeff.init, fp, *o);

			else if(g_ascii_strcasecmp(fp->token->str, "fluidadv") == 0)
	        		assign_val_vars (&lagrangian->fcoeff.fluidadv, fp, *o);
			else{
        			gts_file_error (fp, "Not a valid variable");
   				return;
      			}

		}
	}
	if (fp->type != '}') {
    		gts_file_error (fp, "expecting a closing brace");
    		return;
  	}	
	
	fp->scope_max--;

 	/* do not forget to prepare for next read */
  	gts_file_next_token (fp);
}

static void l_particles_write (GtsObject * o, FILE * fp)
{
  	/* call write method of parent */
  	if (GTS_OBJECT_CLASS (l_particles_class ())->parent_class->write) (* GTS_OBJECT_CLASS (l_particles_class ())->parent_class->write) (o, fp);

  	/* do object specific write here */
	GfsLagrangianParticles * lagrangian = LAGRANGIAN_PARTICLES(o);
	
	fprintf(fp," %s ",lagrangian->name->str);
        if(lagrangian->density)
        	fprintf(fp," %s ",lagrangian->density->name);

	
	fprintf(fp," {");
	if(lagrangian->fcoeff.init == 1)
                fprintf(fp, " init = 1");
        if(lagrangian->fcoeff.fluidadv == 1)
                fprintf(fp, " fluidadv = 1");
	fprintf (fp," } \n");

	//write particles data
	GSList *i = lagrangian->particles;
        Particle * p;
        fputs (" { \n",fp);
        
	while(i) {
                p = (Particle *)(i->data);
                fprintf(fp,"%d %g %g %g %g %g %g %g %g ", p->id, p->pos.x, p->pos.y, p->pos.z,
                p->vel.x, p->vel.y, p->vel.z, p->density, p->volume);
                fprintf(fp,"\n");
                i = i->next;
        }
        fputs (" } \n",fp);

}

static gboolean l_particles_event (GfsEvent * event, GfsSimulation * sim)
{
  	if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (l_particles_class ())->parent_class)->event) (event, sim)) {
    

	/* do object-specific event here */

	//printf("dt=%g\n",sim->advection_params.dt);

    return TRUE;
  }
  return FALSE;
}

static void l_particles_class_init (LParticlesClass * klass)
{
  /* define new methods and overload inherited methods here */

  GFS_EVENT_CLASS (klass)->event = l_particles_event;
  GTS_OBJECT_CLASS (klass)->read = l_particles_read;
  GTS_OBJECT_CLASS (klass)->write = l_particles_write;
}

static void l_particles_init (LParticles * object)
{
  /* initialize object here */
}

LParticlesClass * l_particles_class (void)
{
  static LParticlesClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo l_particles_info = {
      "LParticles",
      sizeof (LParticles),
      sizeof (LParticlesClass),
      (GtsObjectClassInitFunc) l_particles_class_init,
      (GtsObjectInitFunc) l_particles_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &l_particles_info);
  }

  return klass;
}

LParticles * l_particles_new (LParticlesClass * klass)
{
  LParticles * object;

  object = L_PARTICLES (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/* Initialize the lagrangian module*/
const gchar gfs_module_name[] = "lparticles";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void) {

        l_particles_class ();

        return NULL;
}

