#include "lparticles.h"

/* LParticles: Object */


/*Updating the particle velocity due to forces acting on the particle*/
static void compute_particle_velocity (Particle * p, double dt) {

        p->vel.x +=  dt* p->acc.x;
        p->vel.y +=  dt* p->acc.y;
        #if !FTT_2D
                p->vel.z +=  dt* p->acc.z;
        #endif
}

/*Updating the particle acceleration*/
static void compute_particle_acceleration(Particle *p, FttVector force, gdouble cm, gdouble fluid_rho) {

        //Taking a component of the added mass force to the LHS of the momentum equation of the particle
        p->acc.x += force.x/(p->density + fluid_rho*cm);
        p->acc.y += force.y/(p->density + fluid_rho*cm);
        #if !FTT_2D
                p->acc.z += force.z/(p->density + fluid_rho*cm);
        #endif
}

/*Calculation of the lift force acting on the unit volume of the particle*/
static void compute_lift_force(LParticles * lagrangian) {

        Particle *p = lagrangian->pars->p;
        GfsVariable **u = pars->u;
        gdouble fluid_rho = pars->fluid_rho;
        FttVector *force = pars->force;
        ForceCoefficients *fcoeffs = pars->fcoeffs;

        /*Interpolating fluid velocity to obtain the value at the position of the particle*/
        FttVector fluid_vel;
        fluid_vel.x = gfs_interpolate(p->cell, p->pos, u[0]);
        fluid_vel.y = gfs_interpolate(p->cell, p->pos, u[1]);
        #if !FTT_2D
                fluid_vel.z = gfs_interpolate(p->cell, p->pos, u[2]);
        #endif

        /*Subtrating fluid from particle velocity vectors to obtain the relative velocity */
        FttVector relative_vel;
        subs_fttvectors(&fluid_vel, &p->vel, &relative_vel);

        /*Calculation of vorticity */
        FttVector vorticity;
        gfs_vorticity_vector (p->cell, u, &vorticity);

        /*lift coefficient*/
        fcoeffs->cl = 0.5;

        /*lift force calculation*/
        #if FTT_2D
        force->x = fluid_rho*fcoeffs->cl*relative_vel.y*vorticity.z;
        force->y = -fluid_rho*fcoeffs->cl*relative_vel.x*vorticity.z;
        #else
        force->x = fluid_rho*fcoeffs->cl*(relative_vel.y*vorticity.z - relative_vel.z*vorticity.y);
        force->y = fluid_rho*fcoeffs->cl*(relative_vel.z*vorticity.x - relative_vel.x*vorticity.z);
        force->z = fluid_rho*fcoeffs->cl*(relative_vel.x*vorticity.y - relative_vel.y*vorticity.x);
        #endif

    return;
}




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
  	LParticles * lagrangian = L_PARTICLES(*o);
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
			
			else if(g_ascii_strcasecmp(fp->token->str, "RK4") == 0)
	        		assign_val_vars (&lagrangian->fcoeff.RK4, fp, *o);
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
	LParticles * lagrangian = L_PARTICLES(o);
	
	fprintf(fp," %s ",lagrangian->name->str);
        if(lagrangian->density)
        	fprintf(fp," %s ",lagrangian->density->name);

	
	fprintf(fp," {");
	if(lagrangian->fcoeff.init == 1)
                fprintf(fp, " init = 1");
        if(lagrangian->fcoeff.fluidadv == 1)
                fprintf(fp, " fluidadv = 1");
	if(lagrangian->fcoeff.RK4 == 1)
		fprintf(fp, " RK = 1 ");
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



/*l particle object destroy */
static void l_particles_destroy (GtsObject * object) {

        /* do not forget to call destroy method of the parent */
        (* GTS_OBJECT_CLASS (l_particles_class ())->parent_class->destroy) (    object);

        /* do object-specific cleanup here */
        LParticles * lagrangian = L_PARTICLES(object);


        if(lagrangian->particles)
                g_slist_foreach (lagrangian->particles, (GFunc) g_free, NULL);

        g_slist_free(lagrangian->particles);

        g_string_free(lagrangian->name, TRUE);
                

}


/*Updating the particle position using Euler scheme*/
static void advect_particles(Particle * p, double dt) {

        p->pos.x +=  dt* p->vel.x;
        p->pos.y +=  dt* p->vel.y;
        #if !FTT_2D
                p->pos.z +=  dt* p->vel.z;
        #endif
}

/*Updating the particle position using RK4 scheme*/
static void advect_particles_RK4(Particle * p, double dtn, GfsVariable ** u, GfsVariable ** un) {

	double k1, k2, k3, k4; 

	FttVector RK4k1, RK4k2, RK4k3, RK4k4;  

	if(p->cell!=NULL) {

		RK4k4 = RK4k3 = RK4k2 = RK4k1 = p->pos;		

		k1 = gfs_interpolate(p->cell, RK4k1, un[0]);
		RK4k2.x = RK4k1.x + k1*dtn/2;
		k2 = ( gfs_interpolate(p->cell, RK4k2, un[0]) + gfs_interpolate(p->cell, RK4k2, u[0]) )/2;
		RK4k3.x = RK4k1.x + k2*dtn/2;
		k3 = ( gfs_interpolate(p->cell, RK4k3, un[0]) + gfs_interpolate(p->cell, RK4k3, u[0]) )/2;
		RK4k4.x = RK4k1.x + k3*dtn;
		k4 = gfs_interpolate(p->cell, RK4k4, u[0]);
			
		//printf("k1=%g, k2=%g, k3=%g, k4=%g, dtn=%g\n", k1, k2, k3, k4, dtn);	
	
		p->pos.x +=  dtn*(k1 + 2*k2 + 2*k3 + k4)/6;	
		
		RK4k4 = RK4k3 = RK4k2 = RK4k1;

		k1 = gfs_interpolate(p->cell, RK4k1, un[1]);
		RK4k2.y = RK4k1.y + k1*dtn/2;
		k2 = ( gfs_interpolate(p->cell, RK4k2, un[1]) + gfs_interpolate(p->cell, RK4k2, u[1]) )/2;
		RK4k3.y = RK4k1.y + k2*dtn/2;
		k3 = ( gfs_interpolate(p->cell, RK4k3, un[1]) + gfs_interpolate(p->cell, RK4k3, u[1]) )/2;
		RK4k4.y = RK4k1.y + k3*dtn;
		k4 = gfs_interpolate(p->cell, RK4k4, u[1]);
		
        	

		p->pos.y +=  dtn*(k1 + 2*k2 + 2*k3 + k4)/6;


        	#if !FTT_2D
		
			RK4k4 = RK4k3 = RK4k2 = RK4k1;
	
			k1 = gfs_interpolate(p->cell, RK4k1, un[2]);
			RK4k2.z = RK4k1.z + k1*dtn/2;
			k2 = ( gfs_interpolate(p->cell, RK4k2, un[2]) + gfs_interpolate(p->cell, RK4k2, u[2]) )/2;
			RK4k3.z = RK4k1.z + k2*dtn/2;
			k3 = ( gfs_interpolate(p->cell, RK4k3, un[2]) + gfs_interpolate(p->cell, RK4k3, u[2]) )/2;
			RK4k4.z = RK4k1.z + k3*dtn;
			k4 = gfs_interpolate(p->cell, RK4k4, u[2]);
		

                	p->pos.z +=  dtn*(k1 + 2*k2 + 2*k3 + k4)/6;
				
		#endif

	}
}

/*Equating current fluid velocity*/
static void fluidadvect_particles(Particle * p, GfsVariable **u) {

        if(p->cell!=NULL) {
                p->vel.x = gfs_interpolate(p->cell, p->pos, u[0]);
                p->vel.y = gfs_interpolate(p->cell, p->pos, u[1]);
                #if !FTT_2D
                        p->vel.z = gfs_interpolate(p->cell, p->pos, u[2]);
                #endif
        }
}

/*Equating previous fluid velocity*/
static void fluidadvect_particles_RK4(Particle * p, GfsVariable **un) {


        if(p->cell!=NULL) {
                p->vel.x = gfs_interpolate(p->cell, p->pos, un[0]);
                p->vel.y = gfs_interpolate(p->cell, p->pos, un[1]);
                #if !FTT_2D
                        p->vel.z = gfs_interpolate(p->cell, p->pos, un[2]);
                #endif
        }
}

/* Initializes particle velocity*/
static void init_particles(LParticles * l, GfsDomain *domain) {

        GSList *i = l->particles;
        Particle *p;
        GfsVariable ** u = gfs_domain_velocity (domain);

        while(i) {
                p = (Particle *) i->data;
                p->cell = gfs_domain_locate(domain, p->pos, -1, NULL);

                if(p->cell!=NULL) {
                        p->vel.x = gfs_interpolate(p->cell, p->pos, u[0]);
                        p->vel.y = gfs_interpolate(p->cell, p->pos, u[1]);
                        #if !FTT_2D
                                p->vel.z = gfs_interpolate(p->cell, p->pos, u[2]);
                        #endif
                }
        i = i->next;
        }
}


static gboolean l_particles_event (GfsEvent * event, GfsSimulation * sim)
{
  	if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (l_particles_class ())->parent_class)->event) (event, sim)) {
    
		/* do object-specific event here */
		GfsDomain * domain = GFS_DOMAIN (sim);
        	LParticles *lagrangian = L_PARTICLES(event);

                if(lagrangian->first_call) {
			
		//printf("Yeah");
	
                        GSList *i = lagrangian->particles;
                        lagrangian->maxid = 0;
                        while(i) {
                                Particle *p = (Particle *)(i->data);
                                lagrangian->maxid = MAX(lagrangian->maxid, p->id);
                                p->cell = gfs_domain_locate(domain, p->pos, -1, NULL)    ;
                                i = i->next;

                                //remove the particle if outside domain
                                if(!p->cell) {
                                        printf("Particle %d at %g %g is outside domain\n",p->id, p->pos.x, p->pos.y);
                                        lagrangian->particles = g_slist_remove(lagrangian->particles, p);
                                        g_free(p);
                                }
                        }
                        //Initializing particle velocity
                        if(lagrangian->fcoeff.init == 1)
                                init_particles(lagrangian, domain);


                }
		
		//Fetch simulation parameters
                GSList *i = lagrangian->particles;
                Particle * p;
		FttVector force;
 	       	force_pointer f;
        	GSList *fs;


                //Fetch force parameters
                gdouble dt = lagrangian->pars.dt = sim->advection_params.dt;
		gdouble dtn = lagrangian->pars.dtn;
                lagrangian->pars.fcoeffs = &lagrangian->fcoeff;
                GfsVariable ** u = lagrangian->pars.u = gfs_domain_velocity (domain);
		GfsVariable ** un;
		if(lagrangian->first_call) 	
			un = lagrangian->pars.un = u; 	
		else 	
			un = lagrangian->pars.un;	
        	GfsSourceDiffusion *d = source_diffusion_viscosity(u[0]);
		pars->force = &force;


                //Looping over all particles
                i = lagrangian->particles;
                while(i) {
                        p = (Particle *)(i->data);
                        p->cell = gfs_domain_locate(domain, p->pos, -1, NULL);

                        //Remove particle if outside domain
                        if(!p->cell) {
                                lagrangian->particles = g_slist_remove(lagrangian->particles, p);
                                g_free(p);
                        }
                        else {
                                //Making velocity equal to fluid velocity
                                if(lagrangian->fcoeff.fluidadv == 1) {
                                        fluidadvect_particles(p, u);
				//	printf("Using fluid velocity for advection\n");
                                }
				else {
					gdouble viscosity;
         				if(d)
           					viscosity = gfs_diffusion_cell(d->D, p->cell);
         				else
           					viscosity = 0.;

         				gdouble fluid_rho = sim->physical_params.alpha ? 1./
           				gfs_function_value(sim->physical_params.alpha,p->cell) : 1.;

         				p->acc.x = 0.;
         				p->acc.y = 0.;
         				p->acc.z = 0.;

         				lagrangian->pars->p = p;
         				lagrangian->pars->fluid_rho = fluid_rho;
         				lagrangian->pars->viscosity = viscosity;

         				fs = lagrangian->forces;
         				while(fs){
           					f = (force_pointer)fs->data;
           					(f)(lagrangian);
           					compute_particle_acceleration(p, force, lagrangian->fcoeff.cm, fluid_rho);
           					fs = fs->next;
         				}
           			
					//Particle velocity from Newton's equation	
					compute_particle_velocity (p, dt);


					//Computing coupling force
         				//if(p->cell != NULL){
           				//	particles_force(lagrangian);
           				//	compute_coupling_force(p, lagrangian->couplingforce);
					//}
				 }
			}
                        i = i->next;
                }

		//Looping over all particles
                i = lagrangian->particles;
                while(i) {
                        p = (Particle *)(i->data);
                        p->cell = gfs_domain_locate(domain, p->pos, -1, NULL);

                        //Remove particle if outside domain
                        if(!p->cell) {
                                lagrangian->particles = g_slist_remove(lagrangian->particles, p);
                                g_free(p);
                        }
                        else {
				if(lagrangian->fcoeff.RK4 == 1) {
                                	advect_particles_RK4(p, dtn, u, un);
					//printf("Doing RK4 advection\n");
                                }

				else {
                                	//Advect particle according to particle velocity
                                	advect_particles(p, dt);
					//printf("Doing normal advection\n");
                        	}
			}
                    	
			printf("%d %g %g %g %g %g\n",p->id, dtn, sim->time.t-dtn, p->pos.x, p->pos.y, p->pos.z);
                        i = i->next;
                }

		//Storing to be used next time step
		lagrangian->pars.un = u;
		lagrangian->pars.dtn = dt;
                lagrangian->first_call = FALSE;
    		return TRUE;
  	}
  	return FALSE;
}



/*Initializing l particles class*/
static void l_particles_class_init (LParticlesClass * klass)
{
  	/* define new methods and overload inherited methods here */
  	GFS_EVENT_CLASS (klass)->event = l_particles_event;
  	GTS_OBJECT_CLASS (klass)->read = l_particles_read;
  	GTS_OBJECT_CLASS (klass)->write = l_particles_write;
	GTS_OBJECT_CLASS (klass)->destroy = l_particles_destroy;
}


/*Initializing l particles object*/
static void l_particles_init (LParticles * object)
{
  	/* initialize object here */
	object->particles = NULL;
        object->first_call = TRUE;

        object->fcoeff.init = 0;
        object->fcoeff.fluidadv = 0;
	object->fcoeff.RK4 = 0;
}


/*L Particles class*/
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

