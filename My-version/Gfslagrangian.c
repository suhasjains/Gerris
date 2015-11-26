/*Subtration of two vectors*/
static void subs_fttvectors (FttVector *a, FttVector *b, FttVector *result) {
	result->x = a->x - b->x;
	result->y = a->y - b->y;
	result->z = a->z - b->z;
}

/*Calculation of the vorticity vector*/
static void gfs_vorticity_vector (FttCell *cell, GfsVariable **v, FttVector *vort) {
	
	gdouble size;
		
	if (cell == NULL) return;
	if (v == NULL) return;

	size = ftt_cell_size (cell);

	#if FTT_2D
	
	vort->x = 0.;
	vort->y = 0.;
	vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) - gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;

	#else /*FTT_3D*/

	vort->x = (gfs_center_gradient (cell, FTT_Y, v[2]->i) - gfs_center_gradient (cell, FTT_Z, v[1]->i))/size;
   	vort->y = (gfs_center_gradient (cell, FTT_Z, v[0]->i) - gfs_center_gradient (cell, FTT_X, v[2]->i))/size;
   	vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) - gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
	 #endif
} 


  
/*Various forces acting on the particle due to its motion in the fluid*/
typedef struct {
	Particle *p;
	GfsVariable **u, **un;
	FttVector *force;
	ForceCoefficients *fcoeffs;
	gdouble fluid_rho, viscosity, dt;
	GfsLagrangianParticles *lagrangian;
} ForceParams;

/*Calculation of the lift force acting on the particle*/
static void compute_lift_force(ForceParams * pars) {
  
	Particle *p = pars->p;
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

	/*lift force calcualtion*/
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

/*Calculation of drag force acting on the particle */
static void compute_drag_force (ForceParams *pars) {

	Particle *p = pars->p;
	GfsVariable **u = pars->u;
	gdouble fluid_rho = pars->fluid_rho;
	gdouble viscosity = pars->viscosity;
	ForceCoefficients *fcoeffs = pars->fcoeffs;
	FttVector *force = pars->force;	
 
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
			
	gdouble radius = pow(p->volume/M_PI, 1./2.);
	#if !FTT_2D
	radius = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
	#endif


	/* Calculation of the magnitude of the relative velocity*/
 	#if !FTT_2D			
   	gdouble norm_relative_vel = sqrt(relative_vel.x*relative_vel.x + 
					relative_vel.y*relative_vel.y +
                                    	relative_vel.z*relative_vel.z);
 	#else
   	gdouble norm_relative_vel = sqrt(relative_vel.x*relative_vel.x +
                                    	relative_vel.y*relative_vel.y);
 	#endif

	/* Calculation of Reynolds number */	 
	gdouble Re;
   	if(viscosity == 0){
     		force->x = 0.;
     		force->y = 0.;
     		force->z = 0.;
     	return;
   	}
   	else
     	Re = 2.*norm_relative_vel*radius*fluid_rho/viscosity;

		
	

	if(fcoeffs->cdrag){
  
     	GFS_VARIABLE(p->cell, pars->lagrangian->reynolds->i) = Re;
  
     	GFS_VARIABLE(p->cell, pars->lagrangian->urel->i) = relative_vel.x;
  
     	GFS_VARIABLE(p->cell, pars->lagrangian->vrel->i) = relative_vel.y;
 
	#if !FTT_2D
     	GFS_VARIABLE(p->cell, pars->lagrangian->wrel->i) = relative_vel.z;
 	#endif
  
     	GFS_VARIABLE(p->cell, pars->lagrangian->pdia->i) = 2.0*radius;
 
     	fcoeffs->cd = gfs_function_value (fcoeffs->cdrag, p->cell);
     	force->x = fcoeffs->cd*relative_vel.x*fluid_rho;
     	force->y = fcoeffs->cd*relative_vel.y*fluid_rho;
     	force->z = fcoeffs->cd*relative_vel.z*fluid_rho;
 
     	return;
   	}

	
	if(Re < 1e-8){
    		force->x = 0.;
     		force->y = 0.;
     		force->z = 0.;
     		return;
   	}
   	else if(Re < 50.0)
     		fcoeffs->cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
   	else
     		fcoeffs->cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;
 
   		force->x = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.x*fluid_rho;
   		force->y = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.y*fluid_rho;
   		force->z = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.z*fluid_rho;
 }


/*Calculate the buoyant forces*/
static void compute_buoyant_force(ForceParams * pars) {

	Particle *p = pars->p;
   	GfsVariable **u = pars->u;
   	gdouble fluid_rho = pars->fluid_rho;
   	FttVector * force = pars->force;
   	gdouble g[3];
   	FttComponent c;
 
   	for(c = 0; c < FTT_DIMENSION; c++){
     		g[c] = 0.;
     		if (u[c]->sources) {
      	 		GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;
       
       		while (i) {
         		if (GFS_IS_SOURCE (i->data)) {
           			g[c] += gfs_function_value (GFS_SOURCE ((GfsSourceGeneric *) i->data)->intensity,
                                   				p->cell);
         		}
         		i = i->next;
       		}
     		}
   	}

   	force->x = (p->density-fluid_rho)*g[0];
   	force->y = (p->density-fluid_rho)*g[1];
 	#if !FTT_2D
   	force->z = (p->density-fluid_rho)*g[2];
 	#endif

} 


/*Calculate the inertial force*/
static void compute_inertial_force(ForceParams * pars) { 	
	Particle *p = pars->p;
   	GfsVariable **u = pars->u;
   	gdouble fluid_rho = pars->fluid_rho;
   	FttVector * force = pars->force;
 
   	GfsVariable **un = pars->un;
   	gdouble dt = pars->dt;
  
   	gdouble size = ftt_cell_size(p->cell);
 
 
   	FttVector fluid_vel;
   	FttVector fluid_veln;
   	fluid_vel.x = gfs_interpolate(p->cell, p->pos, u[0]);
   	fluid_vel.y = gfs_interpolate(p->cell, p->pos, u[1]);
 
   	fluid_veln.x = gfs_interpolate(p->cell, p->pos, un[0]);
   	fluid_veln.y = gfs_interpolate(p->cell, p->pos, un[1]);
 	#if !FTT_2D
   		fluid_vel.z = gfs_interpolate(p->cell, p->pos, u[2]);
   		fluid_veln.z = gfs_interpolate(p->cell, p->pos, un[2]);
 	#endif
 
 
   	if(dt > 0.){
     		force->x = fluid_rho*(fluid_vel.x - fluid_veln.x)/dt;
     		force->y = fluid_rho*(fluid_vel.y - fluid_veln.y)/dt;
 	#if !FTT_2D
     		force->z = fluid_rho*(fluid_vel.z - fluid_veln.z)/dt;
 	#endif
   	}
 
   	FttComponent c;
   	for(c=0;c < FTT_DIMENSION;c++){
     		force->x += fluid_rho*gfs_center_gradient(p->cell, c, u[0]->i)*
       			GFS_VARIABLE(p->cell, u[c]->i)/size;
     		force->y += fluid_rho*gfs_center_gradient(p->cell, c, u[1]->i)*
       			GFS_VARIABLE(p->cell, u[c]->i)/size;
 	#if !FTT_2D
     		force->z += fluid_rho*gfs_center_gradient(p->cell, c, u[2]->i)*
      	 		GFS_VARIABLE(p->cell, u[c]->i)/size;
 	#endif
   }
 }

/*Calculate the added mass force*/  //particle velocity not subtracted
static void compute_amf_force(ForceParams * pars) {

	Particle *p = pars->p;
   	GfsVariable **u = pars->u;
   	gdouble fluid_rho = pars->fluid_rho;
   	FttVector * force = pars->force;
   	ForceCoefficients * fcoeffs = pars->fcoeffs;
 
   	GfsVariable **un = pars->un;
   	gdouble dt = pars->dt;
 
   	gdouble size = ftt_cell_size(p->cell);
     
   	FttVector fluid_vel;
   	FttVector fluid_veln;
   	fluid_vel.x = gfs_interpolate(p->cell, p->pos, u[0]);
   	fluid_vel.y = gfs_interpolate(p->cell, p->pos, u[1]);
 
   	fluid_veln.x = gfs_interpolate(p->cell, p->pos, un[0]);
   	fluid_veln.y = gfs_interpolate(p->cell, p->pos, un[1]);
 	#if !FTT_2D
   		fluid_vel.z = gfs_interpolate(p->cell, p->pos, u[2]);
   		fluid_veln.z = gfs_interpolate(p->cell, p->pos, un[2]);
 	#endif
   
 
   	if(dt > 0.){
     		force->x = fluid_rho*(fluid_vel.x - fluid_veln.x)/dt;
     		force->y = fluid_rho*(fluid_vel.y - fluid_veln.y)/dt;
 	#if !FTT_2D
     		force->z = fluid_rho*(fluid_vel.z - fluid_veln.z)/dt;
 	#endif
  	}
 
   	FttComponent c;
   	for(c = 0; c < FTT_DIMENSION; c++){
     		force->x += fluid_rho*gfs_center_gradient(p->cell, c, u[0]->i)*
       			GFS_VARIABLE(p->cell, u[c]->i)/size;
     		force->y += fluid_rho*gfs_center_gradient(p->cell, c, u[1]->i)*
       			GFS_VARIABLE(p->cell, u[c]->i)/size;
 	#if !FTT_2D
     		force->z += fluid_rho*gfs_center_gradient(p->cell, c, u[2]->i)*
       			GFS_VARIABLE(p->cell, u[c]->i)/size;
 	#endif
   }
   
   	force->x = fcoeffs->cm*force->x;
   	force->y = fcoeffs->cm*force->y;
   	force->z = fcoeffs->cm*force->z;
}


/* Summation of the forces acting on the particle*/ //Only inertial and buoyancy forces are added
static void particles_force(ForceParams * pars) {
   	Particle *p = pars->p;
   	FttVector * force = pars->force;
	gdouble fluid_rho = pars->fluid_rho;
   	ForceCoefficients * fcoeffs = pars->fcoeffs;
 
   	if(p->cell==NULL)
     		return;
 
 
   	gdouble p3dvolume = p->volume;
 
   	p->phiforce.x = -p->acc.x*(p->density)* p3dvolume;
   	p->phiforce.y = -p->acc.y*(p->density)* p3dvolume;
   	p->phiforce.z = -p->acc.z*(p->density)* p3dvolume;
 
   	if(fcoeffs->buoy == 1){
     		compute_buoyant_force(pars);
     		p->phiforce.x += ((force->x)* p3dvolume);
     		p->phiforce.y += ((force->y)* p3dvolume);
     		p->phiforce.z += ((force->z)* p3dvolume);
   	}
 
   	if(fcoeffs->inertial == 1){
     		compute_inertial_force(pars);
     		p->phiforce.x += ((force->x)* p3dvolume);
     		p->phiforce.y += ((force->y)* p3dvolume);
     		p->phiforce.z += ((force->z)* p3dvolume);
   	}
 
	//recheck - fluid density is being divided 
   	p->phiforce.x /=(pars->fluid_rho);
   	p->phiforce.y /=(pars->fluid_rho);
   	p->phiforce.z /=(pars->fluid_rho);
 
}

/*Smoothening the force using Gaussian smoothing Exp(-x^2/sigma^2)*/
	static void compute_smooth_force(FttCell *cell, gpointer * data)
 	{
   	if(!cell) return;
 
   	Particle *p = (Particle *)(data[0]);
   	GfsVariable **f = data[1];
   	gdouble sigma = *(gdouble *)data[2];
   	FttVector pos;
   	ftt_cell_pos(cell, &pos);
 
 
   	gdouble dist = (p->pos.x - pos.x)*(p->pos.x - pos.x)
     			+ (p->pos.y - pos.y)*(p->pos.y - pos.y);
 
 	#if !FTT_2D
     		dist += (p->pos.z - pos.z)*(p->pos.z - pos.z);
 	#endif
     		dist = exp(-dist/(sigma*sigma))/(2.*M_PI*sigma*sigma);
 	#if !FTT_2D
     		dist /= (pow(2.*M_PI,0.5)*sigma);
 	#endif
    		dist /= (pow(2.*M_PI,0.5)*sigma);
 
   	GFS_VARIABLE(cell, f[0]->i) += p->phiforce.x*dist;
   	GFS_VARIABLE(cell, f[1]->i) += p->phiforce.y*dist;
 	#if !FTT_2D
   		GFS_VARIABLE(cell, f[2]->i) += p->phiforce.z*dist;
 	#endif
}

