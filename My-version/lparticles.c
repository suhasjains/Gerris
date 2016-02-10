#include "lparticles.h"
#include "refine.h"
#include "adaptive.h"
#include "config.h"
#include "solid.h"

#ifdef HAVE_MPI
#include "mpi_boundary.h"
#endif /*HAVE_MPI*/




/* LParticles: Object */

//initialize coupling force
static void couple_force_init(FttCell *c, GfsVariable * f)
{
//  GFS_VARIABLE(c, f->i) = 0.;
	GFS_VALUE(c, f) = 0.;
}

//resetting coupling force
static void reset_couple_force (GfsDomain *domain, GfsVariable **f)
{
  FttComponent c;
  FttDirection d;
  for(c = 0; c < FTT_DIMENSION; c++){
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                              (FttCellTraverseFunc)couple_force_init,
                              f[c]);
  }
}



/*Create force variable for the domain*/
static void couple_force_define(GfsDomain * d, GfsVariable **force, GString * s)
{

  force[0] = gfs_domain_get_or_add_variable(d, g_strconcat(s->str, "x", NULL),
                                            "x-component of the Lagrangian Particle Force");
  force[1] = gfs_domain_get_or_add_variable(d, g_strconcat(s->str, "y", NULL),
                                            "y-component of the Lagrangian Particle Force");
#if !FTT_2D
  force[2] = gfs_domain_get_or_add_variable(d, g_strconcat(s->str, "z", NULL),
                                            "z-component of the Lagrangian Particle Force");
#endif
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

  GFS_VALUE(cell, f[0]) += p->phiforce.x*dist;
  GFS_VALUE(cell, f[1]) += p->phiforce.y*dist;
#if !FTT_2D
  GFS_VALUE(cell, f[2]) += p->phiforce.z*dist;
#endif
}


static gboolean check_stencil(FttCell * cell, FttVector pos0, gdouble sigma, FttDirection *d)
{
  	if(!cell) return FALSE;

  	gdouble size = ftt_cell_size(cell);
  FttVector pos1;
  ftt_cell_pos(cell, &pos1);

  FttDirection d0;
  gdouble dist, dist1;
  gboolean check = TRUE;
  for(d0 = 0; d0 < FTT_DIMENSION; d0++){
    dist = abs(((&pos1.x)[d0] - size/2.0 - (&pos0.x)[d0]));
    dist1 = abs(((&pos1.x)[d0] + size/2.0 - (&pos0.x)[d0]));

    if(dist < dist1)
      d[d0] = 2.0*d0;
    else
      d[d0] = 2.0*d0 + 1;

    if(dist > sigma || dist1 > sigma)
      check = FALSE;
  }

  return check;
}


/*Gaussian Smoothed Two-way Coupling Force: Applied only on 3Sigma surroundings*/
static void compute_coupling_force (Particle *p, GfsVariable **f)
{
  	if(!p->cell) return;

  	gpointer data[3];
  	data[0] = p;
  	data[1] = f;

  	gdouble size = ftt_cell_size(p->cell);
  	gdouble radius = pow(p->volume/M_PI,1./2.);
	#if !FTT_2D
  		radius = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
	#endif

  	gdouble sigma = MAX(2.*radius, size)/2.;
  	data[2] = &sigma;

  	FttCell *cell = p->cell, *neighbor1, *neighbor2;
  	FttDirection d[FTT_DIMENSION];
  	while(!FTT_CELL_IS_ROOT(cell) && ftt_cell_parent(cell) && check_stencil(cell, p->pos, 3.*sigma, d)){
    		cell = ftt_cell_parent(cell);
  	}

	/*   if(FTT_CELL_IS_ROOT(cell)) */
	/*     return; */

  	g_assert(cell!=NULL);

  	ftt_cell_traverse (cell, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                     (FttCellTraverseFunc) compute_smooth_force, data);

  	FttDirection d0, d1;
  	for(d0 = 0; d0 < FTT_DIMENSION; d0++){
    		neighbor1 = ftt_cell_neighbor(cell, d[d0]);
    		if(neighbor1)
      			ftt_cell_traverse (neighbor1, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                         (FttCellTraverseFunc) compute_smooth_force, data);
    	else
      		continue;
  	/*   Do for neighbor if it exists else continue */
    		for(d1 = d0 + 1; d1 < FTT_DIMENSION; d1++){
      			neighbor2 = ftt_cell_neighbor(neighbor1, d[d1]);
      			if(neighbor2){
        		ftt_cell_traverse (neighbor2, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                           (FttCellTraverseFunc) compute_smooth_force, data);
      			}
   		/*    Do for neighbor if it exists else continue */
    		}
  	}
	
	#if !FTT_2D
  	/*Corner Cell*/
  		neighbor1 = ftt_cell_neighbor(cell, d[0]);
  		if(neighbor1)
    			neighbor1 = ftt_cell_neighbor(neighbor1, d[1]);
  		if(neighbor1)
    			neighbor1 = ftt_cell_neighbor(neighbor1, d[2]);
  		if(neighbor1)
    			ftt_cell_traverse (neighbor1, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,(FttCellTraverseFunc) compute_smooth_force, data);
	#endif

}


/*Un-smoothed Two-way Coupling Force*/
static void compute_coupling_force3 (Particle *p, GfsVariable **f)
{

  if(p->cell == NULL) return;
  FttVector force = p->phiforce;

//  GFS_VARIABLE(p->cell, f[0]->i) += force.x;
//  GFS_VARIABLE(p->cell, f[1]->i) += force.y;
//#if !FTT_2D
//  GFS_VARIABLE(p->cell, f[2]->i) += force.z;
//#endif

  GFS_VALUE(p->cell, f[0]) += force.x;
  GFS_VALUE(p->cell, f[1]) += force.y;
#if !FTT_2D
  GFS_VALUE(p->cell, f[2]) += force.z;
#endif
	
}


/*Gaussian Smoothed Two-way Coupling Force: Applied on All cells in the domain*/
static void compute_coupling_force2 (Particle *p, GfsVariable **f, GfsDomain *domain)
{
  if(!p->cell) return;

  gpointer data[3];
  data[0] = p;
  data[1] = f;

  gdouble size = ftt_cell_size(p->cell);
  gdouble radius = pow(p->volume/M_PI,1./2.);
#if !FTT_2D
  radius = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
#endif

  gdouble sigma = MAX(radius, size);
  data[2] = &sigma;

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                            (FttCellTraverseFunc) compute_smooth_force, data);

}



/*Updating the particle velocity due to forces acting on the particle*/
static void compute_particle_velocity (Particle * p, double dt) {

        p->vel.x +=  dt* p->acc.x;
        p->vel.y +=  dt* p->acc.y;
        #if !FTT_2D
                p->vel.z +=  dt* p->acc.z;
        #endif
}

/*Updating the particle acceleration*/
static void compute_particle_acceleration(Particle *p, FttVector force,gdouble cm, gdouble fluid_rho) {
        
        //Taking a component of the added mass force to the LHS of the momentum equation of the particle
        p->acc.x += force.x/(p->density + fluid_rho*cm);
        p->acc.y += force.y/(p->density + fluid_rho*cm);
        #if !FTT_2D
                p->acc.z += force.z/(p->density + fluid_rho*cm);
        #endif
}

/*Same as in source.c used here to obtained viscosity*/
static GfsSourceDiffusion * source_diffusion_viscosity (GfsVariable * v)
{
  	if (v->sources) {
    		GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    		while (i) {
      			GtsObject * o = i->data;
	
      			if (GFS_IS_SOURCE_DIFFUSION (o))
        			return GFS_SOURCE_DIFFUSION (o);
      			i = i->next;
    		}
  	}
  	return NULL;
}


/*Calculation of vorticity vector*/
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


/*Subtraction of two vectors*/
static void subs_fttvectors (FttVector *a, FttVector *b, FttVector *result) {

        result->x  = a->x - b->x;
        result->y  = a->y - b->y;
        result->z  = a->z - b->z;
}


/*Calculation of the lift force acting on the unit volume of the particle*/
static void compute_lift_force(ForceParams *pars) {


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

/*Calculation of drag force acting on the unit volume of the particle */
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

        gdouble radius;
        radius = pow(p->volume/M_PI, 1./2.);
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

	//printf("Re=%g\n",Re);
	if(fcoeffs->cdrag){

        ///GFS_VARIABLE has only 1 arguments - I dont understand this
//        GFS_VARIABLE(p->cell, pars->lagrangian->reynolds->i) = Re;

//        GFS_VARIABLE(p->cell, pars->lagrangian->urel->i) = relative_vel.x;

//        GFS_VARIABLE(p->cell, pars->lagrangian->vrel->i) = relative_vel.y;

        #if !FTT_2D
//        GFS_VARIABLE(p->cell, pars->lagrangian->wrel->i) = relative_vel.z;
        #endif

//        GFS_VARIABLE(p->cell, pars->lagrangian->pdia->i) = 2.0*radius;

        fcoeffs->cd = gfs_function_value (fcoeffs->cdrag, p->cell);
        force->x = fcoeffs->cd*relative_vel.x*fluid_rho;
        force->y = fcoeffs->cd*relative_vel.y*fluid_rho;
//        GFS_VARIABLE(p->cell, pars->lagrangian->pdia->i) = 2.0*radius;

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
	//printf("cd=%g\n",fcoeffs->cd);

        force->x = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.x*fluid_rho;
        force->y = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.y*fluid_rho;
        force->z = 3./(8.*radius)*fcoeffs->cd*norm_relative_vel*relative_vel.z*fluid_rho;

}



/*Calculation of the inertial forces on the unit volume of the particle*/
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

        /*Calculation of local derivative*/
        if(dt > 0.){
                force->x = fluid_rho*(fluid_vel.x - fluid_veln.x)/dt;
                force->y = fluid_rho*(fluid_vel.y - fluid_veln.y)/dt;
        #if !FTT_2D
                force->z = fluid_rho*(fluid_vel.z - fluid_veln.z)/dt;
        #endif
        }


        /*Appending Convective derivative to local derivative to get the total force*/
        FttComponent c;
        for(c=0;c < FTT_DIMENSION;c++){
                force->x += fluid_rho*gfs_center_gradient(p->cell, c, u[0]->i)*GFS_VALUE(p->cell, u[c])/size;
                force->y += fluid_rho*gfs_center_gradient(p->cell, c, u[1]->i)*GFS_VALUE(p->cell, u[c])/size;
        #if !FTT_2D
                force->z += fluid_rho*gfs_center_gradient(p->cell, c, u[2]->i)*GFS_VALUE(p->cell, u[c])/size;
	#endif
        }
}


/*Calculation of the added mass force acting on the unit volume of the particle*/
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

        /*Calculation of local derivative*/
        if(dt > 0.){
                force->x = fluid_rho*(fluid_vel.x - fluid_veln.x)/dt;
                force->y = fluid_rho*(fluid_vel.y - fluid_veln.y)/dt;
        #if !FTT_2D
                force->z = fluid_rho*(fluid_vel.z - fluid_veln.z)/dt;
        #endif
        }


	/*Appending Convective derivative to local derivative to get the total force*/
        FttComponent c;
        for(c=0;c < FTT_DIMENSION;c++){
                force->x += fluid_rho*gfs_center_gradient(p->cell, c, u[0]->i)*GFS_VALUE(p->cell, u[c])/size;
                force->y += fluid_rho*gfs_center_gradient(p->cell, c, u[1]->i)*GFS_VALUE(p->cell, u[c])/size;
        #if !FTT_2D
                force->z += fluid_rho*gfs_center_gradient(p->cell, c, u[2]->i)*GFS_VALUE(p->cell, u[c])/size;
        #endif

   	}


	//coefficient
	fcoeffs->cm = 0.5;

        //Coefficient
        force->x = fcoeffs->cm*force->x;
        force->y = fcoeffs->cm*force->y;
        force->z = fcoeffs->cm*force->z;
}


/*Calculation of the buoyant forces acting on the unit volume of the particle*/
static void compute_buoyant_force(ForceParams * pars) {

        Particle *p = pars->p;
        GfsVariable **u = pars->u;
        gdouble fluid_rho = pars->fluid_rho;
        FttVector * force = pars->force;
        gdouble g[3];
        FttComponent c;

        //getting gravity directly from the sources
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

        ///p->density-fluid_rho is not defined
        force->x = (p->density - fluid_rho)*g[0];
        force->y = (p->density - fluid_rho)*g[1];
	#if !FTT_2D
  		force->z = (p->density-fluid_rho)*g[2];
	#endif
}


static void particles_force(ForceParams * pars)
{
  Particle *p = pars->p;
  FttVector * force = pars->force;
  gdouble fluid_rho = pars->fluid_rho;
  ForceCoefficients * fcoeffs = pars->fcoeffs;

  if(p->cell==NULL)
    return;


/*   gdouble radius = pow(p->volume/M_PI,1./2.); */
/* #if FTT_3D */
/*   radius = pow(3.0*(p->volume)/4.0/M_PI, 1./3.); */
/* #endif */
/*   gdouble p3dvolume = 4.0/3.0*M_PI*radius*radius*radius; */
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


  p->phiforce.x /=(pars->fluid_rho);
  p->phiforce.y /=(pars->fluid_rho);
  p->phiforce.z /=(pars->fluid_rho);

	//for test case
	gdouble viscosity = pars->viscosity;
	p->phiforce.x = 50.0*2.0*M_PI*pow(viscosity,2.0)*fluid_rho;
	p->phiforce.y = 50.0*2.0*M_PI*pow(viscosity,2.0)*fluid_rho;
	p->phiforce.z = 50.0*2.0*M_PI*pow(viscosity,2.0)*fluid_rho;


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

static void assign_val_funcs (GfsFunction * f, GtsFile *fp, LParticles *o) {

  	GfsDomain * domain = GFS_DOMAIN(gfs_object_simulation (o));


  	gts_file_next_token (fp);
  	if (fp->type != '=') {
    		gts_file_error (fp, "expecting `='");
    		return;
  	}

  	gts_file_next_token (fp);
  	gfs_function_read (f, domain, fp);

  	if (fp->type == GTS_ERROR) {
    		gts_file_error (fp, "Trouble reading function for Force Coeffs in GfsLagrangianParticles module");
    		gts_object_destroy (GTS_OBJECT (f));
  	  	return;
  	}
}

static void previous_time_vel(GfsDomain *d, GfsVariable **un) {

  	un[0] = gfs_domain_get_or_add_variable(d, "Un", "x-component of Velocity at previous time step");
  	un[1] = gfs_domain_get_or_add_variable(d, "Vn", "y-component of Velocity at previous time step");
	#if !FTT_2D
  		un[2] = gfs_domain_get_or_add_variable(d, "Wn", "z-component of Velocity at previous time step");
	#endif
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
    		
		lagrangian->couplingforce = g_malloc(sizeof(GfsVariable));
     		couple_force_define(domain, lagrangian->couplingforce, lagrangian->name);


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
	do
		gts_file_next_token (fp);	
	while (fp->type == '\n');
	
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
			
			else if(g_ascii_strcasecmp(fp->token->str, "lift") == 0)
			        assign_val_vars (&lagrangian->fcoeff.lift, fp, *o);

      			else if(g_ascii_strcasecmp(fp->token->str, "drag") == 0)
        			assign_val_vars (&lagrangian->fcoeff.drag, fp, *o);

      			else if(g_ascii_strcasecmp(fp->token->str, "amf") == 0) {
			        assign_val_vars (&lagrangian->fcoeff.amf, fp, *o);}
			
			else if(g_ascii_strcasecmp(fp->token->str, "buoy") == 0)
		        	assign_val_vars (&lagrangian->fcoeff.buoy, fp, *o);

			else if(g_ascii_strcasecmp(fp->token->str, "cdrag") == 0){
			        lagrangian->fcoeff.cdrag = gfs_function_new (gfs_function_class (), 0.);
        			assign_val_funcs (lagrangian->fcoeff.cdrag, fp, lagrangian);
      			}

      			else if(g_ascii_strcasecmp(fp->token->str, "clift") == 0){
        			lagrangian->fcoeff.clift = gfs_function_new (gfs_function_class (), 0.);
        			assign_val_funcs (lagrangian->fcoeff.clift, fp, lagrangian);
      			}
		
			else if(g_ascii_strcasecmp(fp->token->str, "camf") == 0){
			        lagrangian->fcoeff.camf = gfs_function_new (gfs_function_class (), 0.);
        			assign_val_funcs (lagrangian->fcoeff.camf, fp, lagrangian);
      			}

			else if(g_ascii_strcasecmp(fp->token->str, "inertial") == 0)
        			assign_val_vars (&lagrangian->fcoeff.inertial, fp, *o);

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

	lagrangian->forces = NULL;
  	force_pointer f;

  	if(lagrangian->fcoeff.lift == 1){
    		f = &compute_lift_force;
    		lagrangian->forces = g_slist_append(lagrangian->forces, f);
  	}

  	if(lagrangian->fcoeff.drag == 1){
    		f = &compute_drag_force;
    		lagrangian->forces = g_slist_append(lagrangian->forces, f);
  	}

  	if(lagrangian->fcoeff.inertial == 1){
    		f = &compute_inertial_force;
    		lagrangian->forces = g_slist_append(lagrangian->forces, f);
  	}

	if(lagrangian->fcoeff.amf == 1){
		f = &compute_amf_force;
    		lagrangian->forces = g_slist_append(lagrangian->forces, f);
  	}

	if(lagrangian->fcoeff.buoy == 1){
    		f = &compute_buoyant_force;
    		lagrangian->forces = g_slist_append(lagrangian->forces, f);
  	}



	lagrangian->un = g_malloc(sizeof(GfsVariable));
  	previous_time_vel(domain, lagrangian->un);


 	/* do not forget to prepare for next read */
  	gts_file_next_token (fp);
}

static void l_particles_write (GtsObject * o, FILE * fp)
{
  	/* call write method of parent */
  	if (GTS_OBJECT_CLASS (l_particles_class ())->parent_class->write) (* GTS_OBJECT_CLASS (l_particles_class ())->parent_class->write) (o, fp);

  	/* do object specific write here */
	LParticles * lagrangian = L_PARTICLES(o);

	//printf("Hi\n");
	
	fprintf(fp," %s ",lagrangian->name->str);
        if(lagrangian->density)
        	fprintf(fp," %s ",lagrangian->density->name);

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

	
	fprintf(fp," {");
	if(lagrangian->fcoeff.init == 1)
                fprintf(fp, " init = 1");
        if(lagrangian->fcoeff.fluidadv == 1)
                fprintf(fp, " fluidadv = 1");
	if(lagrangian->fcoeff.RK4 == 1)
		fprintf(fp, " RK = 1 ");
	if(lagrangian->fcoeff.lift == 1)
             	fprintf(fp, " lift = 1");
     	if(lagrangian->fcoeff.drag == 1)
             	fprintf(fp, " drag = 1");
	if(lagrangian->fcoeff.inertial == 1)
             fprintf(fp, " inertial = 1");
	if(lagrangian->fcoeff.amf == 1)
                fprintf(fp, " amf = 1");
     	if(lagrangian->fcoeff.buoy == 1)
                 fprintf(fp, " buoy = 1");
	if(lagrangian->fcoeff.cdrag){
                fprintf(fp, " cdrag = ");
        	gfs_function_write(lagrangian->fcoeff.cdrag,fp);
        }
        if(lagrangian->fcoeff.clift){
                fprintf(fp, " clift = ");
                gfs_function_write(lagrangian->fcoeff.clift,fp);
        }
	if(lagrangian->fcoeff.camf){
                fprintf(fp, " camf = ");
                gfs_function_write(lagrangian->fcoeff.camf,fp);
        }
	fprintf (fp," } \n");

}



/*l particle object destroy */
static void l_particles_destroy (GtsObject * object) {

        /* do not forget to call destroy method of the parent */
        (* GTS_OBJECT_CLASS (l_particles_class ())->parent_class->destroy) ( object);

        /* do object-specific cleanup here */
        LParticles * lagrangian = L_PARTICLES(object);


        if(lagrangian->particles)
                g_slist_foreach (lagrangian->particles, (GFunc) g_free, NULL);

        g_slist_free(lagrangian->particles);

        g_string_free(lagrangian->name, TRUE);
        
	if(lagrangian->fcoeff.cdrag)
    		g_free(lagrangian->fcoeff.cdrag);

  	if(lagrangian->fcoeff.clift)
    		g_free(lagrangian->fcoeff.clift);
  
//	if(lagrangian->fcoeff.camf)
//		g_free(lagrangian->fcoeff.camf);

	if(lagrangian->couplingforce)
   		g_free(lagrangian->couplingforce);



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

/*Copying value for previous time*/
static void copy_cell_gfsvariable(FttCell *c, gpointer * data) {
        GfsVariable * v1 = (GfsVariable *) data[0];
        GfsVariable * v2 = (GfsVariable *) data[1];

        GFS_VALUE(c, v1) = GFS_VALUE(c, v2);

}



/*Storing previous time step value*/
static void store_domain_previous_vel(GfsDomain *d, GfsVariable **un)
{
  	GfsVariable ** u = gfs_domain_velocity (d);
  	FttComponent c;
  	for(c = 0; c < FTT_DIMENSION; c++){
    		gpointer data[2];

    		data[0] = un[c];
    		data[1] = u[c];
    		gfs_domain_cell_traverse (d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                              (FttCellTraverseFunc)copy_cell_gfsvariable,
                              data);
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


		store_domain_previous_vel(domain, lagrangian->un);
                }
		
		//Fetch simulation parameters
                GSList *i = lagrangian->particles;
                Particle * p;
		FttVector force;
 	       	force_pointer f;
        	GSList *fs;


                //Fetch force parameters
                ForceParams *pars = g_malloc(sizeof(ForceParams));
		gdouble dt = pars->dt = sim->advection_params.dt;
		//gdouble dtn = pars->dtn;
                pars->fcoeffs = &lagrangian->fcoeff;
                GfsVariable ** u = pars->u = gfs_domain_velocity (domain);
		pars->un = lagrangian->un;
		//GfsVariable ** un;
		//if(lagrangian->first_call) 	
		//	un = lagrangian->pars.un = u; 	
		//else 	
		//	un = lagrangian->pars.un;	
        	GfsSourceDiffusion *d = source_diffusion_viscosity(u[0]);
		pars->force = &force;
		pars->lagrangian = lagrangian;



		reset_couple_force (domain, lagrangian->couplingforce);

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

         				pars->p = p;
         				pars->fluid_rho = fluid_rho;
         				pars->viscosity = viscosity;


         				fs = lagrangian->forces;
         				while(fs){
           					f = (force_pointer)fs->data;
           					(f)(pars);
           					compute_particle_acceleration(p, force, pars->fcoeffs->cm, fluid_rho);
           					fs = fs->next;
         				}
           			
					//Particle velocity from Newton's equation	
					compute_particle_velocity (p, dt);


					//Computing coupling force
         				if(p->cell != NULL){
           					particles_force(pars);
           					compute_coupling_force(p, lagrangian->couplingforce);
						//compute_coupling_force2(p, lagrangian->couplingforce, domain);
						//compute_coupling_force3(p, lagrangian->couplingforce);
					}
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
				//if(lagrangian->fcoeff.RK4 == 1) {
                                //	advect_particles_RK4(p, dtn, u, un);
					//printf("Doing RK4 advection\n");
                                //}

				//else {
                                	//Advect particle according to particle velocity
                                	advect_particles(p, dt);
					//printf("Doing normal advection\n");
                        	//}
			}
                    	
			//printf("%d %g %g %g %g %g\n",p->id, dt, sim->time.t, p->pos.x, p->pos.y, p->pos.z);
                        i = i->next;
                }

		//Storing to be used next time step
		//lagrangian->pars.un = u;
		//lagrangian->pars.dtn = dt;
                store_domain_previous_vel(domain, lagrangian->un); 
		g_free(pars);
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
	object->forces = NULL;
        object->first_call = TRUE;

	object->fcoeff.cl = 0.;
  	object->fcoeff.cd = 0.;

  	object->fcoeff.lift = 0;
  	object->fcoeff.drag = 0;
  	object->fcoeff.inertial = 0;
	object->fcoeff.amf = 0;
	object->fcoeff.buoy = 0;

  	object->fcoeff.cdrag = NULL;
	object->fcoeff.clift = NULL;
	object->fcoeff.camf = NULL;
	
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



/*GfsSourceGfsLagrangianParticles: Object*/

static void source_lagrangian_destroy (GtsObject * o)
{

  	g_string_free(GFS_SOURCE_LAGRANGIAN(o)->name,TRUE);

  	(* GTS_OBJECT_CLASS (gfs_source_lagrangian_class ())->parent_class->destroy) (o);

}

static void source_lagrangian_read (GtsObject ** o, GtsFile * fp)
{

	//parent class read method
  	if (GTS_OBJECT_CLASS (gfs_source_diffusion_explicit_class ())->parent_class->read) (* GTS_OBJECT_CLASS (gfs_source_lagrangian_class ())->parent_class->read) (o, fp);

  	if (fp->type == GTS_ERROR)
    		return;

	//object read method
  	GfsSourceLagrangian * source = GFS_SOURCE_LAGRANGIAN(*o);

  	if (fp->type == GTS_STRING) {
    		source->name = g_string_new(fp->token->str);
    		gts_file_next_token (fp);
  	}
  	else
    	gts_file_error (fp, "expecting a string");

}



static void source_lagrangian_write (GtsObject * o, FILE * fp)
{
  	(* GTS_OBJECT_CLASS (gfs_source_lagrangian_class ())->parent_class->write) (o, fp);

  	GfsSourceLagrangian * source = GFS_SOURCE_LAGRANGIAN(o);

  	fprintf(fp," %s\n",source->name->str);
}

/* Add source term in the Momentum equations*/
static gdouble compute_phi_force_mac_value(GfsSourceGeneric * s,
                                           FttCell * cell,
                                           GfsVariable * v)
{

  	GfsDomain *d = GFS_DOMAIN(gfs_object_simulation(s));
  	GfsVariable *force;
  	gchar * name;
  	GfsSourceLagrangian *sv = GFS_SOURCE_LAGRANGIAN(s);

  	FttCellFace f;
  	f.cell = cell;

  	FttCellNeighbors n;
  	ftt_cell_neighbors (cell, &n);

  	switch(v->component){
  		case FTT_X:
    			name = g_strconcat(sv->name->str, "x", NULL);
    			force = gfs_variable_from_name(d->variables, name);
    			f.d = FTT_RIGHT;
    			f.neighbor = n.c[f.d];
    		break;
  		case FTT_Y:
    			name = g_strconcat(sv->name->str, "y", NULL);
    			force = gfs_variable_from_name(d->variables, name);
    			f.d = FTT_TOP;
    			f.neighbor = n.c[f.d];
    		break;
		#if !FTT_2D
  			case FTT_Z:
    				name = g_strconcat(sv->name->str, "z", NULL);
    				force = gfs_variable_from_name(d->variables, name);
    				f.d = FTT_FRONT;
    				f.neighbor = n.c[f.d];
    			break;
		#endif
  		default: g_assert_not_reached ();
  	}

  	g_free(name);
  	if(force!=NULL)
    		return  gfs_face_interpolated_value(&f, force->i);
  	else
    	g_assert_not_reached ();

  	return 0.;

}


static gdouble compute_phi_force_centered_value(GfsSourceGeneric * s,
                                                FttCell * cell,
                                                GfsVariable * v)
{

  	GfsDomain *d = GFS_DOMAIN(gfs_object_simulation(s));
  	GfsVariable *force;
  	gchar * name;
  	GfsSourceLagrangian *sv = GFS_SOURCE_LAGRANGIAN(s);

  	switch(v->component){
  	case FTT_X:
    		name = g_strconcat(sv->name->str, "x", NULL);
    		force = gfs_variable_from_name(d->variables, name);
    	break;
  	case FTT_Y:
    		name = g_strconcat(sv->name->str, "y", NULL);
    		force = gfs_variable_from_name(d->variables, name);
    	break;
	#if !FTT_2D
  		case FTT_Z:
    			name = g_strconcat(sv->name->str, "z", NULL);
    			force = gfs_variable_from_name(d->variables, name);
    		break;
	#endif
  	default: g_assert_not_reached ();
  	}

  	g_free(name);
  	if(force!=NULL){
    		return GFS_VALUE(cell,force);}
  	else
    	g_assert_not_reached ();

}


static void gfs_source_lagrangian_class_init (GfsSourceGenericClass * klass)
{

  	GTS_OBJECT_CLASS (klass)->destroy = source_lagrangian_destroy;
  	GTS_OBJECT_CLASS (klass)->read = source_lagrangian_read;
  	GTS_OBJECT_CLASS (klass)->write = source_lagrangian_write;

}

static void gfs_source_lagrangian_init (GfsSourceGeneric * s)
{
  	s->mac_value= compute_phi_force_mac_value;
  	s->centered_value = compute_phi_force_centered_value;
}

GfsSourceGenericClass * gfs_source_lagrangian_class (void)
{
  	static GfsSourceGenericClass * klass = NULL;

  	if (klass == NULL) {
    		GtsObjectClassInfo gfs_source_lagrangian_info = {
      		"GfsSourceLagrangian",
      		sizeof (GfsSourceLagrangian),
      		sizeof (GfsSourceGenericClass),
      		(GtsObjectClassInitFunc) gfs_source_lagrangian_class_init,
      		(GtsObjectInitFunc) gfs_source_lagrangian_init,
      		(GtsArgSetFunc) NULL,
      		(GtsArgGetFunc) NULL
    	};
    	klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
                                  &gfs_source_lagrangian_info);
  	}

  	return klass;
}





//typedef struct {
//  FttVector pos,vel;
//  gdouble vol,density;
//} DropletProp;
//
//typedef struct {
//  GfsVariable * tag, * c, *t;
//  Particle * drops;
//  GfsVariable **u;
//  guint * sizes;
//  guint n, min; //n - number of droplets
//  gdouble resetval;
//  gdouble density;
//  GfsFunction *fc;
//} DropletsPar;
//
//
//
//
///*Droplet to Particles: Object*/
//static int greater (const void * a, const void * b)
//{
//  return *((guint *)a) > *((guint *)b) ? -1 : 1;
//}
//
//
//static void reset_small_fraction (FttCell * cell, DropletsPar * p)
//{
//
//  gint i = GFS_VALUE (cell, p->tag);
//  if (i > 0 && p->sizes[i - 1] < p->min)
//    GFS_VALUE (cell, p->c) = p->resetval;
//}
//
//static void compute_droplet_properties (FttCell * cell, DropletsPar * p)
//{
//  g_return_if_fail (cell != NULL);
//
//  gint i = GFS_VALUE (cell, p->tag);
//  gdouble h = ftt_cell_size(cell),vol,avgfact;
//  FttVector pos; 
//  ftt_cell_pos(cell,&pos);
//  GfsVariable **u = p->u;
//
//  if (i > 0){
//    avgfact = (gdouble) p->sizes[i - 1];
//    p->sizes[i - 1]++;
//
//
//#if FTT_2D
//    vol = h*h*GFS_VALUE(cell,p->t);
//#else
//    vol = h*h*h*GFS_VALUE(cell,p->t);
//#endif
//
//    p->drops[i-1].volume += vol*GFS_VALUE(cell,p->c);
//
//
//    p->drops[i-1].pos.x = (p->drops[i-1].pos.x*avgfact+ pos.x)/p->sizes[i-1];
//    p->drops[i-1].pos.y = (p->drops[i-1].pos.y*avgfact+ pos.y)/p->sizes[i-1];
//    p->drops[i-1].vel.x = (p->drops[i-1].vel.x*avgfact+ GFS_VALUE(cell,u[0]))/p->sizes[i-1];
//    p->drops[i-1].vel.y = (p->drops[i-1].vel.y*avgfact+ GFS_VALUE(cell,u[1]))/p->sizes[i-1];
//#if !FTT_2D
//    p->drops[i-1].pos.z = (p->drops[i-1].pos.z*avgfact+pos.z)/p->sizes[i-1];
//    p->drops[i-1].vel.z = (p->drops[i-1].vel.z*avgfact+GFS_VALUE(cell,u[2]))/p->sizes[i-1];
//#endif
//  }
//
//}
//
//
//static void convert_droplets (GfsDomain * domain, 
//			      DropletsPar *p, LParticles * lagrangian)
//{
//
//  GfsSimulation *sim = gfs_object_simulation(lagrangian); 
//
//  g_return_if_fail (p->c != NULL);
//
//  guint i;
//
//  lagrangian->idlast = lagrangian->maxid;
//
//  for(i = 0; i < p->n; i++){
//
//    if(p->sizes[i] < p->min && p->sizes[i] >4){
//      Particle * drop = g_malloc (sizeof(Particle));
//
//      drop->pos = p->drops[i].pos;
//      drop->vel = p->drops[i].vel;
//      drop->volume = p->drops[i].volume;
//      //drop->move = 1;
//      drop->cell = gfs_domain_locate(domain, drop->pos, -1, NULL);
//  
//      if(drop->cell){
//	drop->id = ++lagrangian->idlast;
//	drop->density = sim->physical_params.alpha ? 1./
//	  gfs_function_value(sim->physical_params.alpha,drop->cell) : 1.;
//	drop->density = p->density;
//	lagrangian->particles = g_slist_append(lagrangian->particles,drop);
//	//lagrangian->n++;
//      }       
//    }
//  }
//
//  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
//			    (FttCellTraverseFunc) reset_small_fraction, p);
//}
//
//static void tag_sort_droplets (GfsDomain * domain,
//			       DropletsPar *p,
//			       gint min)
//{
//  g_return_if_fail (domain != NULL);
//  g_return_if_fail (p->c != NULL);
//
//  p->n = gfs_domain_tag_droplets (domain, p->c, p->tag);
//
//  if (p->n > 0 && -min < (gint) p->n) {
//
//    p->sizes = g_malloc0 (p->n*sizeof (guint));
//
//    p->drops = g_malloc0 (p->n*sizeof (Particle));
//  /*Compute Droplet Properties here*/
//    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
//			      (FttCellTraverseFunc) compute_droplet_properties, p);
//
//    if (min >= 0)
//      p->min = min;
//    else {
//      guint * tmp = g_malloc (p->n*sizeof (guint));
//      memcpy (tmp, p->sizes, p->n*sizeof (guint));
//      qsort (tmp, p->n, sizeof (guint), greater);
//      /* fixme: this won't work for parallel jobs */
//      p->min = tmp[-1 - min];
//      g_free (tmp);
//    }
//
//  }
//}
//
//
//static void droplet_to_particles_destroy (GtsObject * o)
//{
//  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->destroy) (o);
//
//  GfsDropletToParticles * d = DROPLET_TO_PARTICLES(o);
//
//  if (d->fc)
//    gts_object_destroy (GTS_OBJECT (d->fc));
//  
//  if (d->fconvert)
//    gts_object_destroy (GTS_OBJECT (d->fconvert));
//  
//
//  if (d->shape){
//    
//    if (d->shape->s)
//      gts_object_destroy (GTS_OBJECT (d->shape->s));
//    
//    if (d->shape->m)
//      gts_matrix_destroy (d->shape->m); 
//    
//    g_free (d->shape);
//    
//  }
//}
//
//static void droplet_to_particles_read (GtsObject ** o, GtsFile * fp)
//{
//
//  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->read) (o, fp);
//
//  if (fp->type == GTS_ERROR)
//    return;
//
//  while(fp->type == '\n')
//    gts_file_next_token (fp);
//
//  if (fp->type != GTS_STRING) {
//    gts_file_error (fp, "expecting a string (variable)");
//    return;
//  }
//
//  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
//  GfsDropletToParticles * r = DROPLET_TO_PARTICLES(*o);
//
//  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
//    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
//    return;
//  }
//  gts_file_next_token (fp);
//
//  if (fp->type != GTS_INT) {
//    gts_file_error (fp, "expecting an integer (min)");
//    return;
//  }
//
//  r->min = atoi (fp->token->str);
//  gts_file_next_token (fp);
//
//  if (fp->type == GTS_FLOAT) {
//    r->resetwith = atof (fp->token->str);
//    gts_file_next_token (fp);
//  }
//
//  if (fp->type == GTS_FLOAT) {
//    r->density = atof (fp->token->str);
//    gts_file_next_token (fp);
//  }
//
//  if (fp->type != '\n' && fp->type != GTS_INT) {
//    r->fc = gfs_function_new (gfs_function_class (), 0.);
//    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
//  }
//
//  if (fp->type == GTS_INT) {
//    r->convert = atoi (fp->token->str);
//    gts_file_next_token (fp);
//
//    if (fp->type == GTS_INT){
//      r->maxlevel = atoi (fp->token->str);
//      gts_file_next_token (fp);
//    }
//
//    while(fp->type == '\n')
//      gts_file_next_token (fp);
//
//    r->shape = GFS_GENERIC_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
//    gfs_generic_surface_read (r->shape, gfs_object_simulation (*o), fp);
//    
//  
//    if (fp->type != '\n') {
//      r->fconvert = gfs_function_new (gfs_function_class (), 0.);
//      gfs_function_read (r->fconvert, gfs_object_simulation (r), fp);
//    }
//  }
//
//}
//
//static void droplet_to_particles_write (GtsObject * o, FILE * fp)
//{
//
//  (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->write) (o, fp);
//
//  GfsDropletToParticles * r = DROPLET_TO_PARTICLES(o);
//
//  fprintf (fp, " %s %d %lf %lf ", r->c->name, r->min, r->resetwith, r->density);
//
//  if (r->fc){
//    gfs_function_write (r->fc, fp);
//  }
//
//  if(r->convert != 0){
//    fprintf (fp," %d %d ", r->convert, r->maxlevel);
//
//    if (r->shape) {
//      fputs (" ( ", fp);
//      gfs_function_write (r->shape->f, fp);
//      fputs (" ) ", fp);
//    }
//
//    if (r->fconvert){
//      gfs_function_write (r->fconvert, fp);
//    }
//  }
//
//}
//
//static void compute_v (FttCell * cell, GfsDropletToParticles * d)
//{
//  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
//}
//
//void save_p_solid (FttCell * cell, gpointer * data)
//{
//  GfsVariable * c = data[0];
//  GfsVariable * save_prev = data[3];
//
//  //GFS_VARIABLE (cell, save_prev->i) = GFS_VARIABLE (cell, c->i);
//  GFS_VALUE (cell, save_prev) = GFS_VALUE (cell, c);
//  
//	//GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i)) = GFS_STATE (cell)->solid;
//	GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, c)) = GFS_STATE (cell)->solid;
//	
//  GFS_STATE (cell)->solid = NULL;
//}
//
//void restore_p_solid (FttCell * cell, gpointer * data)
//{
//  GfsVariable * c = data[0];
//  gboolean * not_cut = data[1];
//  GfsVariable * status = data[2];
//
//  gdouble * resetwith = data[6];
//
//  GfsSolidVector * solid = GFS_STATE (cell)->solid;
//
//  //GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VARIABLE (cell, c->i));
//  GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, c));
//
//  if (solid) {
//
//    //GFS_VARIABLE (cell, c->i) =  solid->a;
//    GFS_VALUE (cell, c) =  solid->a;
//
//    g_free (solid);
//    *not_cut = FALSE;
//  }
//  //else if (GFS_VARIABLE (cell, status->i) == 0.) {
//  else if (GFS_VALUE (cell, status) == 0.) {
//    /* fixme: this can fail for non-contiguous domains (e.g. non-connected GfsBoxes) */
//    g_assert (*not_cut);
//    //GFS_VARIABLE (cell, c->i) = 0.;
//    GFS_VALUE (cell, c) = 0.;
//  }
//  else {
//    //g_assert (GFS_VARIABLE (cell, status->i) == 1. || GFS_VARIABLE (cell, status->i) == 2.);
//    g_assert (GFS_VALUE (cell, status) == 1. || GFS_VALUE (cell, status) == 2.);
//    //GFS_VARIABLE (cell, c->i) =  GFS_VARIABLE (cell, status->i) - 1.;
//    GFS_VALUE (cell, c) =  GFS_VALUE (cell, status) - 1.;
//  }
//}
//
//void add_to_prev_void (FttCell * cell, gpointer * data)
//{
//  GfsVariable * c = data[0];
//  GfsVariable * save_prev = data[3];
//  Particle *p = (Particle *)data[4];
//  GfsDomain *domain = data[5];
//  GfsVariable **u = gfs_domain_velocity(domain);
//
//  if(GFS_VALUE (cell, c) > 0.){
//    guint j;
//    for(j = 0; j < FTT_DIMENSION; j++)
//      GFS_VALUE(cell, u[j]) = (&p->vel.x)[j];
//  }
//
//  //GFS_VARIABLE (cell, c->i) += GFS_VARIABLE (cell, save_prev->i);
//  GFS_VALUE (cell, c) += GFS_VALUE (cell, save_prev);
//
//  if(GFS_VALUE (cell, c) > 1.0)
//    //GFS_VARIABLE (cell, c->i) = 1.0;
//    GFS_VALUE (cell, c) = 1.0;
//  else if (GFS_VALUE (cell, c) < 0.0)
//      //GFS_VARIABLE (cell, c->i) = 0.0;
//      GFS_VALUE (cell, c) = 0.0;
//}
//
//void gfs_domain_assign_fraction (GfsDomain * domain,
//			       GfsGenericSurface * s,
//			       GfsVariable * c, gdouble resetwith, Particle *p)
//{
//  gboolean not_cut = TRUE;
//  gpointer data[7];
//  GfsVariable * status, *save_prev;
//
//  g_return_if_fail (domain != NULL);
//  g_return_if_fail (s != NULL);
//  g_return_if_fail (c != NULL);
//
//  status = gfs_temporary_variable (domain);
//
//  save_prev = gfs_temporary_variable (domain);
//  data[0] = c;
//  data[1] = &not_cut;
//  data[2] = status;
//  data[3] = save_prev;
//  data[4] = p;
//  data[5] = domain;
//  data[6] = &resetwith;
//
//  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
//			    (FttCellTraverseFunc) save_p_solid, data);
//  GSList * l = g_slist_prepend (NULL, s);
//
//  gfs_domain_init_solid_fractions (domain, l, FALSE, NULL, NULL, status);
//  g_slist_free (l);
//
//  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
//			    (FttCellTraverseFunc) restore_p_solid, data);
//
//  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
//			    (FttCellTraverseFunc) add_to_prev_void, data);
//
//  gts_object_destroy (GTS_OBJECT (status));
//  gts_object_destroy (GTS_OBJECT (save_prev));
//
//}
//
//
//static void generate_surface(GfsSurface *surface, GfsSurface *shape, gdouble rad, Particle *p)
//{
//
//  memcpy(surface, shape, sizeof(GfsSurface));
//
//  surface->translate[0] = p->pos.x;
//  surface->translate[1] = p->pos.y;
//  surface->translate[2] = p->pos.z;
//
//  GtsMatrix * m = gts_matrix_translate (NULL, surface->translate);
//
//  surface->scale[0] = rad;
//  surface->scale[1] = rad;
//  surface->scale[2] = rad;
//
//  GtsMatrix * ms = gts_matrix_scale (NULL, surface->scale);
//  if (surface->m)
//    gts_matrix_destroy (surface->m);
//  surface->m = gts_matrix_product (m, ms);
//  gts_matrix_destroy (m);
//  gts_matrix_destroy (ms);
//      
//  if (surface->s) {
//    gts_surface_foreach_vertex (surface->s, (GtsFunc) gts_point_transform, surface->m);
//    gts_matrix_destroy (surface->m);
//    surface->m = NULL;
//    if (surface->flip)
//      gts_surface_foreach_face (surface->s, (GtsFunc) gts_triangle_revert, NULL);
//  }
//  else {
//    GtsMatrix * i = gts_matrix_inverse (surface->m);
//    gts_matrix_destroy (surface->m);
//    surface->m = i;
//  }
//}
//
//static double compute_xyz (FttCell * cell, GfsDropletToParticles * d)
//{
//  if(d->fconvert)
//    return gfs_function_value (d->fconvert, cell);
//  else
//    return -1.0;
//}
//
//#define PROXIMITY_FRAC 0.25
//
//static void check_cut_cells(FttCell *cell, GfsGenericSurface *s, gpointer *data)
//{
//  GfsVariable *v = (GfsVariable *)data[0];
//  GfsVariable **u = (GfsVariable **)data[1];
//  gboolean * check = (gboolean *)data[2];
//  Particle * p = (Particle *)data[3];  
//  
//  if(!(*check)){
//    if(GFS_VALUE(cell, v) > PROXIMITY_FRAC){
//      *check = TRUE;
//    }
//  } 
//}
//
//static gboolean check_proximity(GfsVariable *v, Particle *p, GfsSurface *surface)
//{ 
//
//  GfsDomain *domain = GFS_DOMAIN(gfs_object_simulation(v));
//
//  gboolean check = FALSE;
//  gpointer data[4];
//  data[0] = v;
//  data[1] = gfs_domain_velocity(domain);
//  data[2] = &check;
//  data[3] = p;
//
//  gfs_domain_traverse_cut (domain,
//			   surface,
//			   FTT_PRE_ORDER,
//			   FTT_TRAVERSE_LEAFS,
//			   (FttCellTraverseCutFunc) check_cut_cells,
//			   data);
//  return check;
//}
//
//void refine_implicit_p_cell (FttCell * cell, RefineCut * p)
//{
//  guint maxlevel = p->maxlevel;
//
//  if (ftt_cell_level (cell) < maxlevel && (gfs_cell_is_cut (cell, p->surface, FALSE, maxlevel))){
//    ftt_cell_refine_single (cell, p->domain->cell_init, p->domain->cell_init_data);
//    p->check = TRUE;
//  }
//}
//
//
//static void convert_particles( GfsDomain *domain, GfsDropletToParticles *d)
//{
//
//  LParticles *lagrangian = L_PARTICLES(d);
//  
//  GSList *i = lagrangian->particles;
//  Particle *p;
//
//  while(i){
//    p = (Particle *) (i->data);
//    i = i->next;
//
//    gdouble rad = pow(p->volume/M_PI,1./2.);
//#if !FTT_2D
//    rad = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
//#endif
//    GfsSurface * surface = GFS_GENERIC_SURFACE (gts_object_new 
//						  (GTS_OBJECT_CLASS (gfs_surface_class ())));
//    generate_surface(surface, d->shape, rad, p);
//
//    gboolean check = FALSE;
//
//
//    if(rad > ftt_cell_size(p->cell)){
//      check = check_proximity(d->c, p, surface);
//      if(check)
//	g_warning("Transforming Particle into Droplet-proximity based\n");
//    }
//    else{
//      if(GFS_VALUE(p->cell, d->c) > PROXIMITY_FRAC)
//	check = TRUE;
//      
//      if(check)
//	g_warning("Transforming Particle into Droplet-proximity based\n");
//    }
// 
//    if(compute_xyz(p->cell, d) > 0 || check){
//
//      RefineCut prefine;
//      prefine.domain = domain;
//      prefine.maxlevel = d->maxlevel;
//      prefine.surface = surface;
//      prefine.check = TRUE;
//
//      while(prefine.check){
//      	prefine.check = FALSE;
//	gfs_domain_cell_traverse (domain,
//				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
//				  (FttCellTraverseFunc) refine_implicit_p_cell, &prefine);
//      }
//
//      gfs_domain_assign_fraction (domain,
//				  surface,
//				  d->c, d->resetwith, p);
//
//      lagrangian->particles = g_slist_remove(lagrangian->particles, p);
//      g_free(p);
//    }
//
//    if (surface->s)
//      gts_object_destroy (GTS_OBJECT (surface->s));
//
//    if (surface->m)
//      gts_matrix_destroy (surface->m);
//  }  
//
//  gfs_domain_reshape (domain, gfs_domain_depth(domain));
//
//}
//
//static gboolean droplet_to_particles_event (GfsEvent * event, GfsSimulation * sim)
//{
//
//
//  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class)->event) 
//      (event, sim)) {
//
//    GfsDomain * domain = GFS_DOMAIN (sim);
//    LParticles * lagrangian = L_PARTICLES(event);
//    GfsDropletToParticles *d = DROPLET_TO_PARTICLES(event);
//
//    d->v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
// 
//    DropletsPar p ;
//
//    p.resetval = d->resetwith;
//    p.tag = gfs_temporary_variable (domain);
//    p.u = gfs_domain_velocity(domain);
//    p.density = d->density;
//    p.t = d->c;
//
//    if (d->v){
//      p.c = d->v;
//      tag_sort_droplets (domain, &p, d->min);
//
//      if (p.n > 0 && -d->min < (gint) p.n){
//	p.c = d->c;
//	convert_droplets (domain, &p, lagrangian);
//      }
//    }
//    else {
//
//      d->v = gfs_temporary_variable (domain);
//
//      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
//				(FttCellTraverseFunc) compute_v, d);
//
//      p.c = d->v;
// 
//      tag_sort_droplets (domain, &p, d->min);
//
//       if (p.n > 0 && -d->min < (gint) p.n){
//	p.c = d->c;
//
//	convert_droplets (domain, &p, lagrangian);
//
//      }     
// 
//       gts_object_destroy (GTS_OBJECT (d->v));
//
//    }
//
//    if(d->convert != 0 )
//      convert_particles(domain, d);
//    
//#ifdef HAVE_MPI
//    //mpi_particle_numbering(domain, lagrangian);
//#endif  /*HAVE_MPI*/
//    gts_object_destroy(GTS_OBJECT(p.tag));
//
//    if(p.n > 0){
//      g_free(p.drops);
//      g_free(p.sizes);	
//    }
//
//    return TRUE;
//  }
//  return FALSE;
//}
//
//static void droplet_to_particles_class_init (GfsDropletToParticlesClass * klass)
//{
//  /* define new methods and overload inherited methods here */
//  
//  GFS_EVENT_CLASS (klass)->event = droplet_to_particles_event;
//  GTS_OBJECT_CLASS (klass)->read = droplet_to_particles_read;
//  GTS_OBJECT_CLASS (klass)->write = droplet_to_particles_write;
//  GTS_OBJECT_CLASS (klass)->destroy = droplet_to_particles_destroy;
//}
//
//static void droplet_to_particles_init (GfsDropletToParticles * r)
//{
//  r->resetwith = 0.;
//  r->convert = 0;
//  r->maxlevel = 8;
//
//  r->min = 0;
//}
//
//GfsDropletToParticlesClass * droplet_to_particles_class (void)
//{
//  static GfsDropletToParticlesClass * klass = NULL;
//
//  if (klass == NULL) {
//    GtsObjectClassInfo droplet_to_particles_info = {
//      "GfsDropletToParticles",
//      sizeof (GfsDropletToParticles),
//      sizeof (GfsDropletToParticlesClass),
//      (GtsObjectClassInitFunc) droplet_to_particles_class_init,
//      (GtsObjectInitFunc) droplet_to_particles_init,
//      (GtsArgSetFunc) NULL,
//      (GtsArgGetFunc) NULL
//    };
//    klass = gts_object_class_new (GTS_OBJECT_CLASS (l_particles_class ()),
//				  &droplet_to_particles_info);
//  }
//
//  return klass;
//}
//


/* DropletToParticles: Object */

static void droplet_to_particles_read (GtsObject ** o, GtsFile * fp)
{
  /* call read method of parent */
  if (GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  /* do object-specific read here */

  /* do not forget to prepare for next read */
  gts_file_next_token (fp);
}

static void droplet_to_particles_write (GtsObject * o, FILE * fp)
{
  /* call write method of parent */
  if (GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class->write) 
      (o, fp);

  /* do object specific write here */
}

static gboolean droplet_to_particles_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (droplet_to_particles_class ())->parent_class)->event) (event, sim)) {
    /* do object-specific event here */

    return TRUE;
  }
  return FALSE;
}

static void droplet_to_particles_class_init (DropletToParticlesClass * klass)
{
  /* define new methods and overload inherited methods here */

  GFS_EVENT_CLASS (klass)->event = droplet_to_particles_event;
  GTS_OBJECT_CLASS (klass)->read = droplet_to_particles_read;
  GTS_OBJECT_CLASS (klass)->write = droplet_to_particles_write;
}

static void droplet_to_particles_init (DropletToParticles * object)
{
  /* initialize object here */
}

DropletToParticlesClass * droplet_to_particles_class (void)
{
  static DropletToParticlesClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo droplet_to_particles_info = {
      "DropletToParticles",
      sizeof (DropletToParticles),
      sizeof (DropletToParticlesClass),
      (GtsObjectClassInitFunc) droplet_to_particles_class_init,
      (GtsObjectInitFunc) droplet_to_particles_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (l_particles_class ()),
				  &droplet_to_particles_info);
  }

  return klass;
}

DropletToParticles * droplet_to_particles_new (DropletToParticlesClass * klass)
{
  DropletToParticles * object;

  object = DROPLET_TO_PARTICLES (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}


//\* Initialize modules *\/ */
const gchar gfs_module_name[] = "lparticles";
const gchar * g_module_check_init (void); 

const gchar * g_module_check_init (void) 
{ 
  	l_particles_class (); 
  	gfs_source_lagrangian_class (); 
//	droplet_to_particles_class ();
	droplet_to_particles_class ();
  	return NULL; 
}  

