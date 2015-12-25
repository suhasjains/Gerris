#include<stdio.h>
#include<math.h>
	
	double x = 0.1;
	double y = 0.1;
	double up = 0;
	double vp = 0;
 
	double t = 0;
	double dt = 0.005;
	double end = 2;
	
	double acc_x = 0;
	double acc_y = 0;

	double volume = 0.001; 
	double radius = pow(volume/M_PI, 0.5);
	double visc = 0.00078125; 
	double rho = 1.; 
	double rhop = 2.; 
	double drag_force_x = 0;
	double drag_force_y = 0;
	double lift_force_x = 0;
	double lift_force_y = 0;
	double inertial_force_x = 0;
	double inertial_force_y = 0;	
	double amf_force_x = 0;
	double amf_force_y = 0;
	double buoy_force_x = 0;
	double buoy_force_y = 0;


void compute_drag_force(double x2, double y2, double up2, double vp2) {


	//fluid velocity
	double uf = -sin(M_PI*y2)*cos(M_PI*x2);
	double vf = sin(M_PI*x2)*cos(M_PI*y2);

	//relative velocity
	double rel_u = up2 - uf;
	double rel_v = vp2 - vf;	


			
	double mag_rel = sqrt(rel_u*rel_u + rel_v*rel_v);		

	double Re = 2.*radius*mag_rel*rho/visc;


	double cd;

	if(Re < 1e-8){
                drag_force_x = 0.;
                drag_force_y = 0.;
                return;
        }
        else if(Re < 50.0)
                cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
        else
                cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;

        drag_force_x = -3./(8.*radius)*cd*mag_rel*rel_u*rho;
        drag_force_y = -3./(8.*radius)*cd*mag_rel*rel_v*rho;

}	

void compute_lift_force(double x2, double y2, double up2, double vp2) {

	//fluid velocity
	double uf = -sin(M_PI*y2)*cos(M_PI*x2);
	double vf = sin(M_PI*x2)*cos(M_PI*y2);

	//relative velocity
	double rel_u = up2 - uf;
	double rel_v = vp2 - vf;	

	
	double vort = 2.*M_PI*cos(M_PI*x2)*cos(M_PI*y2);

	double cl = 0.5;
		
	lift_force_x = -rho*cl*rel_v*vort;
	lift_force_y = rho*cl*rel_u*vort;


}


void compute_inertial_force(double x2, double y2) {


	inertial_force_x = -rho*M_PI*sin(2.*M_PI*x2)/2.;
	inertial_force_y = -rho*M_PI*sin(2.*M_PI*y2)/2.;

		
}

void compute_amf_force(double x2, double y2) {

	double cm = 0.5;

	amf_force_x = -rho*M_PI*sin(2.*M_PI*x2)/2.;
	amf_force_y = -rho*M_PI*sin(2.*M_PI*y2)/2.;

	amf_force_x *= cm;
	amf_force_y *= cm;	
}

void compute_buoy_force() {

	double gx = 0;
	double gy = -9.81;

	buoy_force_x = (rhop - rho)*gx;
	buoy_force_y = (rhop - rho)*gy;	

}

void get_acceleration(double x1, double y1, double up1, double vp1) {

	double cm = 0.5;

	
	compute_drag_force(x1,y1,up1,vp1);
	compute_lift_force(x1,y1,up1,vp1);
	compute_inertial_force(x1,y1);
	compute_amf_force(x1,y1);
//	compute_buoy_force();

	acc_x = drag_force_x + lift_force_x + inertial_force_x + amf_force_x + buoy_force_x;
	acc_y = drag_force_y + lift_force_y + inertial_force_y + amf_force_y + buoy_force_y;

	acc_x /= (rhop + rho*cm);
	acc_y /= (rhop + rho*cm);	

}


void advect() {


	double k1_x, k1_y, k2_x, k2_y, k3_x, k3_y, k4_x, k4_y;
	double m1_x, m1_y, m2_x, m2_y, m3_x, m3_y, m4_x, m4_y;
	
	m1_x = up;
	m1_y = vp;
	
	get_acceleration(x,y,m1_x,m1_y);
	
	k1_x = acc_x;
	k1_y = acc_y;

	m2_x = up + k1_x*dt/2.; 
	m2_y = vp + k1_y*dt/2.;

	get_acceleration(x+m1_x*dt/2.,y+m1_y*dt/2.,m2_x,m2_y);

	k2_x = acc_x;
	k2_y = acc_y;

	m3_x = up + k2_x*dt/2.; 
	m3_y = vp + k2_y*dt/2.;

	get_acceleration(x+m2_x*dt/2.,y+m2_y*dt/2.,m3_x,m3_y);

	k3_x = acc_x;
	k3_y = acc_y;
	
	m4_x = up + k3_x*dt/2.; 
	m4_y = vp + k3_y*dt/2.;
	
	get_acceleration(x+m3_x*dt/2.,y+m3_y*dt/2.,m4_x,m4_y);

	k4_x = acc_x;
	k4_y = acc_y;

	up += (k1_x/6. + k2_x/3. + k3_x/3. + k4_x/6.)*dt;	
	vp += (k1_y/6. + k2_y/3. + k3_y/3. + k4_y/6.)*dt;	

	x += (m1_x/6. + m2_x/3. + m3_x/3. + m4_x/6.)*dt;	
	y += (m1_y/6. + m2_y/3. + m3_y/3. + m4_y/6.)*dt;	
		

}


int main() {

	
	
	while(t<end) {
		
		advect();

		printf("1 %g %g %g %g\n", dt, t, x, y);

		t = t + dt;
	}
}

 
