#include "trajectory.hpp"


double
trajectory::timeEvolution(particle &a){

	calculateNonDimension(a);
	double timeStep;

	// use analytical or euler?
	flags->analytical=0;
	if(a.Re<0.01 && a.Mach<0.1 && flags->analytical==1) timeStep=analytical(a);
	else timeStep=euler(a);

	return timeStep;
}

double
trajectory::euler(particle &a){
	double timeStep;
	// calculate time step
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		double v2rand=a.Urand.x[0]*a.Urand.x[0]+a.Urand.x[1]*a.Urand.x[1]+a.Urand.x[2]*a.Urand.x[2];
		if(v2rand>v2) vmag=sqrt(v2rand);
		timeStep=meshScale/vmag;		// 1e-6 is scale of cell
		//dt=a.dt;
	}

	// a.tini is remained life time of dispersion
	// if it is smaller than timeStep
	if(a.tini<timeStep){
		timeStep=a.tini;
		a.update=1;
	}
	for(auto &force : forces) force->computeFD(vars,flags,a);
	for(int i=0; i<3; i++) a.v.x[i]+=timeStep*a.F.x[i];
	for(int i=0; i<3; i++) a.x.x[i]+=timeStep*a.v.x[i];

  	if(a.x.x[plane2D]<Axis && flags->dimensionFlag>0) {
		a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
		a.v.x[plane2D]*=-1.0;
		a.F.x[plane2D]*=-1.0;
		a.reflect*=-1;
	}
	a.tini-=timeStep;
	return timeStep;
}


double
trajectory::analytical(particle &a){
	double timeStep;
	if(flags->autoStep) {
		double v2=a.v.x[0]*a.v.x[0]+a.v.x[1]*a.v.x[1]+a.v.x[2]*a.v.x[2];
		double vmag=sqrt(v2);
		timeStep=meshScale/vmag;
	}

	if(a.tini<timeStep){
		timeStep=a.tini;
		a.update=1;
	}

	double Ux=vars->U[a.cell].x[0]+a.Urand.x[0];
	double Uy=vars->U[a.cell].x[1]+a.Urand.x[1];
	double Uz=vars->U[a.cell].x[2]+a.Urand.x[2];
	double dvx=a.v.x[0]-Ux;
	double dvy=a.v.x[1]-Uy;
	double dvz=a.v.x[2]-Uz;

	double EXP=exp(-timeStep*a.beta);
	double dvxEXP=dvx*EXP;
	double dvyEXP=dvy*EXP;
	double dvzEXP=dvz*EXP;

	a.v.x[0]=Ux+dvxEXP;
	a.v.x[1]=Uy+dvyEXP;
	a.v.x[2]=Uz+dvzEXP;

	a.x.x[0]+=Ux*timeStep-dvxEXP/a.beta;
	a.x.x[1]+=Uy*timeStep-dvyEXP/a.beta;
	a.x.x[2]+=Uz*timeStep-dvzEXP/a.beta;

  if(a.x.x[plane2D]<Axis && flags->dimensionFlag>0) {
    a.x.x[plane2D]=2*Axis-a.x.x[plane2D];
    a.v.x[plane2D]*=-1.0;
    a.F.x[plane2D]*=-1.0;
    a.reflect*=-1;
	}
	a.tini-=timeStep;

	return timeStep;
}
