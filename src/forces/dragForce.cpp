#include "dragForce.hpp"

// Functions for drag force calculation

void
dragForceSM::computeFD(Variables *vars, Flags *flags, particle &par){
	double Cd=computeCd(par);	// Here, this code use another equations for each drag models

	double FD=M_PI*0.125*vars->myu[par.cell]*par.dp*par.Re*Cd/par.m;

	double dUx=vars->U[par.cell].x[0]+par.Urand.x[0]-par.v.x[0];
	double dUy=vars->U[par.cell].x[1]+par.Urand.x[1]-par.v.x[1];
	double dUz=vars->U[par.cell].x[2]+par.Urand.x[2]-par.v.x[2];
	par.F.x[0]=FD*dUx;
	par.F.x[1]=FD*dUy;
	par.F.x[2]=FD*dUz;
}

double
dragForceSM::computeCd(particle &par){
	return 24/par.Re/par.Cc;
}

double
dragForceMA::computeCd(particle &par){
	double CdRE;
	double Re_inv=1/par.Re;
	double Re_inv2=Re_inv*Re_inv;
	if (par.Re < 0.1){CdRE = 24.0 * Re_inv;}
	else if (par.Re < 1.0){	CdRE = 22.73 * Re_inv + 0.0903 * Re_inv2 + 3.69;}
	else if (par.Re < 10.0){CdRE = 29.1667 * Re_inv - 3.8889 * Re_inv2 + 1.222;}
	else if (par.Re < 100.0){CdRE = 46.5 * Re_inv - 116.67 * Re_inv2 + 0.6167;}
	else if (par.Re < 1000.0){CdRE = 98.33 * Re_inv - 2778.0 * Re_inv2 + 0.3644;}
	else if (par.Re < 5000.0){CdRE = 148.62 * Re_inv - 47500.0 * Re_inv2 + 0.357;}
	else if (par.Re < 10000.0){CdRE = -490.546 * Re_inv + 578700.0 * Re_inv2 + 0.46;}
	else if (par.Re < 50000.0){CdRE = -1662.5 * Re_inv - 5416700.0 * Re_inv2 + 0.5191;}
	else{CdRE = 0.4;}
	return CdRE / par.Cc;
}

double
dragForceSingh::computeCd(particle &par){
	double Ts_T, Us_U, Ms, alpha, C1, Br, fKnWr;
	double Cdfm=0;
  	if(par.Mach<1){
		Ts_T=1;
		Us_U=1;
		Ms=par.Mach;
		alpha=1;
		C1=1;
		Br=0;
		fKnWr=1/par.Cc;
	}
	else{
		double MachSQ=par.Mach*par.Mach;
		double s_C1=(gam_m*MachSQ+2);
		double s_C2=(2*gam*MachSQ-gam_m);
		Ts_T=s_C1*s_C2/(gam_pSQ*MachSQ);
		Us_U=s_C1/(gam_p*MachSQ);
		Ms=sqrt(s_C1/s_C2);
		double alpha_inv=alpha0*par.Mach+1-alpha0;
		alpha=1/alpha_inv;
		C1=C1C1/(1-C1C2/par.Mach);
	}

	double Wr=pow(par.Mach,2*omega)/par.Re;
	double WrT=Wr*pow(1+1/Ts_T,omega);
	double Mpow_TOMO=pow(par.Mach,2*omega-1);
	Br=WrT*(Mpow_TOMO+1)/Mpow_TOMO;

	double fKnWr_inv=par.Cc*(1+alpha_hoc*WrT);
	fKnWr=1/fKnWr_inv;

	double s=par.Mach*sqrt(gam*0.5);
	double ss=s*s;
	double ssss=ss*ss;
	Cdfm+=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI));
	Cdfm+=(4*ssss+4*ss-1)*erf(s)/(2*ssss);
	Cdfm+=2/(3*s)*sqrt(M_PI);

	double Ths=pow(1+gam_m*Ms*Ms*0.5,gam/gam_m);
	double Re_=par.Re*pow(1/(alpha*alpha*Ts_T),omega)*pow(Ths,Re_C);
	double C2=1+delta0/sqrt(Re_);
	double Cdc=C1*(1-alpha*Us_U)+C0*Ths*C2*C2;
	double Brita=pow(Br,ita);

	return (Cdc*fKnWr+Cdfm*Brita)/(1+Brita);
}

double
dragForceLoth::computeCd(particle &par){
	double M2=par.Mach*par.Mach;
	double M4=M2*M2;
	double Cd;
	if(par.Re<=45){
		double s=par.Mach*sqrt(gam*0.5);
		double ss=s*s;
		double ssss=ss*ss;
		double Cdfmd=(1+2*ss)*exp(-ss)/(ss*s*sqrt(M_PI))+(4*ssss+4*ss-1)*erf(s)/(2*ssss);
		double Cdfm=Cdfmd+2/(3*s)*sqrt(M_PI);
		double CdfmRe_inv=1+(Cdfmd/1.63-1)*sqrt(par.Re/45.0);
		double CdfmRe=Cdfm/CdfmRe_inv;
		double CdKnRe=24/par.Re*(1+0.15*pow(par.Re,0.687))/par.Cc;
		Cd=(CdKnRe+M4*CdfmRe)/(1+M4);
	}
	if(par.Re>45){
		double CM, GM;
		if(par.Mach>1.45) CM=2.044+0.2*exp(-1.8*log(par.Mach/1.5)*log(par.Mach/1.5));
		if(par.Mach<=1.45) CM=5/3.0+2/3.0*tanh(3*log(par.Mach+0.1));
		if(par.Mach>0.89) GM=0.0002+0.0008*tanh(12.77*(par.Mach-2.02));
		if(par.Mach<=0.89) GM=1-1.525*M4;
		double HM=1-0.258*CM/(1+514*GM);
		Cd=24/par.Re*(1+0.15*pow(par.Re,0.687))*HM;
		Cd+=0.42*CM/(1+42500*GM/pow(par.Re,1.16));
	}
	return Cd;
}
