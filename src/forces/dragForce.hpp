#pragma once
#include "../functions.hpp"
#include "../variables.hpp"
#include "../flags.hpp"



//Stokes-Millikan
class dragForceSM{
	public:
		void computeFD(Variables *vars, Flags *flags, particle &par);
		virtual double computeCd(particle &par);
	private:
};

// Morsi & Alexander (1972)
class dragForceMA : public dragForceSM{
	public:
		double computeCd(particle &par);
	private:
};

// Singh et al. (2021)
class dragForceSingh : public dragForceSM{
	public:
		double computeCd(particle &par);
	private:
		//	Parameters for drag force
		const double delta0=9.4; // or 9.06?
		const double alpha0=0.356;
		const double ita=1.8;
		const double alpha_hoc=1.27;
		const double C0=24/delta0/delta0; // C0*delta0*delta0=24
		const double CdMach=0.9;
		const double omega=0.74;
		const double gam_m=gam-1;
		const double gam_mSQ=gam_m*gam_m;
		const double gam_p=gam+1;
		const double gam_pSQ=gam_p*gam_p;

		const double C1C1=CdMach-C0*pow(1+gam_m*gam_m*0.25/gam,gam/gam_m);
		const double C1C2=gam_m/alpha0/gam_p;
		const double Re_C=gam_p*0.5/gam-gam_m/gam*omega;
};

// Loth model
class dragForceLoth : public dragForceSM{
	public:
		double computeCd(particle &par);
	private:
};
