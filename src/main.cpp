#include "trajectory.hpp"



int main(int argc,char *argv[]){

	trajectory *tra = new trajectory( );		// trajectory class is initialized with constractor (see trajectory.cpp)
	/*particle par;
	for(par.Mach=0.1;par.Mach<5;par.Mach+=0.02){
		par.Re=200;
		par.Kn=par.Mach/par.Re*Kn_coeff;
		par.Cc=1+par.Kn*(A1+A2*exp(-A3/par.Kn));
		double Cd=tra->forces[0]->computeCd(par);
		cout<<par.Re<<" "<<par.Mach<<" "<<Cd<<endl;
	}
	for(float iRe=-4; iRe<4; iRe+=0.1){
		par.Re=pow(10,iRe);
			par.Mach=1.5;
			par.Kn=par.Mach/par.Re*Kn_coeff;
			par.Cc=1+par.Kn*(A1+A2*exp(-A3/par.Kn));
			double Cd=tra->forces[0]->computeCd(par);
			cout<<par.Re<<" "<<par.Mach<<" "<<Cd<<endl;
	}*/
	if (tra->flags->initialBreakFlag==0) tra->run();  // This is main part of this simulaiton
	return 0;
}
