#include "trajectory.hpp"



int main(int argc,char *argv[]){
	trajectory *tra = new trajectory( );		// trajectory class is initialized with constractor (see trajectory.cpp)
	if (tra->flags->initialBreakFlag==0) tra->run();  // This is main part of this simulaiton
	return 0;
}
