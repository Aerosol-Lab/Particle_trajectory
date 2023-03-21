#include "trajectory.hpp"


trajectory::trajectory(void){
// Get number of thread for OpenMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      Nth=omp_get_num_threads();
    }
  }

// set default values
  rho_p=1000; // particle density
  observeTime=1e-5; // output time interval
  totalTime=1;    // total calculation time
  Observe=10000;	// output time steps
  delta_r2_trajectory=1e-4*1e-4;	// output migration distance for the trajectory
  plane2D=-1; //
  startDir="20000"; // CFD simulation directory
  Axis=0; // 0 for 3D simulaiton
  flag=1;
  constTemp=300;  // constant temperature for the incompressible flow
  constRho=1.2;   // constant density for the incompressible flow
  for(int i=0;i<Nth;i++) trapParticle.push_back(0);

// generate class
  vars = new Variables(); // generate variables class
  flags = new Flags();    // initialize flags (see flags.hpp)
  

// initialization
  readCondition();  // read condition file ./particle/condition
  readGeometry();   // read geometry ./constant/polyMesh/
  readCFDresults(); // CFD result ./${startDir}/

  calculateMyu();   // viscosity calculation from field data
  calculatelamda(); // mean free path calculation form field data

  readParticles();  // read particle file ./particle/particleSet
  makeCells();      // make cells data from geometry file
  initialParticle();// initialize particle (see below)
  outputInitial();  // initialization of output files
}

void
trajectory::initialParticle(void){

	findParticle();  // find initial location of particles

  int ps=vars->particles.size();
	for(auto &a:vars->particles){
		int icell=a.cell;
		a.dt=1/(40*a.beta);
		for(int i=0; i<3; i++) a.v.x[i]=vars->U[icell].x[i]*0.99;   // initial velocity is 99% of fluid velocity
  	calculateNonDimension(a);
		for(auto &force : forces) force->computeFD(vars,flags,a);   // force calculation

    double threePiMuDp_Cc=3*M_PI*vars->myu[icell]*a.dp/a.Cc;
    a.Zp=1.6e-19/threePiMuDp_Cc;
		a.beta=threePiMuDp_Cc/a.m;

  // initialize outParticles array (ending information of particles)
    outParticle op;
    op.pid=a.id;
    op.bid=0;
    point dum;
    dum.x[0]=dum.x[1]=dum.x[2]=0;
    op.r=dum;
    op.v=dum;
    outParticles.push_back(op);

    for(auto &pen:penetrates){
      pen.outPositions.push_back(op);
      pen.dx0.push_back(pen.loc-a.x.x[pen.face]);
    }
	}
}
