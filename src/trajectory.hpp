#pragma once
#include "forces/dragForce.hpp"

class trajectory{
	public:
		Variables *vars;	// pointer for variables
		Flags *flags;			// pointer for flags
		std::vector<dragForceSM*> forces;	// pointer array for force functions working on the particle

		int Nth;	// Number of parallel (OpenMP)
		double rho_p;									// particle density
		double observeTime;						// output interval in time scale
		double totalTime;							// total calculation time
		int Observe;									// output interval in time step
		double delta_r2_trajectory;		// squre of migration distance for trajectory output
		int plane2D;									// face of 2D simulation (0:y-z, 1:x-z, 2:x-y)
		string startDir;							// start directory name (20000 is default)
		double Axis;									// axis of 2D axi-symmetric simulation
		double flag;									// NOT USED?
		char filepath[100];						// file path
		int boundaryStartID;					// boundary start face id
		double constTemp;							// constant temperature for incompressible flow
		double constRho;							// constant density for incompressible flow
		double timer;
		FILE*f;
		FILE*file;
		std::vector<outParticle> outParticles; // position & velocity of the particles at the end of the simulation
		std::vector<penetrate> penetrates;
		std::vector<int> trapParticle;				 // number of particles trapped in the curculation

		// The functions for turbulent dispersion (see turbulentDispersion.cpp)
		double get_tini(particle &par);
		void updateDisp(particle &par);

		// Main simulation functions
		void run(void);
		double timeEvolution(particle &a);
		double euler(particle &a);
		double analytical(particle &a);
		int checkCell(int pid);
		int boundAction(int faceID, int pid, point norm);

		// Physical properties calculation functions (see calcProperties.cpp)
		void calculateMyu(void);
		void calculatelamda(void);
		void calculateNonDimension(particle &par);


		// Simualtion geometrical information
		std::vector<face> faces;
		std::vector<point> points;
		std::vector<boundary> boundaries;
		std::vector<cell> cells;
		std::vector<int> neighbors;
		std::vector<int> owners;
		void makeCells(void); // make cells from abouve geometrical arrays

		// Reading functions
		void readGeometry(void);
		void readFaces(void);
		void readBoundaries(void);
		void readPoints(void);
		void readNeighbors(void);
		void readOwners(void);
		void readScalar(char *readFile, std::vector<double> &variable);
		void readScalarDum(char *readFile, std::vector<double> &variable, double value);
		void readVector(char *readFile, std::vector<point> &variable);
		void readVectorDum(char *readFile, std::vector<point> &variable, double value);
		void readParticles(void);
		void readCondition(void);
		void readCFDresults(void);

		// Initialization
		void initialParticle(void);
		void findParticle(void);
		void findParticleFace(std::vector<cell> cells);

		// Export functions
		void output(particle a, double time, int nth);
		void outputTrajectory(particle a);
		void outputInitial(void);
		void outputFinalPosition(void);
		void outputPenetrate(void);
		void outputInitialVelocity(void);

		void initialize(particle &par, int nth){
		if(flags->dispersionFlag) par.update=1;
			else {
				par.tini=1e10;
				par.update=0;
			}
	}

	trajectory(void);
	~trajectory(void);

	private:
};
