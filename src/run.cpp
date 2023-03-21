#include "trajectory.hpp"

void
trajectory::run(void){
	timer=omp_get_wtime();
	const int particleSize=vars->particles.size();

	# pragma omp parallel
	{
		# pragma omp for
		for (int pid=0; pid<particleSize; pid++){
			int nth=omp_get_thread_num();
			particle &a=vars->particles[pid];
			initialize(a,nth);

			int breakFlag=0;
			double preOutTime=0;		// previous output time
			double time=0;					// current time


			while(time<totalTime){
				if(time - preOutTime > observeTime) {
					output(a,time,nth);
					preOutTime=time;
				}

				for(auto &pen:penetrates){
					if(pen.dx0[pid]!=0){
						double dx=pen.loc-a.x.x[pen.face];
						if(pen.dx0[pid]*dx<0){
							pen.outPositions[pid].pid=a.id;
							pen.outPositions[pid].r=a.x;
							pen.outPositions[pid].v=a.v;
							pen.outPositions[pid].bid=1;
							pen.dx0[pid]=0;
						}
					}
				}

				time+=timeEvolution(a);		// see timeEvolution.cpp

				breakFlag=checkCell(pid);

				if(breakFlag==-1) break;

				if(a.update==1&&flags->dispersionFlag==1) updateDisp(a);
			}

			// if particle did not hit on any "bounday" during calculation time
			if(breakFlag==0){
				outParticles[pid].pid=a.id;
				outParticles[pid].r=a.x;
				outParticles[pid].v=a.v;
				outParticles[pid].bid=-1;
				trapParticle[nth]++;
			}
		}
	}

	outputFinalPosition();
	outputPenetrate();
	int totalTrap=0;
	for(int i=0;i<Nth;i++) totalTrap+=trapParticle[i];
  	cout<<totalTrap<<" might be trapped circulation"<<endl;
	cout<<"Trajectory calculation time: "<<omp_get_wtime()-timer<<" sec"<<endl;
}
