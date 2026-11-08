/* Runnable driver for the repository's serial n2 kernel.
 * Uses the same initialization, dt, and leapfrog as rushil/barnes_hut_nbody.c.
 * Build: make -C rushil
 * Run: ./rushil/naive [bodies>=2] [steps>=1]
 * 
 * for 4000 bodies and 2000 steps, the time is 25.4 s
 */

#include "physics.h"
#include "n2_internal.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static void acceleration(Simulation *simulation)
{
    SimulationView state=simulation_view(simulation);
    MutableSimulationView out=simulation_mutable_view(simulation);
    n2_accelerations(&state,out.ax,out.ay,out.az);
}

static void step(Simulation *simulation, double dt)
{
    MutableSimulationView s=simulation_mutable_view(simulation);
    for (size_t i=0; i<s.count; ++i) {
        s.vx[i]+=0.5*dt*s.ax[i];
        s.vy[i]+=0.5*dt*s.ay[i];
        s.vz[i]+=0.5*dt*s.az[i];
        s.x[i]+=dt*s.vx[i];
        s.y[i]+=dt*s.vy[i];
        s.z[i]+=dt*s.vz[i];
    }
    acceleration(simulation);
    for (size_t i=0; i<s.count; ++i) {
        s.vx[i]+=0.5*dt*s.ax[i];
        s.vy[i]+=0.5*dt*s.ay[i];
        s.vz[i]+=0.5*dt*s.az[i];
    }
}

int main(int argc, char **argv)
{
    int n=argc>1 ? atoi(argv[1]) : 40;
    int steps=argc>2 ? atoi(argv[2]) : 20000;
    if (argc>3 || n<2 || steps<1) {
        fprintf(stderr,"usage: %s [bodies>=2] [steps>=1]\n",argv[0]);
        return 2;
    }
    Simulation *simulation=NULL;
    if (!simulation_create(&simulation,(size_t)n)) {
        fprintf(stderr,"allocation failed\n"); return 1;
    }
    srand(42); init_bodies(simulation);
    double initial[6], final[6];
    calc_cm(initial,simulation);
    double start=omp_get_wtime();
    acceleration(simulation);
    for (int k=1; k<=steps; ++k) {
        step(simulation,0.002);
        if (k%100==0) {
            calc_cm(final,simulation);
            printf("Center of Mass: %.9g, %.9g, %.9g\n",final[0],final[1],final[2]);
            printf("Velocity of CM: %.9g, %.9g, %.9g\n",final[3],final[4],final[5]);
        }
    }
    fprintf(stderr,"Execution time: %.3g s\n",omp_get_wtime()-start);
    calc_cm(final,simulation);
    double dr=0, dv=0;
    for (int k=0; k<3; ++k) {
        dr+=(final[k]-initial[k])*(final[k]-initial[k]);
        dv+=(final[k+3]-initial[k+3])*(final[k+3]-initial[k+3]);
    }
    printf("Total position drift: %.9g\nTotal velocity drift: %.9g\n",sqrt(dr),sqrt(dv));
    simulation_destroy(simulation);
    return 0;
}
