#include "physics.h"
#include <stdlib.h>
#include <string.h>

static double uniform(double lo, double hi)
{
    return lo + (hi-lo) * rand() / (double)RAND_MAX;
}

void init_bodies(Simulation *simulation)
{
    MutableSimulationView s=simulation_mutable_view(simulation);
    for (size_t i=0; i<s.count; ++i) {
        s.w[i]=1.0/s.count;
        s.x[i]=uniform(-1.0,1.0);
        s.y[i]=uniform(-1.0,1.0);
        s.z[i]=uniform(-1.0,1.0);
        s.vx[i]=s.vy[i]=s.vz[i]=0;
        s.ax[i]=s.ay[i]=s.az[i]=0;
    }
}

void calc_cm(double cm_vals[6], const Simulation *simulation)
{
    SimulationView s=simulation_view(simulation);
    memset(cm_vals,0,6*sizeof *cm_vals);
    double mass=0;
    for (size_t i=0; i<s.count; ++i) {
        mass+=s.w[i];
        cm_vals[0]+=s.w[i]*s.x[i];
        cm_vals[1]+=s.w[i]*s.y[i];
        cm_vals[2]+=s.w[i]*s.z[i];
        cm_vals[3]+=s.w[i]*s.vx[i];
        cm_vals[4]+=s.w[i]*s.vy[i];
        cm_vals[5]+=s.w[i]*s.vz[i];
    }
    for (int k=0; k<6; ++k) cm_vals[k]/=mass;
}
