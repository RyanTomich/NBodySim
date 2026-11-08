#include "physics.h"
#include "n2_internal.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    Simulation *simulation=NULL;
    assert(simulation_create(&simulation,16));
    srand(42); init_bodies(simulation);
    SimulationView view=simulation_view(simulation);
    srand(42);
    for(size_t i=0;i<view.count;i++) {
        double expected[3];
        for(int k=0;k<3;k++) expected[k]=-1.0+2.0*rand()/(double)RAND_MAX;
        assert(view.x[i]==expected[0] && view.y[i]==expected[1] && view.z[i]==expected[2]);
        assert(view.w[i]==1.0/16);
        assert(view.vx[i]==0 && view.vy[i]==0 && view.vz[i]==0);
        assert(view.ax[i]==0 && view.ay[i]==0 && view.az[i]==0);
    }
    simulation_destroy(simulation);
    assert(simulation_create(&simulation,2));
    MutableSimulationView s=simulation_mutable_view(simulation);
    s.w[0]=2; s.w[1]=3;
    s.x[0]=0; s.x[1]=2; s.y[0]=s.y[1]=s.z[0]=s.z[1]=0;
    s.vx[0]=1; s.vx[1]=-1; s.vy[0]=s.vy[1]=s.vz[0]=s.vz[1]=0;
    double cm[6]={99,99,99,99,99,99};
    calc_cm(cm,simulation);
    assert(fabs(cm[0]-1.2)<1e-15 && fabs(cm[3]+0.2)<1e-15);
    assert(cm[1]==0 && cm[2]==0 && cm[4]==0 && cm[5]==0);
    calc_cm(cm,simulation);
    assert(fabs(cm[0]-1.2)<1e-15);
    view=simulation_view(simulation);
    n2_accelerations(&view,s.ax,s.ay,s.az);
    assert(s.ax[0]==0.75 && s.ax[1]==-0.5);
    assert(s.ay[0]==0 && s.ay[1]==0 && s.az[0]==0 && s.az[1]==0);
    simulation_destroy(simulation);
    puts("Simulation initialization, CM, and force tests passed.");
}
