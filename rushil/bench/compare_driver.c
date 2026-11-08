/* Same leapfrog workload for both implementations. The rushil adapter includes
 * its standalone source to exercise the actual static tree/force functions. */
#ifdef RUSHIL_BACKEND
#define main rushil_original_main
#include "../barnes_hut_nbody.c"
#undef main
#elif defined(NAIVE_BACKEND)
#define main naive_original_main
#include "../naive.c"
#undef main
#else
#include "solver.h"
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#endif
#include <limits.h>
#include <string.h>

static void evaluate(Simulation *state
#if !defined(RUSHIL_BACKEND) && !defined(NAIVE_BACKEND)
                     , Solver *solver
#endif
)
{
#ifdef RUSHIL_BACKEND
    (void)state;
    (void)accelerations();
#elif defined(NAIVE_BACKEND)
    acceleration(state);
#else
    SimulationView view=simulation_view(state);
    MutableSimulationView out=simulation_mutable_view(state);
    solver->acceleration(&view,out.ax,out.ay,out.az,solver->context);
#endif
}

/* Identical integration arithmetic and scheduling for both force backends. */
static void integrate(Simulation *state, double dt, int steps
#if !defined(RUSHIL_BACKEND) && !defined(NAIVE_BACKEND)
                      , Solver *solver
#endif
)
{
    MutableSimulationView v=simulation_mutable_view(state);
    double *r[]={v.x,v.y,v.z}, *vel[]={v.vx,v.vy,v.vz}, *a[]={v.ax,v.ay,v.az};
    evaluate(state
#if !defined(RUSHIL_BACKEND) && !defined(NAIVE_BACKEND)
             ,solver
#endif
    );
    for(int s=0;s<steps;s++) {
        #pragma omp parallel for
        for(size_t i=0;i<v.count;i++) for(int k=0;k<3;k++) {
            vel[k][i]+=0.5*dt*a[k][i];
            r[k][i]+=dt*vel[k][i];
        }
        evaluate(state
#if !defined(RUSHIL_BACKEND) && !defined(NAIVE_BACKEND)
                 ,solver
#endif
        );
        #pragma omp parallel for
        for(size_t i=0;i<v.count;i++) for(int k=0;k<3;k++)
            vel[k][i]+=0.5*dt*a[k][i];
    }
}

int main(int argc,char **argv)
{
    if(argc!=9) {
        fprintf(stderr,"usage: driver INPUT DT STEPS SAMPLES WARMUPS THETA FINAL BACKEND\n");
        return 2;
    }
    double dt=strtod(argv[2],NULL), opening=strtod(argv[6],NULL);
    int steps=atoi(argv[3]), samples=atoi(argv[4]), warmups=atoi(argv[5]);
    if(!isfinite(dt)||dt<=0||!isfinite(opening)||opening<=0||steps<1||samples<2||warmups<0)
        return 2;
    FILE *input=fopen(argv[1],"r");
    size_t count;
    if(!input || fscanf(input,"%zu",&count)!=1 || count<2 || count>(INT_MAX-1)/16) return 2;
    Simulation *state=NULL;
    if(!simulation_create(&state,count)) return 1;
    MutableSimulationView v=simulation_mutable_view(state);
    double *fields[]={v.w,v.x,v.y,v.z,v.vx,v.vy,v.vz,v.ax,v.ay,v.az};
    double *saved=malloc(10*count*sizeof *saved);
    if(!saved) return 1;
    for(size_t i=0;i<count;i++) for(int k=0;k<10;k++) {
        if(fscanf(input,"%lf",&fields[k][i])!=1 || !isfinite(fields[k][i])) return 2;
    }
    fclose(input);
    for(int k=0;k<10;k++) memcpy(saved+k*count,fields[k],count*sizeof(double));
#ifdef RUSHIL_BACKEND
    (void)argv[8];
    simulation=state; bodies=v; n=(int)count; theta=opening;
    position[0]=v.x; position[1]=v.y; position[2]=v.z;
    velocity[0]=v.vx; velocity[1]=v.vy; velocity[2]=v.vz;
    acceleration[0]=v.ax; acceleration[1]=v.ay; acceleration[2]=v.az;
    capacity=16*n+1;
    tree=calloc((size_t)capacity,sizeof *tree);
    if(!tree) return 1;
#elif defined(NAIVE_BACKEND)
    (void)argv[8];
#else
    Solver solver;
    const char *error;
    SolverOption option={"theta",argv[6]};
    int is_bh=strcmp(argv[8],"barnes-hut")==0;
    if(!solver_create(&solver,argv[8],count,is_bh?&option:NULL,is_bh?1:0,&error)) {
        fprintf(stderr,"%s\n",error); return 1;
    }
#endif
    puts("sample,seconds");
    for(int sample=-warmups;sample<samples;sample++) {
        for(int k=0;k<10;k++) memcpy(fields[k],saved+k*count,count*sizeof(double));
#ifdef RUSHIL_BACKEND
        nodes=0; /* Every run starts with a fresh topology, with allocation reused. */
#endif
        double start=omp_get_wtime();
        integrate(state,dt,steps
#if !defined(RUSHIL_BACKEND) && !defined(NAIVE_BACKEND)
                  ,&solver
#endif
        );
        double duration=omp_get_wtime()-start;
        for(int k=0;k<10;k++) for(size_t i=0;i<count;i++)
            if(!isfinite(fields[k][i])) { fprintf(stderr,"nonfinite final state\n"); return 1; }
        if(sample>=0) printf("%d,%.17g\n",sample+1,duration);
    }
    FILE *output=fopen(argv[7],"w");
    if(!output) return 1;
    fprintf(output,"%zu\n",count);
    for(size_t i=0;i<count;i++) {
        for(int k=0;k<10;k++) fprintf(output,"%s%.17g",k?" ":"",fields[k][i]);
        fputc('\n',output);
    }
    if(fclose(output)) return 1;
#ifdef RUSHIL_BACKEND
    free(tree);
#elif defined(NAIVE_BACKEND)
    /* No prepared solver workspace is needed by n2. */
#else
    solver_destroy(&solver);
#endif
    free(saved); simulation_destroy(state);
    return 0;
}
