/*
 * Barnes-Hut N-body simulation with unsoftened Newtonian gravity (G = 1).
 *
 * clang -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp barnes_hut_nbody.c physics.c state.c -O3 -o barnes_hut
 * 
 * 14.1 s of execution time for 40 bodies and 20000 steps
 * 12.7 s of execution time for 40 bodies and 20000 steps when using reverse AD
 * 
 * The issue with making this fully general is that you would need to know the accelerations of the other objects
 * in your system in order to compute the current acceleration of your object. So, you would need to solve a system of
 * linear equations, basically making the reverse AD pointless.
 * 
 * It's now 5.21 s for 40 bodies and 20000 steps
 * It's 13.3 s for 100 bodies and 20000 steps
 * 
 * It's 8.18 s for 100 bodies and 20000 steps without using reverse AD
 * It's 21 s for 100 bodies and 20000 steps with reverse AD
 * 
 * 4.33 s for calc time
 * 2.81 s for traversing tree from root up
 * 
 * I think the recursion helps, it's just that the base case needs to be pretty large, like greater than 10000
 * I tried making the check for the relative error, but it's just really slow because calculating the accelerations is a slow process and not
 * needed. I think it's still pretty accurate at least for gravity sims, but maybe for other potentials which don't decay as inverse square, you
 * may notice a difference
 * 
 * 
 * For 4000 bodies and 2000 steps, the time is 10.7s
 */

#include "physics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* Override at build time for task-cutoff benchmarks. */
#ifndef SUMMARIZE_BASE_CASE
#define SUMMARIZE_BASE_CASE 256
#endif
enum { OUTPUT_EVERY = 100 };
typedef struct {
    double lo[3], hi[3], mass, com[3], vcom[3];
    int child[8], body, count;
} Node;

// We index the child by the octant index within the parent

static Simulation *simulation;
static MutableSimulationView bodies;
/* Coordinate pointer tables borrow the Simulation arrays; no data is copied. */
static double *position[3], *velocity[3], *acceleration[3];
static Node *tree;
static int n, nodes, capacity;
static double theta = 0.6;


static int new_node(const double lo[3], const double hi[3]) {
    if (nodes == capacity) { fprintf(stderr, "tree capacity exhausted\n"); exit(1); }
    Node *q = &tree[nodes];
    *q = (Node){0};
    for (int k = 0; k < 3; ++k) { q->lo[k] = lo[k]; q->hi[k] = hi[k]; }
    for (int k = 0; k < 8; ++k) q->child[k] = -1;
    q->body = -1;
    return nodes++;
}

static int octant(const Node *q, int bi) {
    double mid[3] = {(q->lo[0]+q->hi[0])/2, (q->lo[1]+q->hi[1])/2, (q->lo[2]+q->hi[2])/2};
    return (position[0][bi] >= mid[0]) | ((position[1][bi] >= mid[1]) << 1) | ((position[2][bi] >= mid[2]) << 2);
}

static int make_child(int parent, int o) {
    double lo[3], hi[3]; Node *p = &tree[parent];
    for (int k = 0; k < 3; ++k) {
        double mid = (p->lo[k] + p->hi[k]) / 2;
        lo[k] = (o & (1<<k)) ? mid : p->lo[k];
        hi[k] = (o & (1<<k)) ? p->hi[k] : mid;
    }

    int c = new_node(lo, hi);
    tree[parent].child[o] = c;
    return c;
}

// We insert body bi into node ni
static void insert(int ni, int bi, int depth) {
    Node *q = &tree[ni]; q->count++;
    if (q->body < 0) {
        int has_child = 0; for (int k=0;k<8;k++) has_child |= q->child[k] >= 0;
        if (!has_child) { q->body = bi; return; }
    }
    if (q->body >= 0) {
        int old = q->body; q->body = -1;
        int oo = octant(q, old); if (q->child[oo] < 0) make_child(ni, oo);
        insert(q->child[oo], old, depth+1);
    }
    int o = octant(q, bi);

    if (q->child[o] < 0){
        make_child(ni, o); // Because child doesn't exist if nothing is in there
    }

    insert(q->child[o], bi, depth+1);
}

static void summarize(int ni) {
    Node *q = &tree[ni];
    q->mass = 0;
    for (int k=0;k<3;k++) q->com[k]=q->vcom[k]=0;

    // If it stores an element, then it doesn't have any children by definition
    if (q->body >= 0) {
        int bi=q->body; q->mass=bodies.w[bi];
        for(int k=0;k<3;k++) { q->com[k]=position[k][bi]; q->vcom[k]=velocity[k][bi]; }
        return;
    }

    for (int c=0;c<8;c++) if(q->child[c]>=0) {
        int child=q->child[c];
        if (q->count > SUMMARIZE_BASE_CASE) {
            #pragma omp task firstprivate(child)
            summarize(child);
        } else {
            summarize(child);
        }
    }

    #pragma omp taskwait

    for (int c=0;c<8;c++) if(q->child[c]>=0) {
        Node *s=&tree[q->child[c]];
        q->mass += s->mass;
        for(int k=0;k<3;k++) { q->com[k]+=s->mass*s->com[k]; q->vcom[k]+=s->mass*s->vcom[k]; }
    }
    for(int k=0;k<3;k++) { q->com[k]/=q->mass; q->vcom[k]/=q->mass; }
}

static void source_acceleration(int i, double mass, const double x[3],
                                const double v[3], double out_a[3]) {
    (void)v; /* Newtonian gravity depends only on positions and source mass. */
    double d[3], r2=0;
    for (int k=0; k<3; ++k) {
        d[k]=x[k]-position[k][i];
        r2+=d[k]*d[k];
    }
    /* Match n2.c, including its convention for coincident positions. */
    double scale=r2 == 0 ? 0 : mass/(r2*sqrt(r2));
    for (int k=0; k<3; ++k) out_a[k]=scale*d[k];
}

static void add_acceleration(int i, const double a[3]) {
    for (int k=0; k<3; ++k) acceleration[k][i] += a[k];
}

static void add_source(int i, double mass, const double x[3], const double v[3]) {
    double a[3];
    source_acceleration(i, mass, x, v, a);
    add_acceleration(i, a);
}

static void walk(int i, int ni) {
    Node *q=&tree[ni]; if(!q->mass) return;
    if(q->body == i) return;

    double d2=0, width=q->hi[0]-q->lo[0];
    for(int k=0;k<3;k++){ double d=q->com[k]-position[k][i]; d2+=d*d; }

    // If the node is a leaf, then we just add the acceleration from the leaf to the object
    if (q->body >= 0) {
        add_source(i,q->mass,q->com,q->vcom); return;
    }

    int contains = position[0][i]>=q->lo[0] && position[0][i]<=q->hi[0] && position[1][i]>=q->lo[1] && position[1][i]<=q->hi[1] && position[2][i]>=q->lo[2] && position[2][i]<=q->hi[2];

    if (!contains && width*width < theta*theta*d2) {
        add_source(i,q->mass,q->com,q->vcom);
        return;
    }

    for(int c=0;c<8;c++) if(q->child[c]>=0) walk(i,q->child[c]);
}

/* The topology remains valid while every body stays inside its leaf box.
 * Checking leaves also ensures bodies remain inside all ancestor boxes. */
static int tree_is_valid(void) {
    if (!nodes) return 0;
    for (int ni=0; ni<nodes; ++ni) {
        const Node *q=&tree[ni];
        if (q->body < 0) continue;
        for (int k=0; k<3; ++k)
            if (!(position[k][q->body]>=q->lo[k] &&
                  position[k][q->body]<=q->hi[k])) return 0;
    }
    return 1;
}

static void rebuild_tree(void) {
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};

    for(int i=0;i<n;i++) for(int k=0;k<3;k++){
        if(position[k][i]<lo[k])lo[k]=position[k][i];
        if(position[k][i]>hi[k])hi[k]=position[k][i];
    }

    double side=hi[0]-lo[0];
    for(int k=1;k<3;k++) if(hi[k]-lo[k]>side)side=hi[k]-lo[k];

    side = side ? side*1.001 : 1; // added a boundary around the octant for fp ops
    
    for(int k=0;k<3;k++){
        double mid=(lo[k]+hi[k])/2;
        lo[k]=mid-side/2; hi[k]=mid+side/2;
    }

    nodes=0;
    int root=new_node(lo,hi);
    for(int i=0;i<n;i++) insert(root,i,0);
}

static double accelerations(void) {
    if (!tree_is_valid()) rebuild_tree();

    #pragma omp parallel
    {
        #pragma omp single
        summarize(0);
    }
    double start_time = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic) // for uneven calculation of updates
    for(int i=0;i<n;i++) {
        for(int k=0;k<3;k++) acceleration[k][i]=0;
        walk(i,0);
    }

    double end_time = omp_get_wtime();

    return end_time-start_time;
}

static double step(double dt) {
    /* Newtonian forces depend only on position, so reuse the last acceleration. */
    double calc_time = 0;
    if (!nodes) calc_time += accelerations();
    #pragma omp parallel for
    for(int i=0;i<n;i++) {
        for(int k=0;k<3;k++){
            velocity[k][i]+=0.5*dt*acceleration[k][i];
            position[k][i]+=dt*velocity[k][i];
        }
    }
    calc_time += accelerations();
    #pragma omp parallel for
    for(int i=0;i<n;i++){
        for(int k=0;k<3;k++){
            velocity[k][i]+=0.5*dt*acceleration[k][i];
        }
    }
    return calc_time;
}

int main(int argc,char **argv) {
    double start_time = omp_get_wtime();

    n = argc>1 ? atoi(argv[1]) : 500; int steps=argc>2 ? atoi(argv[2]) : 1000;
    if(argc>3) theta=atof(argv[3]);
    if(argc>4 || n<2 || steps<1 || theta<=0) {
        fprintf(stderr,"usage: %s [bodies>=2] [steps>=1] [theta>0]\n",argv[0]); return 2;
    }
    if (!simulation_create(&simulation,(size_t)n)) {
        fprintf(stderr,"allocation failed\n"); return 1;
    }
    bodies=simulation_mutable_view(simulation);
    position[0]=bodies.x; position[1]=bodies.y; position[2]=bodies.z;
    velocity[0]=bodies.vx; velocity[1]=bodies.vy; velocity[2]=bodies.vz;
    acceleration[0]=bodies.ax; acceleration[1]=bodies.ay; acceleration[2]=bodies.az;
    capacity=16*n+1; // we arbitrarily choose this bc we dont know the size of the tree beforehand
    tree=calloc((size_t)capacity,sizeof(*tree));
    if(!tree) { simulation_destroy(simulation); return 1; }

    srand(42); init_bodies(simulation);
    puts("step,id,x,y,z,vx,vy,vz");
    
    double execution_time = omp_get_wtime();

    double initial_loc[6] = {0, 0, 0, 0, 0, 0};
    double final_loc[6] = {0, 0, 0, 0, 0, 0};
    double tot_calc_time = 0;

    for (int s=1; s<=steps; ++s) {
        if (s == 1){
            calc_cm(initial_loc, simulation);
        }
        double cm_vals[6] = {0, 0, 0, 0, 0, 0};

        tot_calc_time += step(0.002);

        calc_cm(cm_vals, simulation);
        if (s == steps){
            calc_cm(final_loc, simulation);
        }
        if (s % OUTPUT_EVERY == 0) {
            printf("Center of Mass: %.9g, %.9g, %.9g\n", cm_vals[0], cm_vals[1], cm_vals[2]);
            printf("Velocity of CM: %.9g, %.9g, %.9g\n", cm_vals[3], cm_vals[4], cm_vals[5]);
        }
    }

    double total_drift = 0;
    double total_velocity_drift = 0;
    for (int k = 0; k < 3; k++){
        total_drift += (final_loc[k] - initial_loc[k]) * (final_loc[k] - initial_loc[k]);
        total_velocity_drift += (final_loc[k + 3] - initial_loc[k + 3]) * (final_loc[k + 3] - initial_loc[k + 3]);
    }

    printf("Total position drift: %.9g\n", sqrt(total_drift));
    printf("Total velocity drift: %.9g\n", sqrt(total_velocity_drift));

    double end_time = omp_get_wtime();

    fprintf(stderr,"Initialization time: %.3g s\n",execution_time-start_time);
    fprintf(stderr,"Execution time: %.3g s\n",end_time-execution_time);
    fprintf(stderr,"Total calculation time: %.3g s\n",tot_calc_time);
    free(tree); simulation_destroy(simulation); return 0;
}
