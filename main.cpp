// We are going to break this up into parts
// Part A: 
// We will make a struct for a single particle and its position
// We will then define a "rotation" function, to rotate a vector about origin
// Part B: 
// We will make an efficient 1D indexable array for all particles
// after intialising them.
// We will create another 1D indexable array for the "bins" of size a.
// Part C:
// We will make a main loop that runs for time_steps time steps.
// Inside each time step: 
// First we will do the streaming step. Update all particle velocities.
// Next we will do the binning step, where we update what particles are in
// each bin according to their position.
// Then we will loop over each bin, and do the collision step.
// The collision step is detailed on the wikipedia, essentially consisting
// of updating the velocities based off a rotation matrix from the COM,
// where the rotation is random in 3D.
// The process loop then repeats.
// Part D:
// This will be the visualisation step, where we use C++ and open-gl to 
// visualise what is happening.

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159365258979
#define grid_size 100
#define particle_number 100
#define time_steps 1000

typedef struct {
    double position[2];
    double velocity[2];
} particle;

double distance(particle const& a, particle const& b) {
    double dx = a.position[0] - b.position[0];
    double dy = a.position[1] - b.position[1];
    return sqrt(dx*dx + dy*dy);
};

typedef struct {
    particle *particles=NULL;
    size_t pnumber=0;
    size_t capacity=0;
} bin;

void add_particle(bin *bina, particle particla) {
    // Adds the second argument (particle) to the first argument (bin)
    // If the particle number currently equals the capacity, capacity
    // is extended. 
    if(bina->capacity == bina->pnumber){
        size_t newcap = bina->capacity ? bina->capacity*2 : 20;
        particle *newparticles= (particle*)realloc(bina->particles, newcap*sizeof(particle));
        bina->particles = newparticles;
        bina->capacity = newcap;
    }
    bina->particles[bina->pnumber++]=particla;
};

// We use a fixed array here.
bin bin_grid[grid_size*grid_size];
#define bin_grid(x,y) bin_grid[x + y*grid_size]
// because x goes from 0 to grid_size - 1...

// Use malloc here because this will be big.
bin *bin_grid_storage = (bin *)malloc(grid_size*grid_size*time_steps*sizeof(bin));
// A shortcut for accessing the bin grid storage.
#define bin_grid_storage(l,x,y) bin_grid_storage[(l*grid_size + y)*grid_size +x]

int main() {
    particle billiebob;
    add_particle(bin_grid_storage(1,0,0), billiebob);
    return 0;
;}

