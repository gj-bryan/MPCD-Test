// We are going to break this up into parts
// Part A: 
// We will make a struct for a single particle and its position, velocity.
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

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define PI 3.14159365258979
#define grid_size 30
#define particle_number 20
#define time_steps 1000
#define time_dt 0.1

typedef struct {
    double x;
    double y;
} vec2;

typedef struct {
    vec2 position;
    vec2 velocity;
} particle;

// A basic distance function. Playing with C++ syntax here -- & is similar to * as in it is a pointer, but one can use the "." instead of the "->" and it is read only. Fun stuff I suppose, advantage is that it won't break because I modified it stupidly, and I don't have to pay to copy it. Note that to make this backwards compatible with C I'll need to rewrite this as * instead, and remove "const", and make "." into "->".
double distance(particle const &a, particle const &b) {
    double dx = a.position.x - b.position.x;
    double dy = a.position.y - b.position.y;
    return sqrt(dx*dx + dy*dy);
};

// Takes a given vector (vec2 object) and returns the same vector rotated randomly.
vec2 rotate_random(vec2 orig_vector) {
    vec2 new_vector;
    double random_angle = 2*PI*(double)rand()/((double)RAND_MAX);
    double new_x = cos(random_angle)*orig_vector.x - sin(random_angle)*orig_vector.y;
    double new_y = cos(random_angle)*orig_vector.y + sin(random_angle)*orig_vector.x;
    new_vector.x = new_x;
    new_vector.y = new_y;
    return new_vector;
}

#define TEMPERATURE 0.5

// Draws a random number from a Boltzmann distribution
double boltzmann_number() {
    // a random probability
    double r_p = (double)rand()/((double)RAND_MAX);
    double norm_factor = 1/(1-exp(-1/TEMPERATURE));
    // Draws from the boltzmann distribution. This is done by reverse engineering the Boltzmann distribution to take a probability and return the boltzmann number (energy) corresponding to that probability.
    double b_number= TEMPERATURE*log((1/r_p)*norm_factor);
    return b_number;
}

// Draws the velocity corresponding to a Boltzmann distribution
double boltzmann_velocity() {
    double b_number = boltzmann_number();
    double mass = 1;
    return sqrt(2*b_number/mass);
}

// The so-called "Bin" struct. Has particle ids, particle number, and particle capacity as an attribute. Particle capacity is automatically extended if more particle ids are added here than expected. The structure of "*particle_ids" is a pointer to a malloc.
typedef struct {
    size_t *particle_ids=NULL;
    size_t pnumber=0;
    size_t capacity=0;
    vec2 position;
    vec2 size;
} bin;

// We use a fixed array containing all the bins at "present"
bin bin_grid[grid_size*grid_size];
// A shortcut for accessing the bin grid fixed array. N.B. we must use the "&" prefix when accessing this as pointer.
#define bin_grid(x,y) bin_grid[x + y*grid_size]
// because x goes from 0 to grid_size - 1...

// We use a fixed array containing all the particles at "present." In this way, each particle will have a unique ID. We can ccess the particle itself from the bins (which include the indices) or directly here.
particle particle_array[grid_size*grid_size*particle_number];

// A shortcut for accessing the particle array.
#define particle_array(i) particle_array[i]

// Adds a particle ID to a bin. Note that this takes a pointer so be sure to use the "&" when passing something in here.
void add_particle(bin *bina, size_t particle_ida) {
    if(bina->capacity == bina->pnumber){
        size_t newcap = bina->capacity ? bina->capacity*2 : 20;
        size_t *new_particles= (size_t*)realloc(bina->particle_ids, newcap*sizeof(size_t));
        bina->particle_ids = new_particles;
        bina->capacity = newcap;
    }
    bina->particle_ids[bina->pnumber++]=particle_ida;
};

// Removes a particle ID from a bin.
void remove_particle(bin *bina, size_t particle_ida) {
    for(int j=0;j<bina->pnumber;j++) {
        if (particle_ida == bina->particle_ids[j]) {
            // We use the "swap and pop" trick. Essentially we take the last entry and put it here, and decrease the length of the entries (pnumber).
            bina->particle_ids[j]=bina->particle_ids[--bina->pnumber];
            return;
        }
    }
}

// Determines if a particle is physically located inside a proposed bin. Returns true/false if particle is in/outside the domain of the bin.
bool is_in_bin(bin *bina, particle particla) {
    if(particla.position.x >= bina->position.x && particla.position.x < bina->position.x+bina->size.x && particla.position.y >= bina->position.y && particla.position.y < bina->position.y+bina->size.y) {
        return true;
    }
    else {
        return false;
    }
}

int main() {

    setbuf(stdout, NULL);

	srand(time(NULL));

    FILE *log = fopen("debug.log", "w");

    // A malloc of bins representing the history of all bins. 
    bin *bin_grid_storage = (bin*)malloc(grid_size*grid_size*time_steps*sizeof(bin));
    // A shortcut for accessing the "bin grid storage".
    #define BGS(l,x,y) bin_grid_storage[(l*grid_size + y)*grid_size +x]

    particle *particle_storage = (particle*)malloc(grid_size*grid_size*particle_number*time_steps*sizeof(particle));
    #define PS(l,n) particle_storage[(l*grid_size*grid_size*particle_number) + n]

    // INITIALISATION OF THE SIMULATION
    //
    // We will first initialise all the bins so that they are sorted properly according to position. N.B: The bins are initialised so that the top left corner is their position, and the bin sizes are normalised so that the sides of the total box have length 1.
    for (int x = 0; x < grid_size; x++) {
        for (int y= 0; y< grid_size; y++) {
            bin *b = &bin_grid(x,y);
            b->position.x = (double)x/(double)grid_size;
            b->position.y = (double)y/(double)grid_size;
            b->size.x= 1/(double)grid_size;
            b->size.y=1/(double)grid_size;
        }
    }

    // We will now initialise the particles per bin. We should do this by using a Boltzmann distribution of velocities.

    for (int x = 0; x < grid_size; x++) {
        for (int y =0;y<grid_size; y++) {
            for (int j=0; j< particle_number; j++) {
                size_t ID = (size_t)((grid_size*x+y)*particle_number + j);
                particle *p = &particle_array(ID);
                bin *b = &bin_grid(x,y);
                // Initialise a random position in the bin.
                double random_number_x = (double)rand()/((double)RAND_MAX);
                double random_number_y = (double)rand()/((double)RAND_MAX);
                p->position.x = b->position.x + random_number_x/(double)grid_size;
                p->position.y= b->position.y + random_number_y/(double)grid_size;
                // Create random velocity vector from the Boltzmann distribution.
                vec2 new_velocity;
                new_velocity.x = boltzmann_velocity();
                new_velocity.y = 0.0;
                vec2 rot_new_velocity = rotate_random(new_velocity);
                p->velocity.x = rot_new_velocity.x;
                p->velocity.y = rot_new_velocity.y;
                // Add the particle ID to the bin.
                add_particle(b,ID);
            }
        }
    }

    // MAIN LOOP
    //
    //

    for (int t = 0; t < time_steps; t++) {
        // STREAMING STEP
        // 
        // UPDATE PARTICLE VELOCITIES:

        for (int j=0; j< grid_size*grid_size*particle_number; j++) {
            particle *p = &particle_array(j);
            // Stream the velocities so the particles are possibly outside the bin.
            p->position.x += p->velocity.x*time_dt;
            p->position.y += p->velocity.y*time_dt;
            p->position.x = fmod(fmod(p->position.x, 1.0) + 1.0, 1.0);
            p->position.y = fmod(fmod(p->position.y, 1.0) + 1.0, 1.0);  
        }

        // RE-BINNING STEP

        for (int x = 0; x < grid_size; x++) {
            for (int y =0;y<grid_size; y++) {
                bin *b = &bin_grid(x,y);
                for (int j=0; j< b->pnumber; j++) {
                    // The particle is assumed to in the old bin *b
                    size_t ID = b->particle_ids[j];
                    particle *p = &particle_array(ID);
                    // If the particle is not in the bin:
                    if(!is_in_bin(b, *p)) {
                        remove_particle(b, ID);
                        int new_x = (int)(p->position.x*grid_size);
                        int new_y = (int)(p->position.y*grid_size);
                        bin *new_b = &bin_grid(new_x,new_y);
                        add_particle(new_b,ID);
                        j--;
                    }
                }
            }
        } 

        // RECORD PARTICLES

        for(int j=0; j<particle_number*grid_size*grid_size;j++) {
            PS(t,j) = particle_array(j);
        }

        // RECORD BINS

        for(int x=0; x<grid_size; x++) {
            for(int y=0; y<grid_size; y++) {
                BGS(t,x,y) = bin_grid(x,y);
            }
        }
    }
    
    #define WindowSize 200

	// Create folder
	mkdir("frames", 0777);

	// Allocate image buffer
	unsigned char *image = (unsigned char*)malloc(WindowSize * WindowSize * 3);
	int jindex=0;
	// Loop over saved timesteps
	for (int t = 0; t < time_steps; t++) {
		if((t%40)!= 1){continue;}
		// Fill image from stored order parameter
		jindex=jindex+1;
		for (int y = 0; y < WindowSize; y++) {
			for (int x = 0; x < WindowSize; x++) {

				int gx = (int)(x * grid_size / WindowSize);
				int gy = (int)(y * grid_size / WindowSize);

				int density = BGS(t,gx,gy).pnumber;
                double brightness = (double)density/40;

				int index = (y * WindowSize + x) * 3;
				unsigned char b_val = (unsigned char)(brightness * 255 < 255 ? brightness * 255 : 255);
                image[index + 0] = b_val;
                image[index + 1] = b_val;
                image[index + 2] = b_val;
			}
		}

		char filename[64];
		sprintf(filename, "frames/frame_%05d.png", jindex);

		stbi_write_png(filename, WindowSize, WindowSize, 3,
					image, WindowSize * 3);

		printf("Wrote %s\n", filename);
	}

    system("ffmpeg -framerate 6 "
           "-i frames/frame_%05d.png "
           "-c:v libx264 -pix_fmt yuv420p output.mp4");
           
	free(image);

    return 0;
;}

