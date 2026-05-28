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
