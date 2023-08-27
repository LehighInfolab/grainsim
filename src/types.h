#pragma once

// The number of neighbors a voxel can have.
const char NEIGH_COUNT = 26;

// The data type that holds a grain's spin (ID).
typedef unsigned int spin_t;
// The data type that holds an activity (probability).
typedef double activ_t;
// The data type that holds a single dimension within the lattice (must be signed to allow for wrapping).
typedef int coord_t;