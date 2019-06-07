#ifndef GAUSSIAN_NOISE_H
#define GAUSSIAN_NOISE_H

#include <memory>
#include <random>
#include "lattice_object.h"
#include "boundary_conditions.h"
#ifdef PAR_MESODYN
  #include <thrust/random/linear_congruential_engine.h>
  #include <thrust/random/normal_distribution.h>
#endif

class Gaussian_noise {
  //Makes sure that we keep generating new numbers, instead of the same over and over.

public:
  Gaussian_noise(shared_ptr<Boundary1D>, Real, Real); // Not seeded (32 bits of randomness)
  Gaussian_noise(shared_ptr<Boundary1D>, Real, Real, size_t); // Seeded
  int generate(size_t);
  int add_noise(stl::device_vector<Real>&) const;
  int add_noise(Lattice_object<Real>&) const;
  stl::device_vector<Real> noise;

private:
  const seed_seq seed;
  stl::minstd_rand prng;
  stl::normal_distribution<Real> dist;
  shared_ptr<Boundary1D> boundary;
};

#endif