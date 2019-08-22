#ifndef GAUSSIAN_NOISE_H
#define GAUSSIAN_NOISE_H

#include <memory>
#include <random>
#include "lattice_object.h"
#include "boundary_conditions.h"
#ifdef PAR_MESODYN
  #include <cstdlib>
  #include <thrust/transform.h>
  #include <thrust/device_vector.h>
  #include <thrust/device_ptr.h>
  #include <cuda.h>
  #include <curand.h>
#endif


class Gaussian_noise {
  //Makes sure that we keep generating new numbers, instead of the same over and over.

public:
  Gaussian_noise(shared_ptr<Boundary1D>, Real, Real, size_t); // Not seeded (32 bits of randomness)
  Gaussian_noise(shared_ptr<Boundary1D>, Real, Real, size_t, size_t); // Seeded
  ~Gaussian_noise();
  int generate(size_t);
  int add_noise(stl::device_vector<Real>&) const;
  int add_noise(Lattice_object<Real>&) const;
  stl::device_vector<Real> noise;

  private:
  const seed_seq seed;
  std::minstd_rand  prng;
  std::normal_distribution<Real> dist;
  Real variance;
#ifdef PAR_MESODYN
  curandGenerator_t gen;
#endif
  shared_ptr<Boundary1D> boundary;
};

#endif