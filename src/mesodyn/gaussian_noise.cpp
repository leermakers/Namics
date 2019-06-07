#include "gaussian_noise.h"

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean, Real stddev)
: noise(0), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary}{}

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean, Real stddev, size_t seed)
: noise(0), prng(seed), dist(mean, stddev), boundary{boundary}
{}

int Gaussian_noise::generate(size_t system_size) {
  stl::host_vector<Real> tmp_noise(system_size);
  for (Real& value : tmp_noise)
    value = dist(prng);

  noise = tmp_noise;

  boundary->update_boundaries(noise);

  return 0;
}

int Gaussian_noise::add_noise(stl::device_vector<Real>& target) const{
  assert(noise.size() == target.size());
  stl::transform(noise.begin(), noise.end(), target.begin(), target.begin(), stl::plus<Real>());
  return 0;
}


int Gaussian_noise::add_noise(Lattice_object<Real>& target) const {
  assert(noise.size() == target.size());
  stl::transform(noise.begin(), noise.end(), target.begin(), target.begin(), stl::plus<Real>());
  return 0;
}