#include "gaussian_noise.h"


#ifdef PAR_MESOYN
Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean, Real stddev)
: noise(0), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary}
{
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MT19937);
  curandSetPseudoRandomGeneratorSeed(gen, rand());
}

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean_, Real stddev_, size_t seed)
: noise(0), prng(seed), dist(mean_, stddev_), boundary{boundary}
{
  mean = mean_; stddev=stddev_;
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MT19937);
  curandSetPseudoRandomGeneratorSeed(gen, seed);
}

Gaussian_noise::~Gaussian_noise() {
  curandDestroyGenerator(gen);
}

int Gaussian_noise::generate(size_t system_size) {
  curandGenerateNormalDouble(gen, thrust::raw_pointer_cast(noise.data()), system_size, mean, stddev);

  //boundary->update_boundaries(noise);

  return 0;
}
#else
Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean, Real stddev)
: noise(0), prng{std::random_device{}()}, dist(mean, stddev), boundary{boundary}{}

Gaussian_noise::Gaussian_noise(shared_ptr<Boundary1D> boundary, Real mean, Real stddev, size_t seed)
: noise(0), prng(seed), dist(mean, stddev), boundary{boundary}
{}

Gaussian_noise::~Gaussian_noise() {}

int Gaussian_noise::generate(size_t system_size) {
  stl::host_vector<Real> tmp_noise(system_size);
  for (Real& value : tmp_noise)
    value = dist(prng);

  noise = tmp_noise;

  //boundary->update_boundaries(noise);

  return 0;
}
#endif

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