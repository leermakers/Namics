#ifndef H_PERTURBATION
#define H_PERTURBATION

#include "lattice_accessor.h"
#include "gaussian_noise.h"
#include "lattice_object.h"
#include "stl_typedef.h"

#include "../tools.h"

#include <memory>

class Range_perturbation;

namespace Perturbation {
  typedef Factory_template<Range_perturbation, std::string, const Range&, Real, Real> Factory;
}

class IPerturbation {
    public:
        IPerturbation();
        ~IPerturbation();
        virtual void perturb(Lattice_object<Real>& object_) = 0;
        virtual void start() = 0;
        virtual void stop() = 0;
        virtual void next() = 0;

    protected:

    private:
};

class Range_perturbation : public IPerturbation {
    public:
        Range_perturbation(const Range& range_);
        Range_perturbation(const Range_perturbation& other);
        ~Range_perturbation();

        virtual void start();
        virtual void stop();
        virtual void next() = 0;

    protected:
        Range m_range;
        bool m_started;

    private:
};

class Step_perturbation : public Range_perturbation {
    public:
        Step_perturbation(const Range& range_, Real intensity_, Real dummy_);
        virtual ~Step_perturbation();
        Step_perturbation(const Step_perturbation&);
        void perturb(Lattice_object<Real>& object_);
        void next();

    protected:
        Real m_intensity;

    private:
};

class Spatial_wave_perturbation : public Range_perturbation {
    public:
        Spatial_wave_perturbation(const Range& range_, Real intensity_, Real dummy_);
        Spatial_wave_perturbation(const Spatial_wave_perturbation&);
        void perturb(Lattice_object<Real>& object_);
        void next();

    protected:
        Real m_amplitude;
        Real m_wavelength;
        std::vector<Real> m_phase_x;
        std::vector<Real> m_phase_y;
        size_t m_current_step;
        stl::device_vector<Real> m_perturbation;

    private:
};

class Sine_perturbation : public Range_perturbation {
    public:
        Sine_perturbation(const Range& range_, Real amplitude, Real wavelength);
        void perturb(Lattice_object<Real>& object_);
        void next();

    protected:
        Real m_amplitude;
        Real m_wavelength;
        size_t m_current_step;

    private:
};


class Gaussian_perturbation : public Range_perturbation {
    public:
        Gaussian_perturbation(const Range& range_, std::shared_ptr<Gaussian_noise> noise_);
        virtual ~Gaussian_perturbation();
        Gaussian_perturbation(const Gaussian_perturbation&);
        void perturb(Lattice_object<Real>& object_);
        void next();

    protected:
        std::shared_ptr<Gaussian_noise> m_noise;

    private:
};

#endif