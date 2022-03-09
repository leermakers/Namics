#include "perturbation.h"
#include "value_index_pair.h"

Register_class<Range_perturbation, Step_perturbation, std::string, const Range&, Real, Real> Step_perturbation("step");
Register_class<Range_perturbation, Sine_perturbation, std::string, const Range&, Real, Real> Sine_perturbation("sine");
Register_class<Range_perturbation, Spatial_wave_perturbation, std::string, const Range&, Real, Real> Spatial_wave_perturbation("spatial_sine");

IPerturbation::IPerturbation()
{}

IPerturbation::~IPerturbation()
{}

Range_perturbation::Range_perturbation(const Range& location_)
: IPerturbation(), m_range{location_}, m_started{true}
{}

Range_perturbation::~Range_perturbation()
{}

Range_perturbation::Range_perturbation(const Range_perturbation& other)
: m_range{other.m_range}, m_started{other.m_started}
{}

void Range_perturbation::start()
{
    m_started = true;
}

void Range_perturbation::stop()
{
    m_started = false;
}

Range Range_perturbation::get_range()
{
    return m_range;
}

Step_perturbation::Step_perturbation(const Range& location_, Real intensity_, Real dummy_)
: Range_perturbation(location_), m_intensity{intensity_}
{}

Step_perturbation::~Step_perturbation()
{}

Step_perturbation::Step_perturbation(const Step_perturbation& other)
: Range_perturbation(other.m_range), m_intensity{other.m_intensity}
{}

void Step_perturbation::perturb(Lattice_object<Real>& object_)
{
    if (m_started)
    {
        Value_index_pair<Real> pair(object_.m_data, m_range.get_indices(object_));
        Real intensity = m_intensity; //workaround for template deduction errors

        stl::for_each(pair.begin(), pair.end(), [intensity] DEVICE_LAMBDA (Real& x) mutable { x += intensity; });
    }
}

void Step_perturbation::next() {
    m_started ? stop() : start();
}

Sine_perturbation::Sine_perturbation(const Range& location_, Real amplitude, Real wavelength)
: Range_perturbation(location_), m_amplitude{amplitude}, m_wavelength{wavelength}, m_current_step{0}
{}

void Sine_perturbation::next() {
    ++m_current_step;
}

void Sine_perturbation::perturb(Lattice_object<Real>& object_)
{
    if (m_started)
    {
        Value_index_pair<Real> pair(object_.m_data, m_range.get_indices(object_));
        Real intensity = m_amplitude*sin(2*PIE*(1/m_wavelength)*m_current_step);
        stl::for_each(pair.begin(), pair.end(), [intensity] DEVICE_LAMBDA (Real& x) mutable { x += intensity; });
    }
}

Spatial_wave_perturbation::Spatial_wave_perturbation(const Range& location_, Real amplitude, Real wavelength)
: Range_perturbation(location_), m_amplitude{amplitude}, m_wavelength{wavelength}, m_phase_x{}, m_phase_y{}, m_phase_z{}, m_current_step{0}
{
    const Coordinate size = location_.as_box();

    const size_t half_x_size = size[Dimension::X];
    const size_t half_y_size = size[Dimension::Y];
    const size_t half_z_size = size[Dimension::Z];

    m_phase_x.resize(half_x_size);
    m_phase_y.resize(half_y_size);
    m_phase_z.resize(half_z_size);

    std::for_each(m_phase_x.begin(),m_phase_x.end(), [half_x_size](Real& a) {a += rand() % (half_x_size) + 1;} );
    std::for_each(m_phase_y.begin(),m_phase_y.end(), [half_y_size](Real& a) {a += rand() % (half_y_size) + 1;} );
    std::for_each(m_phase_z.begin(),m_phase_z.end(), [half_z_size](Real& a) {a += rand() % (half_z_size) + 1;} );

    m_range.loop([&](size_t x, size_t y, size_t z) mutable {

            Real value = m_amplitude*sin(2*PIE*(1/(m_wavelength))*x)*sin(2*PIE*(1/(m_wavelength))*y)*sin(2*PIE*(1/(m_wavelength))*z);

            for (size_t i{2} ; i < half_x_size ; ++i)
            {
                value += m_amplitude*sin(2*PIE*(1/(m_wavelength/i))*x+m_phase_x[i]);
            }

            for (size_t i{2}; i < half_y_size; ++i)
            {
                value += m_amplitude*sin(2*PIE*(1/(m_wavelength/i))*y+m_phase_y[i]);
            }

            for (size_t i{2}; i < half_z_size; ++i)
            {
                value += m_amplitude*sin(2*PIE*(1/(m_wavelength/i))*z+m_phase_z[i]);
            }

            m_perturbation.push_back(value);
    });
}

void Spatial_wave_perturbation::next() {
    ++m_current_step;
}

void Spatial_wave_perturbation::perturb(Lattice_object<Real>& object_)
{
    if (m_started)
    {
        Value_index_pair<Real> pair(object_.m_data, m_range.get_indices(object_) );
        stl::transform( m_perturbation.begin(), m_perturbation.end(), pair.begin(), pair.begin(), stl::plus<Real>());
    }
}

Gaussian_perturbation::Gaussian_perturbation(const Range& location_, std::shared_ptr<Gaussian_noise> noise_)
: Range_perturbation(location_), m_noise{noise_}
{}

Gaussian_perturbation::~Gaussian_perturbation()
{}

Gaussian_perturbation::Gaussian_perturbation(const Gaussian_perturbation& other)
: Range_perturbation(other.m_range), m_noise{other.m_noise}
{}

void Gaussian_perturbation::perturb(Lattice_object<Real>& object_)
{
    if (m_started)
    {
        Value_index_pair<Real> noise_range(m_noise->noise, m_range.get_indices(object_));
        Value_index_pair<Real> target_range(object_.m_data, m_range.get_indices(object_));

        stl::transform(noise_range.begin(), noise_range.end(), target_range.begin(), target_range.begin(), stl::plus<Real>() );
    }
}

void Gaussian_perturbation::next() {
    m_noise->generate(m_range.size());
}
