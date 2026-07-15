#include "utils.h"

#include <array>
#include <cmath>
#include <iostream>

namespace
{
double radicalInverse(unsigned int index, unsigned int base)
{
    double value = 0.;
    double factor = 1. / static_cast<double>(base);
    while (index > 0)
    {
        value += factor * static_cast<double>(index % base);
        index /= base;
        factor /= static_cast<double>(base);
    }
    return value;
}
}

int main()
{
    constexpr int samples = 65536;
    constexpr double deuterium_mass = 3.34449469e-27;
    constexpr double ti_eV = 10.;
    constexpr double te_eV = 5.;
    constexpr double sheath_factor = 3.;
    const std::array<double, 3> drift{0., 0., 0.};

    double base_energy_sum = 0.;
    double sheath_energy_sum = 0.;
    double cosine_sum = 0.;
    double surface_fraction_sum = 0.;
    double maximum_sheath_increment_error = 0.;
    for (int index = 1; index <= samples; ++index)
    {
        const double xi_normal = radicalInverse(index, 2);
        const double xi_radius = radicalInverse(index, 3);
        const double xi_angle = radicalInverse(index, 5);
        const auto base = Tools::SampleIncidentFlux(
            ti_eV, te_eV, deuterium_mass, drift,
            1., 0., 0., xi_normal, xi_radius, xi_angle);
        const auto sheath = Tools::SampleIncidentFlux(
            ti_eV, te_eV, deuterium_mass, drift,
            1., 0., sheath_factor, xi_normal, xi_radius, xi_angle);

        const double speed = std::sqrt(
            Tools::sqr(base.velocity[0]) + Tools::sqr(base.velocity[1]) +
            Tools::sqr(base.velocity[2]));
        const double wallward_speed = -base.velocity[1];
        if (!(speed > 0.) || !(wallward_speed > 0.) ||
            !std::isfinite(base.energy_eV) || !std::isfinite(sheath.energy_eV))
        {
            std::cerr << "Invalid incident flux sample at index " << index << '\n';
            return 1;
        }
        base_energy_sum += base.energy_eV;
        sheath_energy_sum += sheath.energy_eV;
        cosine_sum += wallward_speed / speed;
        surface_fraction_sum += Tools::SampleAxisymmetricSurfaceFraction(
            1., 2., (static_cast<double>(index) - 0.5) / samples);
        maximum_sheath_increment_error = std::max(
            maximum_sheath_increment_error,
            std::abs(
                sheath.energy_eV - base.energy_eV - sheath_factor * te_eV));
    }

    const double base_mean_energy = base_energy_sum / samples;
    const double sheath_mean_energy = sheath_energy_sum / samples;
    const double mean_cosine = cosine_sum / samples;
    const double mean_surface_fraction = surface_fraction_sum / samples;
    const double expected_surface_fraction = 5. / 9.;

    constexpr double wall_temperature_eV = 0.1;
    double emitted_energy_sum = 0.;
    double emitted_cosine_sum = 0.;
    std::vector<double> direction(3, 0.);
    for (int index = 0; index < samples; ++index)
    {
        const double launch_speed = Tools::MaxwellianFluxSpeed(
            wall_temperature_eV, 2. * deuterium_mass);
        Tools::calculateReflectionVelocity(direction, 0.6, 0.8, 0);
        emitted_energy_sum +=
            0.5 * 2. * deuterium_mass * launch_speed * launch_speed / qe;
        emitted_cosine_sum += -0.8 * direction[0] + 0.6 * direction[1];
    }
    const double mean_emitted_energy = emitted_energy_sum / samples;
    const double mean_emitted_cosine = emitted_cosine_sum / samples;
    std::cout << "base_mean_energy_eV=" << base_mean_energy << '\n'
              << "expected_base_mean_energy_eV=" << 2. * ti_eV << '\n'
              << "sheath_mean_energy_eV=" << sheath_mean_energy << '\n'
              << "expected_sheath_mean_energy_eV="
              << 2. * ti_eV + sheath_factor * te_eV << '\n'
              << "mean_cosine=" << mean_cosine << '\n'
              << "expected_mean_cosine=" << 2. / 3. << '\n'
              << "mean_surface_fraction=" << mean_surface_fraction << '\n'
              << "expected_surface_fraction=" << expected_surface_fraction << '\n'
              << "mean_D2_emitted_energy_eV=" << mean_emitted_energy << '\n'
              << "expected_D2_emitted_energy_eV="
              << 2. * wall_temperature_eV << '\n'
              << "mean_D2_inward_cosine=" << mean_emitted_cosine << '\n'
              << "expected_D2_inward_cosine=" << 2. / 3. << '\n'
              << "max_sheath_increment_error_eV="
              << maximum_sheath_increment_error << '\n';

    if (std::abs(base_mean_energy - 2. * ti_eV) > 0.02 ||
        std::abs(sheath_mean_energy - (2. * ti_eV + sheath_factor * te_eV)) > 0.02 ||
        std::abs(mean_cosine - 2. / 3.) > 0.002 ||
        std::abs(mean_surface_fraction - expected_surface_fraction) > 1.e-5 ||
        std::abs(mean_emitted_energy - 2. * wall_temperature_eV) > 0.003 ||
        std::abs(mean_emitted_cosine - 2. / 3.) > 0.006 ||
        maximum_sheath_increment_error > 1.e-10)
        return 1;
    return 0;
}
