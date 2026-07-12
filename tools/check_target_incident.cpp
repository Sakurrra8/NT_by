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
    const std::array<double, 3> wallward_drift{0., -15000., 0.};

    double base_energy_sum = 0.;
    double sheath_energy_sum = 0.;
    double cosine_sum = 0.;
    double eirene_base_energy_sum = 0.;
    double eirene_weight_sum = 0.;
    double drift_energy_sum = 0.;
    double eirene_drift_energy_sum = 0.;
    double eirene_drift_weight_sum = 0.;
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
        const auto eirene = Tools::SampleEireneIncidentFlux(
            ti_eV, te_eV, deuterium_mass, drift,
            1., 0., 0., xi_normal, xi_radius, xi_angle);
        const auto drifted = Tools::SampleIncidentFlux(
            ti_eV, te_eV, deuterium_mass, wallward_drift,
            1., 0., 0., xi_normal, xi_radius, xi_angle);
        const auto eirene_drifted = Tools::SampleEireneIncidentFlux(
            ti_eV, te_eV, deuterium_mass, wallward_drift,
            1., 0., 0., xi_normal, xi_radius, xi_angle);

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
        eirene_base_energy_sum += eirene.statistical_weight * eirene.energy_eV;
        eirene_weight_sum += eirene.statistical_weight;
        drift_energy_sum += drifted.energy_eV;
        eirene_drift_energy_sum +=
            eirene_drifted.statistical_weight * eirene_drifted.energy_eV;
        eirene_drift_weight_sum += eirene_drifted.statistical_weight;
        maximum_sheath_increment_error = std::max(
            maximum_sheath_increment_error,
            std::abs(
                sheath.energy_eV - base.energy_eV - sheath_factor * te_eV));
    }

    const double base_mean_energy = base_energy_sum / samples;
    const double sheath_mean_energy = sheath_energy_sum / samples;
    const double mean_cosine = cosine_sum / samples;
    const double eirene_base_mean_energy = eirene_base_energy_sum / eirene_weight_sum;
    const double drift_mean_energy = drift_energy_sum / samples;
    const double eirene_drift_mean_energy =
        eirene_drift_energy_sum / eirene_drift_weight_sum;
    std::cout << "base_mean_energy_eV=" << base_mean_energy << '\n'
              << "expected_base_mean_energy_eV=" << 2. * ti_eV << '\n'
              << "sheath_mean_energy_eV=" << sheath_mean_energy << '\n'
              << "expected_sheath_mean_energy_eV="
              << 2. * ti_eV + sheath_factor * te_eV << '\n'
              << "mean_cosine=" << mean_cosine << '\n'
              << "expected_mean_cosine=" << 2. / 3. << '\n'
              << "eirene_base_mean_energy_eV=" << eirene_base_mean_energy << '\n'
              << "eirene_mean_weight=" << eirene_weight_sum / samples << '\n'
              << "drift_mean_energy_eV=" << drift_mean_energy << '\n'
              << "eirene_drift_mean_energy_eV=" << eirene_drift_mean_energy << '\n'
              << "eirene_drift_mean_weight="
              << eirene_drift_weight_sum / samples << '\n'
              << "max_sheath_increment_error_eV="
              << maximum_sheath_increment_error << '\n';

    if (std::abs(base_mean_energy - 2. * ti_eV) > 0.02 ||
        std::abs(sheath_mean_energy - (2. * ti_eV + sheath_factor * te_eV)) > 0.02 ||
        std::abs(mean_cosine - 2. / 3.) > 0.002 ||
        std::abs(eirene_base_mean_energy - 2. * ti_eV) > 0.02 ||
        std::abs(eirene_weight_sum / samples - 1.) > 1.e-12 ||
        std::abs(eirene_drift_mean_energy - drift_mean_energy) > 0.05 ||
        std::abs(eirene_drift_weight_sum / samples - 1.) > 0.002 ||
        maximum_sheath_increment_error > 1.e-10)
        return 1;
    return 0;
}
