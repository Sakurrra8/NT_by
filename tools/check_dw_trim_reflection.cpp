#include "Reflect.h"

#include <algorithm>
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

int main(int argc, char **argv)
{
    const char *database = argc > 1 ? argv[1] : "Inputfile/database/D_on_W";
    DWTrimReflection reflection;
    if (!reflection.Load(database))
        return 1;

    constexpr int samples = 131072;
    constexpr double incident_energy_eV = 20.;
    constexpr double incident_angle_deg = 40.;
    constexpr double minimum_energy_eV = 0.2;
    double energy_sum = 0.;
    double minimum_sampled_energy = incident_energy_eV;
    for (int index = 1; index <= samples; ++index)
    {
        const double xi_energy = radicalInverse(index, 2);
        const double xi_polar = radicalInverse(index, 3);
        const double xi_azimuth = radicalInverse(index, 5);
        const auto sample = reflection.Sample(
            incident_energy_eV, incident_angle_deg,
            xi_energy, xi_polar, xi_azimuth, minimum_energy_eV);
        if (!std::isfinite(sample.energy_eV) ||
            sample.energy_eV < minimum_energy_eV ||
            sample.energy_eV > incident_energy_eV ||
            sample.cos_polar < 0. || sample.cos_polar > 1. ||
            sample.cos_azimuth < -1. || sample.cos_azimuth > 1.)
        {
            std::cerr << "Invalid D-on-W TRIM sample at index " << index << '\n';
            return 1;
        }
        energy_sum += sample.energy_eV;
        minimum_sampled_energy = std::min(minimum_sampled_energy, sample.energy_eV);
    }

    const double sampled_mean_energy = energy_sum / samples;
    const double tabulated_mean_energy =
        reflection.MeanReflectedEnergy(incident_energy_eV, incident_angle_deg);
    std::cout << "sampled_mean_energy_eV=" << sampled_mean_energy << '\n'
              << "tabulated_mean_energy_eV=" << tabulated_mean_energy << '\n'
              << "minimum_sampled_energy_eV=" << minimum_sampled_energy << '\n';

    if (std::abs(sampled_mean_energy - tabulated_mean_energy) > 0.02)
        return 1;

    const auto grazing_sample = reflection.Sample(
        500., 85., 1., 0., 0.5, minimum_energy_eV);
    if (grazing_sample.cos_polar < 0.08716 ||
        grazing_sample.cos_polar > 0.999999)
    {
        std::cerr << "D-on-W polar-angle sample is outside EIRENE REFLC1 bounds\n";
        return 1;
    }

    constexpr double behrisch_grid_energy_eV = 10. * 58. / 57.;
    constexpr double expected_normal_probability = 0.23256897420786862;
    const double normal_probability = EireneDFeReflection::ReflectionProbability(
        behrisch_grid_energy_eV, 0.);
    const double oblique_probability = EireneDFeReflection::ReflectionProbability(
        behrisch_grid_energy_eV, 60.);
    const double expected_oblique_probability =
        1. - 0.5 * (1. - expected_normal_probability);
    if (EireneDFeReflection::ReflectionProbability(1., 0.) != 0. ||
        std::abs(normal_probability - expected_normal_probability) > 1.e-12 ||
        std::abs(oblique_probability - expected_oblique_probability) > 1.e-12)
    {
        std::cerr << "EIRENE Behrisch reflection probability mismatch\n";
        return 1;
    }

    double behrisch_energy_sum = 0.;
    for (int index = 1; index <= samples; ++index)
    {
        const double energy = EireneDFeReflection::SampleReflectedEnergy(
            20., 35., radicalInverse(index, 2));
        if (!std::isfinite(energy) || energy < 0. || energy > 20.)
        {
            std::cerr << "Invalid EIRENE Behrisch energy sample\n";
            return 1;
        }
        behrisch_energy_sum += energy;
    }
    const double behrisch_sampled_mean = behrisch_energy_sum / samples;
    const double behrisch_expected_mean =
        EireneDFeReflection::MeanReflectedEnergy(20., 35.);
    std::cout << "behrisch_sampled_mean_energy_eV="
              << behrisch_sampled_mean << '\n'
              << "behrisch_expected_mean_energy_eV="
              << behrisch_expected_mean << '\n';
    if (std::abs(behrisch_sampled_mean - behrisch_expected_mean) > 0.002)
        return 1;
    return 0;
}
