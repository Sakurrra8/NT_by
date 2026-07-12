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
    return 0;
}
