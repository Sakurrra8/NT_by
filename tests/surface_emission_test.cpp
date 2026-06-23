#include "utils.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
    constexpr double temperature_eV = 0.1;
    constexpr double deuterium_molecule_mass = 6.690487694e-27;
    constexpr int samples = 200000;

    double energy_sum_eV = 0.0;
    for (int i = 0; i < samples; ++i)
    {
        const double speed =
            Tools::MaxwellianFluxSpeed(temperature_eV, deuterium_molecule_mass);
        assert(std::isfinite(speed));
        assert(speed > 0.0);
        energy_sum_eV +=
            0.5 * deuterium_molecule_mass * speed * speed / qe;
    }

    const double mean_energy_eV = energy_sum_eV / samples;
    const double expected_mean_energy_eV = 2.0 * temperature_eV;
    assert(std::abs(mean_energy_eV - expected_mean_energy_eV) <
           0.03 * expected_mean_energy_eV);

    std::cout << "mean thermal-flux energy = " << mean_energy_eV << " eV\n";
    return 0;
}
