#include "EIRENE.h"

#include <cmath>
#include <iomanip>
#include <iostream>

int main()
{
    EIRENE potential(0, 989);
    EIRENE cross_section(1, 1334);

    const double epsilon = potential.data(0);
    const double r_minimum = potential.data(3);
    const double r_zero = potential.data(4);
    const bool potential_ok =
        std::abs(potential.potential(r_minimum) + epsilon) < 1.e-10 &&
        std::abs(potential.potential(r_zero)) < 5.e-4;
    const bool limits_ok =
        std::abs(cross_section.data(15) - 1.5598) < 1.e-8 &&
        std::abs(cross_section.data(16) - 61.8164) < 1.e-8;

    bool cross_sections_ok = true;
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "V(rm)_eV " << potential.potential(r_minimum) << '\n';
    std::cout << "V(r0)_eV " << potential.potential(r_zero) << '\n';
    std::cout << "ELABMIN_eV " << cross_section.data(15) << '\n';
    std::cout << "ELABMAX_eV " << cross_section.data(16) << '\n';
    for (double energy : {0.1, 1.5598, 10.0, 61.8164, 1000.0})
    {
        const double sigma = cross_section.cal(0., energy);
        std::cout << "sigma_m2 " << energy << ' ' << sigma << '\n';
        cross_sections_ok = cross_sections_ok && sigma > 0. && std::isfinite(sigma);
    }

    if (!potential_ok || !limits_ok || !cross_sections_ok)
    {
        std::cerr << "EIRENE H.0/H.1 elastic-data check failed\n";
        return 1;
    }
    return 0;
}
