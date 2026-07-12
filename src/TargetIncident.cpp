#include "TargetIncident.h"

#include "Global.h"
#include "Particle.h"

#include <algorithm>
#include <array>
#include <stdexcept>

Tools::IncidentFluxSample SampleDTargetIncidentFlux(
    int target, double xi_normal,
    double xi_gaussian_radius, double xi_gaussian_angle)
{
    if (target < 0 || target >= 2 * N_radial)
        throw std::out_of_range("D target incident sample has invalid target index");

    const int i = target < N_radial ? 1 : poloidalLastIndex();
    const int j = target < N_radial ? target : target - N_radial;
    const auto tangent = Grid4.InwardTargetTangent(target);
    const std::array<double, 3> drift{
        ua_D_1[i][j] * B[i][j][0],
        ua_D_1[i][j] * B[i][j][1],
        ua_D_1[i][j] * B[i][j][2]};
    return Tools::SampleIncidentFlux(
        Ti[i][j], Te[i][j], Dmass, drift,
        tangent.first, tangent.second, DTargetSheathFactor,
        xi_normal, xi_gaussian_radius, xi_gaussian_angle);
}

double DTargetFastReflectionProbability(
    const Tools::IncidentFluxSample &incident)
{
    if (incident.energy_eV < DWTrimERMIN)
        return 0.;
    return std::min(
        std::max(0., coeff_recyc_target),
        D_W_Trim.ReflectionProbability(
            incident.energy_eV, incident.angle_deg));
}
