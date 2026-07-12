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
    return SampleDIncidentFluxAtSurface(
        i, j, tangent.first, tangent.second,
        xi_normal, xi_gaussian_radius, xi_gaussian_angle);
}

Tools::IncidentFluxSample SampleDIncidentFluxAtSurface(
    int grid_i, int grid_j,
    double inward_tangent_cos, double inward_tangent_sin,
    double xi_normal, double xi_gaussian_radius,
    double xi_gaussian_angle)
{
    if (grid_i < 0 || grid_i >= N_poloidal ||
        grid_j < 0 || grid_j >= N_radial)
        throw std::out_of_range("D incident sample has invalid B2 index");

    const std::array<double, 3> drift{
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][0],
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][1],
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][2]};
    return Tools::SampleIncidentFlux(
        Ti[grid_i][grid_j], Te[grid_i][grid_j], Dmass, drift,
        inward_tangent_cos, inward_tangent_sin, DTargetSheathFactor,
        xi_normal, xi_gaussian_radius, xi_gaussian_angle);
}

Tools::IncidentFluxSample SampleDPlasmaBoundaryOutflow(
    int grid_i, int grid_j,
    double inward_tangent_cos, double inward_tangent_sin,
    double xi_normal, double xi_gaussian_radius,
    double xi_gaussian_angle)
{
    if (grid_i < 0 || grid_i >= N_poloidal ||
        grid_j < 0 || grid_j >= N_radial)
        throw std::out_of_range("D boundary outflow sample has invalid B2 index");

    const std::array<double, 3> drift{
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][0],
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][1],
        ua_D_1[grid_i][grid_j] * B[grid_i][grid_j][2]};
    return Tools::SampleIncidentFlux(
        Ti[grid_i][grid_j], Te[grid_i][grid_j], Dmass, drift,
        inward_tangent_cos, inward_tangent_sin, 0.,
        xi_normal, xi_gaussian_radius, xi_gaussian_angle);
}

double DFastReflectionProbability(
    const Tools::IncidentFluxSample &incident,
    double recycling_coefficient)
{
    if (incident.energy_eV < DWTrimERMIN)
        return 0.;
    return std::min(
        std::max(0., recycling_coefficient),
        D_W_Trim.ReflectionProbability(
            incident.energy_eV, incident.angle_deg));
}

double DTargetFastReflectionProbability(
    const Tools::IncidentFluxSample &incident)
{
    return DFastReflectionProbability(incident, coeff_recyc_target);
}
