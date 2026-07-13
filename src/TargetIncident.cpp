#include "TargetIncident.h"

#include "Global.h"
#include "Particle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

double DTargetSheathFactorAtCell(int grid_i, int grid_j)
{
    if (grid_i < 0 || grid_i >= N_poloidal ||
        grid_j < 0 || grid_j >= N_radial)
        throw std::out_of_range("D target sheath factor has invalid B2 index");
    if (DTargetSheathFactor > 0.)
        return DTargetSheathFactor;

    constexpr double electron_mass = 5.448e-4 * 1.6606e-27;
    constexpr double two_pi = 6.28318530717958647692;
    const double electron_temperature = Te[grid_i][grid_j];
    const double ion_speed = std::abs(ua_D_1[grid_i][grid_j]);
    if (!(electron_temperature > 0.) || !(ion_speed > 0.))
        return 2.8;

    const double electron_thermal_speed =
        std::sqrt(qe * electron_temperature / electron_mass);
    const double sheath_argument =
        std::sqrt(two_pi) * ion_speed / electron_thermal_speed;
    if (!(sheath_argument > 0.) || !std::isfinite(sheath_argument))
        return 2.8;
    return std::max(0., -std::log(sheath_argument));
}

double EireneTargetIncidentEnergy(int target)
{
    if (target < 0 || target >= 2 * N_radial)
        throw std::out_of_range("D target incident energy has invalid target index");

    const int grid_i = target < N_radial ? 0 : N_poloidal - 1;
    const int grid_j = target < N_radial ? target : target - N_radial;
    // INFCOP assigns 3 Ti + 0.5 Te to the main poloidal target flux.
    return std::max(
        0., 3. * Ti[grid_i][grid_j] + 0.5 * Te[grid_i][grid_j] +
                DTargetSheathFactorAtCell(grid_i, grid_j) *
                    Te[grid_i][grid_j]);
}

double EireneRadialBoundaryIncidentEnergy(int grid_i, int grid_j)
{
    // INFCOP assigns 2 Ti to the main radial PFR/outer-boundary flux.
    return std::max(
        0., 2. * Ti[grid_i][grid_j] +
                DTargetSheathFactorAtCell(grid_i, grid_j) *
                    Te[grid_i][grid_j]);
}

Tools::IncidentFluxSample SampleDTargetIncidentFlux(
    int target, double xi_normal,
    double xi_gaussian_radius, double xi_gaussian_angle)
{
    if (target < 0 || target >= 2 * N_radial)
        throw std::out_of_range("D target incident sample has invalid target index");

    const int i = target < N_radial ? 0 : N_poloidal - 1;
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
        inward_tangent_cos, inward_tangent_sin,
        DTargetSheathFactorAtCell(grid_i, grid_j),
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
