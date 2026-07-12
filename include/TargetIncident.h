#ifndef TARGET_INCIDENT_H
#define TARGET_INCIDENT_H

#include "utils.h"

Tools::IncidentFluxSample SampleDTargetIncidentFlux(
    int target, double xi_normal,
    double xi_gaussian_radius, double xi_gaussian_angle);
Tools::IncidentFluxSample SampleDIncidentFluxAtSurface(
    int grid_i, int grid_j,
    double inward_tangent_cos, double inward_tangent_sin,
    double xi_normal, double xi_gaussian_radius,
    double xi_gaussian_angle);
Tools::IncidentFluxSample SampleDPlasmaBoundaryOutflow(
    int grid_i, int grid_j,
    double inward_tangent_cos, double inward_tangent_sin,
    double xi_normal, double xi_gaussian_radius,
    double xi_gaussian_angle);
double DFastReflectionProbability(
    const Tools::IncidentFluxSample &incident,
    double recycling_coefficient);
double DTargetFastReflectionProbability(
    const Tools::IncidentFluxSample &incident);

#endif
