#ifndef TARGET_INCIDENT_H
#define TARGET_INCIDENT_H

#include "utils.h"

Tools::IncidentFluxSample SampleDTargetIncidentFlux(
    int target, double xi_normal,
    double xi_gaussian_radius, double xi_gaussian_angle);
double DTargetFastReflectionProbability(
    const Tools::IncidentFluxSample &incident);

#endif
