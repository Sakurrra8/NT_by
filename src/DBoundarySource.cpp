#include "DBoundarySource.h"

#include "Global.h"
#include "SOLPSField.h"
#include "TargetIncident.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

DBoundarySourceModel D_BoundarySource;

namespace
{
constexpr double two_pi = 6.28318530717958647692;

double RadicalInverse(unsigned int index, unsigned int base)
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

void BuildCdf(const std::vector<DBoundarySourceSegment> &segments,
              bool molecular, std::vector<double> &cdf, double &sum)
{
    cdf.resize(segments.size());
    sum = 0.;
    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const double source = molecular ? segments[i].d2_source : segments[i].d_source;
        if (source > 0.)
            sum += source;
        cdf[i] = sum;
    }
}

int SampleCdf(const std::vector<double> &cdf, double sum)
{
    if (cdf.empty() || !(sum > 0.))
        throw std::runtime_error("Attempted to sample an empty D boundary source");
    const double value = Tools::Random() * sum;
    const auto found = std::lower_bound(cdf.begin(), cdf.end(), value);
    return found == cdf.end()
               ? static_cast<int>(cdf.size() - 1)
               : static_cast<int>(found - cdf.begin());
}

int RegionIndex(DBoundarySourceRegion region)
{
    return static_cast<int>(region);
}
}

const char *DBoundarySourceRegionName(DBoundarySourceRegion region)
{
    switch (region)
    {
    case DBoundarySourceRegion::PFRSide1:
        return "PFRSide1";
    case DBoundarySourceRegion::PFRSide2:
        return "PFRSide2";
    case DBoundarySourceRegion::OuterSide:
        return "OuterSide";
    }
    return "Unknown";
}

void DBoundarySourceModel::Load(const std::string &path)
{
    const SOLPSIndexedField field =
        ReadSOLPSIndexedField(path, "Radial particle flux");
    radial_flux_.assign(
        N_poloidal, std::vector<double>(N_radial, 0.));
    std::vector<std::vector<bool>> seen(
        N_poloidal, std::vector<bool>(N_radial, false));

    for (const SOLPSIndexedValue &row : field.values)
    {
        const int i = row.i + 1;
        const int j = row.j + 1;
        if (i < 0 || i >= N_poloidal || j < 0 || j >= N_radial)
            throw std::runtime_error("FNIY index is outside the configured B2 grid");
        if (row.component != 0)
            throw std::runtime_error("FNIY file contains an unexpected component index");
        if (seen[i][j])
            throw std::runtime_error("FNIY file contains a duplicate B2 index");
        seen[i][j] = true;
        radial_flux_[i][j] = row.value;
    }

    if (field.values.size() != static_cast<std::size_t>(gridCellCount()))
    {
        std::ostringstream message;
        message << "FNIY row count " << field.values.size()
                << " does not match B2 grid size " << gridCellCount();
        throw std::runtime_error(message.str());
    }
    for (int i = 0; i < N_poloidal; ++i)
        for (int j = 0; j < N_radial; ++j)
            if (!seen[i][j])
                throw std::runtime_error("FNIY file does not cover every B2 index");

    loaded_ = true;
    std::cout << "Loaded D+ radial boundary flux: " << path << std::endl;
}

void DBoundarySourceModel::Prepare()
{
    if (!loaded_)
        throw std::runtime_error("D boundary source was enabled but FNIY was not loaded");
    if (K_DWTrimReflection != 1 || K_Wallelement != 1)
        throw std::runtime_error(
            "D boundary source currently requires D-on-W TRIM reflection");

    using Key = std::pair<int, int>;
    std::map<Key, std::vector<DBoundarySourceSegment>> candidates;
    for (std::size_t triangle = 0; triangle < Grid4.tris_.size(); ++triangle)
    {
        if (Grid4.if_in_plasmagrid(static_cast<int>(triangle)) != 1)
            continue;
        const int grid_i = Grid4.b2_index(static_cast<int>(triangle), 0);
        const int grid_j = Grid4.b2_index(static_cast<int>(triangle), 1);

        DBoundarySourceRegion region;
        bool selected = false;
        if (grid_j == 1 && grid_i >= 1 && grid_i <= 24)
        {
            region = DBoundarySourceRegion::PFRSide1;
            selected = true;
        }
        else if (grid_j == 1 && grid_i >= 73 && grid_i <= poloidalLastIndex())
        {
            region = DBoundarySourceRegion::PFRSide2;
            selected = true;
        }
        else if (grid_j == radialLastIndex() &&
                 grid_i >= 1 && grid_i <= poloidalLastIndex())
        {
            region = DBoundarySourceRegion::OuterSide;
            selected = true;
        }
        if (!selected)
            continue;

        const auto &tri = Grid4.tris_[triangle];
        for (int edge = 0; edge < 3; ++edge)
        {
            if (Grid4.lines_info(static_cast<int>(triangle), edge) != 2)
                continue;
            const auto &a = Grid4.nodes_[tri.v[edge]];
            const auto &b = Grid4.nodes_[tri.v[(edge + 1) % 3]];
            const double dx = b.r - a.r;
            const double dy = b.z - a.z;
            const double length = std::hypot(dx, dy);
            if (!(length > 0.))
                continue;

            DBoundarySourceSegment segment;
            segment.region = region;
            segment.grid_i = grid_i;
            segment.grid_j = grid_j;
            segment.file_i = grid_i - 1;
            segment.file_j = region == DBoundarySourceRegion::OuterSide ? 36 : -1;
            segment.triangle = static_cast<int>(triangle);
            segment.edge = edge;
            segment.x0 = a.r;
            segment.y0 = a.z;
            segment.x1 = b.r;
            segment.y1 = b.z;
            segment.tangent_cos = dx / length;
            segment.tangent_sin = dy / length;

            const double mid_x = 0.5 * (a.r + b.r);
            const double mid_y = 0.5 * (a.z + b.z);
            const double centroid_x =
                (Grid4.nodes_[tri.v[0]].r + Grid4.nodes_[tri.v[1]].r +
                 Grid4.nodes_[tri.v[2]].r) /
                3.;
            const double centroid_y =
                (Grid4.nodes_[tri.v[0]].z + Grid4.nodes_[tri.v[1]].z +
                 Grid4.nodes_[tri.v[2]].z) /
                3.;
            const double inward_x = -segment.tangent_sin;
            const double inward_y = segment.tangent_cos;
            if ((centroid_x - mid_x) * inward_x +
                    (centroid_y - mid_y) * inward_y <
                0.)
            {
                segment.tangent_cos = -segment.tangent_cos;
                segment.tangent_sin = -segment.tangent_sin;
            }
            segment.surface_area = two_pi * std::max(0., mid_x) * length;
            candidates[{RegionIndex(region), grid_i}].push_back(segment);
        }
    }

    segments_.clear();
    ion_flux_sum_ = 0.;
    const auto append_faces = [&](DBoundarySourceRegion region,
                                  int first_i, int last_i,
                                  int flux_j, double outward_sign)
    {
        for (int grid_i = first_i; grid_i <= last_i; ++grid_i)
        {
            const double face_flux =
                std::max(0., outward_sign * radial_flux_[grid_i][flux_j]);
            if (!(face_flux > 0.))
                continue;
            const auto found = candidates.find({RegionIndex(region), grid_i});
            if (found == candidates.end() || found->second.empty())
            {
                std::ostringstream message;
                message << "No plasma-boundary triangle edge for "
                        << DBoundarySourceRegionName(region)
                        << " B2 i=" << grid_i << " with outward FNIY=" << face_flux;
                throw std::runtime_error(message.str());
            }

            double area_sum = 0.;
            for (const auto &segment : found->second)
                area_sum += segment.surface_area;
            if (!(area_sum > 0.))
                throw std::runtime_error("D boundary source has zero surface area");
            for (auto segment : found->second)
            {
                segment.ion_flux = face_flux * segment.surface_area / area_sum;
                segments_.push_back(segment);
            }
            ion_flux_sum_ += face_flux;
        }
    };

    append_faces(DBoundarySourceRegion::PFRSide1, 1, 24, 0, -1.);
    append_faces(DBoundarySourceRegion::PFRSide2, 73, poloidalLastIndex(), 0, -1.);
    append_faces(DBoundarySourceRegion::OuterSide, 1, poloidalLastIndex(),
                 N_radial - 1, 1.);
    if (segments_.empty() || !(ion_flux_sum_ > 0.))
        throw std::runtime_error("D boundary source has no positive outward FNIY");

    std::array<double, 3> region_ion{};
    std::array<double, 3> region_d{};
    std::array<double, 3> region_d2{};
    double weighted_energy = 0.;
    double weighted_angle = 0.;
    double weighted_fast = 0.;
    const double recycling = std::max(0., coeff_ercyc_wall);
    for (auto &segment : segments_)
    {
        double energy_sum = 0.;
        double angle_sum = 0.;
        double fast_sum = 0.;
        double reflected_energy_sum = 0.;
        if (DBoundaryLaunchModel == 1)
        {
            for (int sample_index = 1;
                 sample_index <= DTargetIncidentSamples; ++sample_index)
            {
                const auto incident = SampleDPlasmaBoundaryOutflow(
                    segment.grid_i, segment.grid_j,
                    segment.tangent_cos, segment.tangent_sin,
                    RadicalInverse(sample_index, 2),
                    RadicalInverse(sample_index, 3),
                    RadicalInverse(sample_index, 5));
                energy_sum += incident.energy_eV;
                angle_sum += incident.angle_deg;
            }
        }
        else if (DTargetIncidentModel == 0)
        {
            Tools::IncidentFluxSample incident;
            incident.energy_eV =
                2. * Ti[segment.grid_i][segment.grid_j] +
                3. * Te[segment.grid_i][segment.grid_j];
            incident.angle_deg = Tools::CalBFieldToWallNormalAngle(
                B[segment.grid_i][segment.grid_j][0],
                B[segment.grid_i][segment.grid_j][1],
                B[segment.grid_i][segment.grid_j][2],
                segment.tangent_cos, segment.tangent_sin);
            const double fast =
                DFastReflectionProbability(incident, recycling);
            energy_sum = incident.energy_eV * DTargetIncidentSamples;
            angle_sum = incident.angle_deg * DTargetIncidentSamples;
            fast_sum = fast * DTargetIncidentSamples;
            if (fast > 0.)
                reflected_energy_sum = fast * DTargetIncidentSamples *
                    D_W_Trim.MeanReflectedEnergy(
                        incident.energy_eV, incident.angle_deg);
        }
        else
        {
            for (int sample_index = 1;
                 sample_index <= DTargetIncidentSamples; ++sample_index)
            {
                const auto incident = SampleDIncidentFluxAtSurface(
                    segment.grid_i, segment.grid_j,
                    segment.tangent_cos, segment.tangent_sin,
                    RadicalInverse(sample_index, 2),
                    RadicalInverse(sample_index, 3),
                    RadicalInverse(sample_index, 5));
                const double fast =
                    DFastReflectionProbability(incident, recycling);
                energy_sum += incident.energy_eV;
                angle_sum += incident.angle_deg;
                fast_sum += fast;
                if (fast > 0.)
                    reflected_energy_sum += fast *
                        D_W_Trim.MeanReflectedEnergy(
                            incident.energy_eV, incident.angle_deg);
            }
        }
        const double inverse_samples = 1. / DTargetIncidentSamples;
        segment.mean_incident_energy = energy_sum * inverse_samples;
        segment.mean_incident_angle = angle_sum * inverse_samples;
        segment.fast_probability = fast_sum * inverse_samples;
        if (fast_sum > 0.)
            segment.mean_reflected_energy = reflected_energy_sum / fast_sum;
        if (DBoundaryLaunchModel == 1)
        {
            segment.d_source = segment.ion_flux;
            segment.d2_source = 0.;
        }
        else
        {
            segment.d_source = segment.ion_flux * segment.fast_probability;
            segment.d2_source =
                segment.ion_flux *
                std::max(0., recycling - segment.fast_probability) / 2.;
        }

        const int region = RegionIndex(segment.region);
        region_ion[region] += segment.ion_flux;
        region_d[region] += segment.d_source;
        region_d2[region] += segment.d2_source;
        weighted_energy += segment.ion_flux * segment.mean_incident_energy;
        weighted_angle += segment.ion_flux * segment.mean_incident_angle;
        weighted_fast += segment.ion_flux * segment.fast_probability;
    }

    BuildCdf(segments_, false, d_cdf_, d_source_sum_);
    BuildCdf(segments_, true, d2_cdf_, d2_source_sum_);
    std::cout << std::scientific << std::setprecision(8)
              << "D boundary source (launch model=" << DBoundaryLaunchModel
              << ", incident model=" << DTargetIncidentModel
              << "): segments=" << segments_.size()
              << ", ion=" << ion_flux_sum_
              << ", D=" << d_source_sum_
              << ", D2=" << d2_source_sum_
              << ", D nuclei=" << d_source_sum_ + 2. * d2_source_sum_
              << std::endl;
    for (int region = 0; region < 3; ++region)
        std::cout << "  "
                  << DBoundarySourceRegionName(
                         static_cast<DBoundarySourceRegion>(region))
                  << ": ion=" << region_ion[region]
                  << ", D=" << region_d[region]
                  << ", D2=" << region_d2[region] << std::endl;
    std::cout << "  incident averages: E="
              << weighted_energy / ion_flux_sum_
              << " eV, angle=" << weighted_angle / ion_flux_sum_
              << " deg, fast reflection=" << weighted_fast / ion_flux_sum_
              << std::defaultfloat << std::endl;
}

void DBoundarySourceModel::WriteSummary(const std::string &path) const
{
    if (!loaded_)
        return;
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("Cannot write D boundary source summary: " + path);
    output << "region,grid_i,grid_j,file_i,file_j,triangle,edge,x0,y0,x1,y1,"
              "surface_area_m2,ion_flux_s-1,mean_incident_energy_eV,"
              "mean_incident_angle_deg,fast_probability,"
              "mean_reflected_energy_eV,D_source_s-1,D2_source_s-1\n";
    output << std::scientific << std::setprecision(12);
    for (const auto &segment : segments_)
        output << DBoundarySourceRegionName(segment.region) << ','
               << segment.grid_i << ',' << segment.grid_j << ','
               << segment.file_i << ',' << segment.file_j << ','
               << segment.triangle << ',' << segment.edge << ','
               << segment.x0 << ',' << segment.y0 << ','
               << segment.x1 << ',' << segment.y1 << ','
               << segment.surface_area << ',' << segment.ion_flux << ','
               << segment.mean_incident_energy << ','
               << segment.mean_incident_angle << ','
               << segment.fast_probability << ','
               << segment.mean_reflected_energy << ','
               << segment.d_source << ',' << segment.d2_source << '\n';
}

const DBoundarySourceSegment &DBoundarySourceModel::Segment(
    std::size_t index) const
{
    if (index >= segments_.size())
        throw std::out_of_range("D boundary source segment index is invalid");
    return segments_[index];
}

int DBoundarySourceModel::SampleDSource() const
{
    return SampleCdf(d_cdf_, d_source_sum_);
}

int DBoundarySourceModel::SampleD2Source() const
{
    return SampleCdf(d2_cdf_, d2_source_sum_);
}

double DBoundarySourceModel::DSourceSum() const { return d_source_sum_; }
double DBoundarySourceModel::D2SourceSum() const { return d2_source_sum_; }
double DBoundarySourceModel::IonFluxSum() const { return ion_flux_sum_; }
std::size_t DBoundarySourceModel::SegmentCount() const { return segments_.size(); }
