#ifndef D_BOUNDARY_SOURCE_H
#define D_BOUNDARY_SOURCE_H

#include <cstddef>
#include <string>
#include <vector>

enum class DBoundarySourceRegion
{
    PFRSide1 = 0,
    PFRSide2,
    OuterSide
};

const char *DBoundarySourceRegionName(DBoundarySourceRegion region);

struct DBoundarySourceSegment
{
    DBoundarySourceRegion region{DBoundarySourceRegion::PFRSide1};
    int grid_i{-1};
    int grid_j{-1};
    int file_i{0};
    int file_j{0};
    int triangle{-1};
    int edge{-1};
    double x0{0.};
    double y0{0.};
    double x1{0.};
    double y1{0.};
    double tangent_cos{0.};
    double tangent_sin{0.};
    double surface_area{0.};
    double ion_flux{0.};
    double mean_incident_energy{0.};
    double mean_incident_angle{0.};
    double fast_probability{0.};
    double mean_reflected_energy{0.};
    double d_source{0.};
    double d2_source{0.};
};

class DBoundarySourceModel
{
public:
    void Load(const std::string &path);
    void Prepare();
    void WriteSummary(const std::string &path) const;

    const DBoundarySourceSegment &Segment(std::size_t index) const;
    int SampleDSource() const;
    int SampleD2Source() const;
    double DSourceSum() const;
    double D2SourceSum() const;
    double IonFluxSum() const;
    std::size_t SegmentCount() const;

private:
    std::vector<std::vector<double>> radial_flux_;
    std::vector<DBoundarySourceSegment> segments_;
    std::vector<double> d_cdf_;
    std::vector<double> d2_cdf_;
    double d_source_sum_{0.};
    double d2_source_sum_{0.};
    double ion_flux_sum_{0.};
    bool loaded_{false};
};

extern DBoundarySourceModel D_BoundarySource;

#endif
