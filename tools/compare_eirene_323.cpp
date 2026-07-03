#include "EIRENE.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct SettingPaths
{
    int n_poloidal = 98;
    int n_radial = 38;
    std::string case_path;
    std::string database_path = "Inputfile/database/";
    std::string output_path = "Outputfile/test/data/";
};

SettingPaths readSettingPaths(const std::string &setting_file)
{
    SettingPaths paths;
    std::ifstream in(setting_file);
    if (!in.is_open())
        throw std::runtime_error("Cannot open setting file: " + setting_file);

    std::string line;
    while (std::getline(in, line))
    {
        const auto comment = line.find('#');
        if (comment != std::string::npos)
            line = line.substr(0, comment);
        std::istringstream fields(line);
        std::string key;
        std::string value;
        if (!(fields >> key >> value))
            continue;
        if (key == "N_poloidal")
            paths.n_poloidal = std::stoi(value);
        else if (key == "N_radial")
            paths.n_radial = std::stoi(value);
        else if (key == "Casepath")
            paths.case_path = value;
        else if (key == "Databasepath")
            paths.database_path = value;
        else if (key == "Outputpath")
            paths.output_path = value;
    }
    if (paths.case_path.empty())
        throw std::runtime_error("Casepath is missing in setting file: " + setting_file);
    return paths;
}

void setDatabaseEnvironment(const std::string &database_path)
{
#ifdef _WIN32
    _putenv_s("NT_DATABASE_PATH", database_path.c_str());
#else
    setenv("NT_DATABASE_PATH", database_path.c_str(), 1);
#endif
}

std::filesystem::path firstExistingPath(const std::vector<std::filesystem::path> &candidates)
{
    for (const auto &candidate : candidates)
    {
        if (std::filesystem::exists(candidate))
            return candidate;
    }
    std::ostringstream message;
    message << "None of these input files exist:";
    for (const auto &candidate : candidates)
        message << "\n  " << candidate.string();
    throw std::runtime_error(message.str());
}

std::vector<std::vector<double>> readB2Field(const std::filesystem::path &path, int n_poloidal, int n_radial)
{
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open B2 field: " + path.string());

    std::vector<std::vector<double>> field(n_poloidal, std::vector<double>(n_radial, 0.0));
    std::string header;
    std::getline(in, header);
    for (int cell = 0; cell < n_poloidal * n_radial; ++cell)
    {
        int a = 0;
        int b = 0;
        int c = 0;
        double value = 0.0;
        if (!(in >> a >> b >> c >> value))
            throw std::runtime_error("Unexpected B2 field format in: " + path.string());
        const int i = a + 1;
        const int j = b + 1;
        if (i >= 0 && i < n_poloidal && j >= 0 && j < n_radial)
            field[i][j] = value;
    }
    return field;
}

struct Stats
{
    long long count = 0;
    long long finite_count = 0;
    double h2_min = std::numeric_limits<double>::infinity();
    double h2_max = 0.0;
    double h3_min = std::numeric_limits<double>::infinity();
    double h3_max = 0.0;
    double ratio_min = std::numeric_limits<double>::infinity();
    double ratio_max = 0.0;
    double ratio_sum = 0.0;
    double abs_rel_sum = 0.0;
};
}

int main(int argc, char **argv)
{
    const std::string setting_file = argc > 1 ? argv[1] : "Inputfile/settingfile/setting_Trimesh_D_5.log";
    const SettingPaths paths = readSettingPaths(setting_file);
    setDatabaseEnvironment(paths.database_path);

    const std::filesystem::path case_path(paths.case_path);
    const auto ne_path = firstExistingPath({
        case_path / "2D_data" / "ne_2D.dat",
        case_path / "2D_data" / "ne_2D.data",
        case_path / "ne_2D.dat",
        case_path / "ne_2D.data",
    });
    const auto te_path = firstExistingPath({
        case_path / "2D_data" / "te_2D.dat",
        case_path / "2D_data" / "te_2D.data",
        case_path / "te_2D.dat",
        case_path / "te_2D.data",
    });
    const auto ne = readB2Field(ne_path, paths.n_poloidal, paths.n_radial);
    const auto te = readB2Field(te_path, paths.n_poloidal, paths.n_radial);

    EIRENE h2(2, 2253, "amjuel");
    EIRENE h3(3, 3841, "amjuel");

    const std::filesystem::path out_dir = std::filesystem::path(paths.output_path) / "eirene_323_compare";
    std::filesystem::create_directories(out_dir);
    const auto csv_path = out_dir / "R3_2_3_H2_vs_H3_B2.csv";
    const auto summary_path = out_dir / "summary.txt";

    std::ofstream csv(csv_path);
    if (!csv.is_open())
        throw std::runtime_error("Cannot write CSV: " + csv_path.string());
    csv << std::setprecision(12)
        << "i,j,ne_m-3,Te_eV,H2_arg1_ne_m-3,H2_arg2_Te_eV,H2_value_m3_s,"
        << "H3_arg1_projectile_energy_eV,H3_arg2_background_temperature_eV,H3_value_m3_s,"
        << "H3_over_H2,relative_difference_H3_minus_H2\n";

    Stats stats;
    for (int i = 0; i < paths.n_poloidal; ++i)
    {
        for (int j = 0; j < paths.n_radial; ++j)
        {
            const double h2_value = h2.cal(ne[i][j], te[i][j]);
            const double h3_value = h3.cal(te[i][j], te[i][j]);
            const double ratio = h2_value > 0.0 ? h3_value / h2_value : std::numeric_limits<double>::quiet_NaN();
            const double rel = h2_value > 0.0 ? (h3_value - h2_value) / h2_value : std::numeric_limits<double>::quiet_NaN();

            ++stats.count;
            stats.h2_min = std::min(stats.h2_min, h2_value);
            stats.h2_max = std::max(stats.h2_max, h2_value);
            stats.h3_min = std::min(stats.h3_min, h3_value);
            stats.h3_max = std::max(stats.h3_max, h3_value);
            if (std::isfinite(ratio))
            {
                ++stats.finite_count;
                stats.ratio_min = std::min(stats.ratio_min, ratio);
                stats.ratio_max = std::max(stats.ratio_max, ratio);
                stats.ratio_sum += ratio;
                stats.abs_rel_sum += std::abs(rel);
            }

            csv << i << ',' << j << ','
                << ne[i][j] << ',' << te[i][j] << ','
                << ne[i][j] << ',' << te[i][j] << ',' << h2_value << ','
                << te[i][j] << ',' << te[i][j] << ',' << h3_value << ','
                << ratio << ',' << rel << '\n';
        }
    }

    std::ofstream summary(summary_path);
    if (!summary.is_open())
        throw std::runtime_error("Cannot write summary: " + summary_path.string());
    summary << std::setprecision(12)
            << "setting=" << setting_file << '\n'
            << "ne_path=" << ne_path.string() << '\n'
            << "te_path=" << te_path.string() << '\n'
            << "H2=AMJUEL H.2 3.2.3 line 2253, cal(ne, Te)\n"
            << "H3=AMJUEL H.3 3.2.3 line 3841, cal(Te, Te)\n"
            << "rows=" << stats.count << '\n'
            << "finite_ratio_rows=" << stats.finite_count << '\n'
            << "H2_min=" << stats.h2_min << '\n'
            << "H2_max=" << stats.h2_max << '\n'
            << "H3_min=" << stats.h3_min << '\n'
            << "H3_max=" << stats.h3_max << '\n'
            << "ratio_min=" << stats.ratio_min << '\n'
            << "ratio_max=" << stats.ratio_max << '\n'
            << "ratio_mean=" << (stats.finite_count ? stats.ratio_sum / stats.finite_count : 0.0) << '\n'
            << "mean_abs_relative_difference=" << (stats.finite_count ? stats.abs_rel_sum / stats.finite_count : 0.0) << '\n';

    std::cout << "Wrote " << csv_path << '\n';
    std::cout << "Wrote " << summary_path << '\n';
    std::cout << "ratio_mean=" << (stats.finite_count ? stats.ratio_sum / stats.finite_count : 0.0)
              << " ratio_min=" << stats.ratio_min
              << " ratio_max=" << stats.ratio_max << '\n';
    return 0;
}
