#include "EIRENE.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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

struct ReactionSpec
{
    std::string name;
    std::string database;
    int fit;
    int line;
    std::string eirene_table;
    std::string reaction;
    std::string argument1_name;
    std::string argument2_name;
};

std::string trim(const std::string &text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

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

std::string csvEscape(const std::string &text)
{
    if (text.find_first_of(",\"\n\r") == std::string::npos)
        return text;
    std::string escaped = "\"";
    for (const char c : text)
    {
        if (c == '"')
            escaped += "\"\"";
        else
            escaped += c;
    }
    escaped += '"';
    return escaped;
}

std::vector<double> logGrid(double min_value, double max_value, int count)
{
    std::vector<double> values;
    values.reserve(count);
    const double log_min = std::log10(min_value);
    const double log_max = std::log10(max_value);
    for (int i = 0; i < count; ++i)
    {
        const double fraction = count == 1 ? 0.0 : static_cast<double>(i) / (count - 1);
        values.push_back(std::pow(10.0, log_min + fraction * (log_max - log_min)));
    }
    return values;
}

std::string sanitizeFileName(const std::string &text)
{
    std::string name = text;
    for (char &c : name)
    {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-' && c != '.')
            c = '_';
    }
    return name;
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

    std::vector<std::vector<double>> field(
        n_poloidal, std::vector<double>(n_radial, 0.0));
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

std::vector<ReactionSpec> reactionSpecs()
{
    const std::string ne = "electron_density_m-3";
    const std::string te = "electron_temperature_eV";
    const std::string eproj = "projectile_energy_eV";
    const std::string tback = "background_temperature_eV";
    return {
        {"R3_2_3r_H4", "amjuel", 4, 6710, "H.4 3.2.3r", "MAR via H2+", ne, te},
        {"R3_2_3d_H4", "amjuel", 4, 6767, "H.4 3.2.3d", "MAD via H2+", ne, te},
        {"R3_2_3i_H4", "amjuel", 4, 6827, "H.4 3.2.3i", "MAI via H2+", ne, te},
        {"R2_2_17r_H4", "amjuel", 4, 7007, "H.4 2.2.17r", "MAR via H-", ne, te},
        {"R2_2_17d_H4", "amjuel", 4, 7061, "H.4 2.2.17d", "MAD via H-", ne, te},
        {"R2_2_5_H4", "amjuel", 4, 4511, "H.4 2.2.5", "e + H2 -> e + H + H", ne, te},
        {"R2_2_14_H2", "amjuel", 2, 2157, "H.2 2.2.14", "e + H2+ -> H + H", ne, te},
        {"R2_2_14_H4", "amjuel", 4, 4844, "H.4 2.2.14", "e + H2+ -> H + H", ne, te},
        {"R2_2_14_H8", "amjuel", 8, 7147, "H.8 2.2.14", "energy loss for e + H2+(v)", ne, te},
        {"R2_2_9_H4", "amjuel", 4, 4628, "H.4 2.2.9", "e + H2 -> 2e + H2+", ne, te},
        {"R2_1_5_H4", "amjuel", 4, 4156, "H.4 2.1.5", "e + H -> H+ + 2e", ne, te},
        {"R2_1_5_H10", "amjuel", 10, 7843, "H.10 2.1.5", "e + H -> H+ + 2e", ne, te},
        {"R2_1_8_H4", "amjuel", 4, 4329, "H.4 2.1.8", "H+ + e -> H(1s)", ne, te},
        {"R2_1_8_H10", "amjuel", 10, 7962, "H.10 2.1.8", "H+ + e -> H(1s)", ne, te},
        {"R2_2_5g_H4", "amjuel", 4, 4571, "H.4 2.2.5g", "e + H2 -> e + H + H", ne, te},
        {"R2_2_10_H4", "amjuel", 4, 4679, "H.4 2.2.10", "e + H2 -> 2e + H + H+", ne, te},
        {"R3_2_3_H2", "amjuel", 2, 2253, "H.2 3.2.3", "p + H2 -> H + H2+", ne, te},
        {"R3_2_3_H3", "amjuel", 3, 3841, "H.3 3.2.3", "p + H2(v) -> H + H2+", eproj, tback},
        {"R3_1_8_H3", "amjuel", 3, 3683, "H.3 3.1.8", "p + H(1s) -> H(1s) + p", eproj, tback},
        {"R0_3T_H3", "amjuel", 3, 2972, "H.3 0.3T", "p + H2 elastic", eproj, tback},
        {"R0_3D_H3", "amjuel", 3, 3026, "H.3 0.3D", "p + H2 elastic", eproj, tback},
        {"R2_2_12_H4", "amjuel", 4, 4788, "H.4 2.2.12", "e + H2+ -> e + H + H+", ne, te},
        {"R2_2_11_H4", "amjuel", 4, 4732, "H.4 2.2.11", "e + H2+ -> 2e + H+ + H+", ne, te},
        {"R2_2_17_H2", "amjuel", 2, 2188, "H.2 2.2.17", "e + H2+ -> H + H-", ne, te},
        {"R7_2_3a_H4", "amjuel", 4, 6888, "H.4 7.2.3a", "H- channel a", ne, te},
        {"R7_2_3b_H4", "amjuel", 4, 6948, "H.4 7.2.3b", "H- channel b", ne, te},
        {"R5_3_1_H3", "hydhel", 3, 5533, "H.3 5.3.1", "He+ + He -> He + He+", eproj, tback},
        {"R_H_H", "ammonx_elsa", 2, 1026, "H.2 H-H", "H + H elastic", ne, te},
        {"R_H2_H2", "ammonx_elsa", 2, 1033, "H.2 H2-H2", "H2 + H2 elastic", ne, te},
        {"R_H_H2", "ammonx_elsa", 2, 1040, "H.2 H-H2", "H + H2 elastic", ne, te},
        {"R_H2_H", "ammonx_elsa", 2, 1050, "H.2 H2-H", "H2 + H elastic", ne, te},
    };
}

void writeReadme(const std::filesystem::path &path)
{
    std::ofstream out(path);
    out << "EIRENE reaction coefficient export\n"
        << "Generated by tools/export_eirene_reactions.cpp using the same EIRENE::cal() evaluator as the NT code.\n\n"
        << "Output files: eirene_reaction_b2/<reaction>.csv\n"
        << "Units: value_m3_s is the code-side EIRENE fit coefficient after EIRENE::cal() conversion.\n"
        << "Most exported quantities are <sigma v> rate coefficients, not raw microscopic sigma(E).\n"
        << "H.8 entries are energy-loss coefficients and use the corresponding H.2 value as helper input.\n"
        << "For electron reactions, each B2 cell uses ne_2D and te_2D from the setting Casepath.\n"
        << "For H.3 heavy-particle fits, the B2 Te value is used as the local energy/temperature argument so the output keeps B2 shape.\n";
}
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

    std::filesystem::create_directories(paths.output_path);
    const std::filesystem::path output_dir(paths.output_path);
    const std::filesystem::path b2_output_dir = output_dir / "eirene_reaction_b2";
    std::filesystem::create_directories(b2_output_dir);
    const std::filesystem::path summary_path = output_dir / "eirene_reaction_b2_summary.csv";
    const std::filesystem::path readme_path = output_dir / "eirene_reaction_coefficients_readme.txt";

    std::ofstream summary(summary_path);
    if (!summary.is_open())
        throw std::runtime_error("Cannot write summary CSV: " + summary_path.string());
    summary << "name,database,fit,line,eirene_table,reaction,file,rows\n";

    for (const ReactionSpec &spec : reactionSpecs())
    {
        EIRENE reaction(spec.fit, spec.line, spec.database);
        const std::filesystem::path reaction_path =
            b2_output_dir / (sanitizeFileName(spec.name + "_" + spec.eirene_table) + ".csv");
        std::ofstream out(reaction_path);
        if (!out.is_open())
            throw std::runtime_error("Cannot write reaction CSV: " + reaction_path.string());
        out << std::setprecision(12)
            << "i,j,ne_m-3,Te_eV,argument1_name,argument1_value,argument2_name,argument2_value,value_m3_s\n";

        long long rows = 0;
        for (int i = 0; i < paths.n_poloidal; ++i)
        {
            for (int j = 0; j < paths.n_radial; ++j)
            {
                const double arg1 = spec.fit == 3 ? te[i][j] : ne[i][j];
                const double arg2 = te[i][j];
                double value = 0.0;
                if (spec.name == "R2_2_14_H8")
                {
                    EIRENE helper(2, 2157, "amjuel");
                    double h2_value = helper.cal(ne[i][j], te[i][j]);
                    value = reaction.cal(ne[i][j], te[i][j], &h2_value);
                }
                else
                {
                    value = reaction.cal(arg1, arg2);
                }
                out << i << ',' << j << ','
                    << ne[i][j] << ',' << te[i][j] << ','
                    << csvEscape(spec.argument1_name) << ',' << arg1 << ','
                    << csvEscape(spec.argument2_name) << ',' << arg2 << ','
                    << value << '\n';
                ++rows;
            }
        }
        summary << csvEscape(spec.name) << ','
                << csvEscape(spec.database) << ','
                << spec.fit << ','
                << spec.line << ','
                << csvEscape(spec.eirene_table) << ','
                << csvEscape(spec.reaction) << ','
                << csvEscape(reaction_path.string()) << ','
                << rows << '\n';
    }

    writeReadme(readme_path);
    std::cout << "Read " << ne_path << '\n';
    std::cout << "Read " << te_path << '\n';
    std::cout << "Wrote reaction files under " << b2_output_dir << '\n';
    std::cout << "Wrote " << summary_path << '\n';
    std::cout << "Wrote " << readme_path << '\n';
    return 0;
}
