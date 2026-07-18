#define _CRT_SECURE_NO_WARNINGS
#include "Read.h"
#include "Particle.h"
#include "DBoundarySource.h"

#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>

namespace
{
    std::ifstream fp;
    std::string Path_Grid,
        Path_Te, Path_Ti, Path_ne, Path_ni, Path_Vol, Path_ua, Path_NumPar_l, Path_NumPar_r,
        Path_n_D2_0, Path_n_D_0, Path_T_D_0, Path_D2_1, Path_T_D2_0,
        Path_n_T2_0, Path_n_T_0, Path_T_T_0, Path_T2_1, Path_T_T2_0,
        Path_Ei_Dion_l, Path_Ei_Dion_r, Path_bb, Path_bpol,
        Path_brad, Path_btor, Path_Epol, Path_Erad, Path_Dn;
    int a, b, c;

    void TargetSectionRead(const std::string &path_l,
                           const std::string &path_r,
                           const std::string &section,
                           double values[])
    {
        const auto read_side = [&](const std::string &path, int offset)
        {
            std::ifstream input(path);
            if (!input)
                throw std::runtime_error("Cannot open target profile: " + path);

            std::string line;
            bool found = false;
            int count = 0;
            while (std::getline(input, line))
            {
                if (!found)
                {
                    found = line.find(section) != std::string::npos;
                    continue;
                }
                if (!line.empty() && line.front() == '#')
                    break;

                std::istringstream row(line);
                double coordinate = 0.0;
                double value = 0.0;
                if (!(row >> coordinate >> value))
                    continue;
                values[offset + count] = std::abs(value);
                if (++count == N_radial)
                    return;
            }

            if (!found)
                throw std::runtime_error(
                    "Target profile section '" + section + "' is missing in " + path);
            throw std::runtime_error(
                "Target profile section '" + section + "' is incomplete in " + path);
        };

        read_side(path_l, 0);
        read_side(path_r, N_radial);
    }
}

void Read()
{
    /// Grid reading
    if (K_GRID)
    {
        Grid1.GRID_Read(N_poloidal, N_radial, Casepath);
    }
    // std::cout << Grid1.PLasma_Grid_Boundry_num() << endl;

    Path_Grid = Casepath + "shapedata/com_nx1";
    GridRead(Path_Grid, 0);

    Path_Grid = Casepath + "shapedata/com_nx2";
    GridRead(Path_Grid, 1);

    Path_Grid = Casepath + "shapedata/com_nx3";
    GridRead(Path_Grid, 2);

    Path_Grid = Casepath + "shapedata/com_nx4";
    GridRead(Path_Grid, 3);

    Path_Grid = Casepath + "shapedata/com_ny1";
    GridRead(Path_Grid, 4);

    Path_Grid = Casepath + "shapedata/com_ny2";
    GridRead(Path_Grid, 5);

    Path_Grid = Casepath + "shapedata/com_ny3";
    GridRead(Path_Grid, 6);

    Path_Grid = Casepath + "shapedata/com_ny4";
    GridRead(Path_Grid, 7);

    /// read the triangle mesh
    Grid4.load_eirene_mesh(Casepath + "solps_output/fort.33", Casepath + "solps_output/fort.34", Casepath + "solps_output/fort.35", Casepath + "shapedata/wall_segment.txt");

    /// process the triangle information
    Grid4.mesh_find(N_poloidal, N_radial);

    // Grid4.judge_orientation(); //calculate the orientation, but all of triangle mesh are Clockwise

    /// read the eirene data
    // Grid4.read_b2fgmtry(Casepath + "solps_output/b2fgmtry");
    // Grid4.read_fort44(Casepath + "solps_output/fort.44");
    // auto dab2_tri = Grid4.map_fort44_to_triangles("dab2", 3); // 3邻居IDW

    /*for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            std::cout << Grid1.Plasma_Grid(i, j).Grid_Point(0).x() - Grid[i][j][0] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(1).x() - Grid[i][j][1] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(2).x() - Grid[i][j][2] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(3).x() - Grid[i][j][3] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(0).y() - Grid[i][j][4] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(1).y() - Grid[i][j][5] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(2).y() - Grid[i][j][6] << '\t'
                      << Grid1.Plasma_Grid(i, j).Grid_Point(3).y() - Grid[i][j][7] << endl;
        }
    }*/

    /// Tokamak wall position reading(Clockwise)

    /// test Grid.Find()
    /*double x_temp[3];
    int Index_temp, XY_temp[2], Zone_temp;
    D.SetX_new(0, 2.149370057975554);
    D.SetX_new(1, -0.3580073909978209);
    D.Find();
    x_temp[0] = 2.149370057975554;
    x_temp[1] = -0.3580073909978209;
    Grid1.Find(x_temp, &Zone_temp, XY_temp);
    std::cout << "test: " << XY_temp[0] << ',' << XY_temp[1] << '\t' << Zone_temp << endl;
    std::cout << "old: " << D.XY(0) << ',' << D.XY(1) << '\t' << D.Zone() << endl;*/
    /*for (int i = 0; i < 10000; i++)
    {
        Zone_temp = 0;
        x_temp[0] = Tools::Random() * (Grid1.Rectangular_End(0) - Grid1.Rectangular_Start(0)) + Grid1.Rectangular_Start(0);
        x_temp[1] = Tools::Random() * (Grid1.Rectangular_End(1) - Grid1.Rectangular_Start(1)) + Grid1.Rectangular_Start(1);
        D.SetX_new(0, x_temp[0]);
        D.SetX_new(1, x_temp[1]);
        D.Find();
        // std::cout << x_temp[0] << ',' << x_temp[1] << endl;
        Grid1.Find(x_temp, &Zone_temp, XY_temp);
        // std::cout << D.Zone() << ',' << D.XY(0) << ',' << D.XY(1) << endl;
        // std::cout << endl;
        if (XY_temp[0] != D.XY(0) || XY_temp[1] != D.XY(1) || Zone_temp != D.Zone())
        {
            std::cout << "error of XY" << endl;
            std::cout << x_temp[0] << ',' << x_temp[1] << endl;
            std::cout << D.Zone() << ',' << D.XY(0) << ',' << D.XY(1) << endl;
            std::cout << endl;
        }
    }*/
    /*for (int i = 1; i < 38; i++)
    {
        std::cout << Grid[1][i][0] << '\t' << Grid[1][i][4] << endl;
    }
    std::cout << endl;
    for (int i = 1; i < 38; i++)
    {
        std::cout << Grid[96][38 - i][1] << '\t' << Grid[96][38 - i][5] << endl;
    }*/

    /// Read the volume
    Path_Vol = Casepath + "2D_data/vol_2D.dat";
    PlasmaRead(Path_Vol, Volume);

    ////Plasma parameter reading
    Path_Ti = Casepath + "2D_data/ti_2D.dat";
    PlasmaRead(Path_Ti, Ti);
    Path_ne = Casepath + "2D_data/ne_2D.dat";
    PlasmaRead(Path_ne, ne);
    Path_Te = Casepath + "2D_data/te_2D.dat";
    PlasmaRead(Path_Te, Te);
    Path_ua = Casepath + "2D_data/ua_Dion_2D.dat";
    PlasmaRead(Path_ua, ua_D_1);

    Path_ni = Casepath + "2D_data/nDion_2D.data";
    PlasmaRead(Path_ni, n_D_1);
    Path_T_D_0 = Casepath + "2D_data/tDatom_2D.data";
    PlasmaRead(Path_T_D_0, T_D_0, T_D_0_Tri);
    Path_T_D2_0 = Casepath + "2D_data/tDmolecule_2D.data";
    PlasmaRead(Path_T_D2_0, T_D2_0, T_D2_0_Tri);
    Path_n_D2_0 = Casepath + "2D_data/nDmolecule_2D.data";
    PlasmaRead(Path_n_D2_0, n_D2_0, n_D2_0_Tri);
    Path_n_D_0 = Casepath + "2D_data/nDatom_2D.data";
    PlasmaRead(Path_n_D_0, n_D_0, n_D_0_Tri);
    Path_D2_1 = Casepath + "2D_data/nDmoleculeion_2D.data";
    PlasmaRead(Path_D2_1, n_D2_1);

    if (K_T)
    {
        Path_ni = Casepath + "2D_data/nTion_2D.data";
        PlasmaRead(Path_ni, n_T_1);
        Path_T_T_0 = Casepath + "2D_data/tTatom_2D.data";
        PlasmaRead(Path_T_T_0, T_T_0, T_T_0_Tri);
        Path_T_T2_0 = Casepath + "2D_data/tTmolecule_2D.data";
        PlasmaRead(Path_T_T2_0, T_T2_0, T_T2_0_Tri);
        Path_n_T2_0 = Casepath + "2D_data/nTmolecule_2D.data";
        PlasmaRead(Path_n_T2_0, n_T2_0, n_T2_0_Tri);
        Path_n_T_0 = Casepath + "2D_data/nTatom_2D.data";
        PlasmaRead(Path_n_T_0, n_T_0, n_T_0_Tri);
        Path_T2_1 = Casepath + "2D_data/nTmoleculeion_2D.data";
        PlasmaRead(Path_T2_1, n_T2_1);
        Path_ua = Casepath + "2D_data/ua_Tion_2D.dat";
        PlasmaRead(Path_ua, ua_T_1);
    }
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        if (!Grid4.if_in_plasmagrid(i))
        {
            T_D_0_Tri[i] = 0.25;
            T_D2_0_Tri[i] = 0.25;
            n_D_0_Tri[i] = 1e10;
            n_D2_0_Tri[i] = 1e10;
            if (K_T)
            {
                T_T_0_Tri[i] = 0.25;
                T_T2_0_Tri[i] = 0.25;
                n_T_0_Tri[i] = 1e10;
                n_T2_0_Tri[i] = 1e10;
            }
        }
    }

    Path_bb = Casepath + "2D_data/bb_2D.dat";
    PlasmaRead(Path_bb, Btot);
    Path_bpol = Casepath + "2D_data/bpol_2D.dat";
    PlasmaRead(Path_bpol, Bpol);
    Path_btor = Casepath + "2D_data/btor_2D.dat";
    PlasmaRead(Path_btor, Btor);
    Path_brad = Casepath + "2D_data/brad_2D.dat";
    PlasmaRead(Path_brad, Brad);
    Path_Epol = Casepath + "2D_data/Ep_2D.dat";
    PlasmaRead(Path_Epol, Epol);
    Path_Erad = Casepath + "2D_data/Er_2D.dat";
    PlasmaRead(Path_Erad, Erad);

    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            Etor[i][j] = 0.;
        }
    }

    T_0_read_Tri(Casepath + "D_0_n_Tri", n_D_0_Tri);
    T_0_read_Tri(Casepath + "D2_0_n_Tri", n_D2_0_Tri);
    T_0_read_Tri(Casepath + "D_0_T_Tri", T_D_0_Tri);
    T_0_read_Tri(Casepath + "D2_0_T_Tri", T_D2_0_Tri);
    V_0_read_Tri(Casepath + "D_0_V_Tri", ua_D_0_Tri);
    V_0_read_Tri(Casepath + "D2_0_V_Tri", ua_D2_0_Tri);
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        // Preserve the exported neutral temperature for collision-partner sampling.
        // Reaction fits apply their own validity-range floor when evaluated.
        if (T_D_0_Tri[i] < 0.)
            T_D_0_Tri[i] = 0.;
        if (T_D2_0_Tri[i] < 0.)
            T_D2_0_Tri[i] = 0.;
    }

    // Read the gamma from target recycle
    if (K_D)
    {
        Path_NumPar_l = Casepath + "profiles_data/recycled_neutral_flux_D_l.data";
        Path_NumPar_r = Casepath + "profiles_data/recycled_neutral_flux_D_r.data";
        TargetRead(Path_NumPar_l, Path_NumPar_r, NumPar_wall_D);
    }
    if (K_T)
    {
        Path_NumPar_l = Casepath + "profiles_data/recycled_neutral_flux_T_l.data";
        Path_NumPar_r = Casepath + "profiles_data/recycled_neutral_flux_T_r.data";
        TargetRead(Path_NumPar_l, Path_NumPar_r, NumPar_wall_T);
    }

    if (K_Ei == 1)
    {
        Path_Ei_Dion_l = Casepath + "2D_data/ei_Dion_l.data";
        Path_Ei_Dion_r = Casepath + "2D_data/ei_Dion_r.data";
        TargetSectionRead(
            Path_Ei_Dion_l, Path_Ei_Dion_r, "D@S@+1@N@", Ei_Dion);
    }

    // Dn read
    if (K_dn == 2)
    {
        Path_Dn = Casepath + "dn3da.data";
        radial_Read(Path_Dn, Dn);
    }
    if (K_Reflect == 1)
    {
        H_W.SetReflect(H_W_n_a, H_W_E_a, H_W_n_b, H_W_E_b, H_W_epsilon_L);
        D_W.SetReflect(D_W_n_a, D_W_E_a, D_W_n_b, D_W_E_b, D_W_epsilon_L);
        T_W.SetReflect(T_W_n_a, T_W_E_a, T_W_n_b, T_W_E_b, T_W_epsilon_L);
        H_C.SetReflect(H_C_n_a, H_C_E_a, H_C_n_b, H_C_E_b, H_C_epsilon_L);
        D_C.SetReflect(D_C_n_a, D_C_E_a, D_C_n_b, D_C_E_b, D_C_epsilon_L);
        T_C.SetReflect(T_C_n_a, T_C_E_a, T_C_n_b, T_C_E_b, T_C_epsilon_L);
        H_Be.SetReflect(H_Be_n_a, H_Be_E_a, H_Be_n_b, H_Be_E_b, H_Be_epsilon_L);
        D_Be.SetReflect(D_Be_n_a, D_Be_E_a, D_Be_n_b, D_Be_E_b, D_Be_epsilon_L);
        T_Be.SetReflect(T_Be_n_a, T_Be_E_a, T_Be_n_b, T_Be_E_b, T_Be_epsilon_L);
    }
    if (K_Reflect == 2) // 2 for Trim database;
    {
        H_W.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_H_on_W_1.00", Inputpath + "database/RE_H_on_W_1.00");
        D_W.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_D_on_W_1.00", Inputpath + "database/RE_D_on_W_1.00");
        T_W.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_T_on_W_1.00", Inputpath + "database/RE_T_on_W_1.00");
        H_C.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_H_on_C_1.00", Inputpath + "database/RE_H_on_C_1.00");
        D_C.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_D_on_C_1.00", Inputpath + "database/RE_D_on_C_1.00");
        T_C.ReadTrimData(K_Reflect, 12, 7, Inputpath + "database/Rn_T_on_C_1.00", Inputpath + "database/RE_T_on_C_1.00");
    }
    if (K_Reflect == 3) // 3 for calculate from Trim
    {
        H_W.ReadTrimData(K_Reflect, 25, 11, Inputpath + "database/hw_rn.dat", Inputpath + "database/hw_re.dat");
        D_W.ReadTrimData(K_Reflect, 25, 11, Inputpath + "database/dw_rn_1.dat", Inputpath + "database/dw_re_1.dat");
        T_W.ReadTrimData(K_Reflect, 25, 11, Inputpath + "database/tw_rn.dat", Inputpath + "database/tw_re.dat");
        std::cout << "the C wall is lack" << endl;
    }
    if (K_DWTrimReflection == 1)
    {
        if (!D_W_Trim.Load(Inputpath + "database/" + DWTrimDatabase))
            throw std::runtime_error("Failed to load D-on-W reflection database: " + DWTrimDatabase);
        std::cout << "D-on-W TRIM distribution model enabled: " << DWTrimDatabase << endl;
    }
    if (K_DBoundarySource)
    {
        if (!K_D)
            throw std::runtime_error("K_DBoundarySource requires K_D=1");
        std::filesystem::path flux_path(DBoundaryFluxFile);
        if (flux_path.is_relative())
            flux_path = std::filesystem::path(Casepath) / flux_path;
        D_BoundarySource.Load(flux_path.string());
    }
    /*for (double i = 1; i < 1000; i++)
    {
        std::cout << i << '\t' << D_W.n_RefCoeff(K_Reflect, i) << '\t' << T_W.n_RefCoeff(K_Reflect, i) << "\t";
        std::cout << D_W.n_RefCoeff(K_Reflect, i, 45) << '\t' << T_W.n_RefCoeff(K_Reflect, i, 45) << "\t";
        std::cout << D_W.n_RefCoeff(K_Reflect, i, 60) << '\t' << T_W.n_RefCoeff(K_Reflect, i, 60) << endl;
    }*/
    /*for (int k = 0; k < Grid4.num_tris(); k++)
    {

        /// b2 Grid Addressing
        double B2_[2] = {1. / 3 * (Grid4.nodes_[Grid4.tris_[k].v[0]].r + Grid4.nodes_[Grid4.tris_[k].v[1]].r + Grid4.nodes_[Grid4.tris_[k].v[2]].r),
                         1. / 3 * (Grid4.nodes_[Grid4.tris_[k].v[0]].z + Grid4.nodes_[Grid4.tris_[k].v[1]].z + Grid4.nodes_[Grid4.tris_[k].v[2]].z)};

        double secx, secy;
        bool found = false;
        Grid4.SetB2Index_find(k, -1, -1);
        for (int m = 0; m < N_poloidal; m++)
        {
            for (int n = 0; n < N_radial; n++)
            {
                if (Grid1.Plasma_Grid(m, n).Ifingrid(B2_[0], B2_[1]))
                {
                    Grid4.SetB2Index_find(k, m, n);
                    found = true;
                    break;
                }
            }
            if (found == true)
                break;
        }
    }*/
    /*ofstream Out_temp("doc/B2Index_fixed.txt");
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        Out_temp << i << "\t" << Grid4.b2_index(i, 0) << "\t" << Grid4.b2_index(i, 1) << endl;
    }
    Out_temp.close();*/
}

void GridRead(string Path, int n)
{
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for grid have some problem!!!\n";
    }
    for (int i = 0; i < N_poloidal; i++)
        for (int j = 0; j < N_radial; j++)
        {
            fp >> Grid[i][j][n];
        }
    fp.close();
}

void radial_Read(string Path, double Var[])
{
    string line;
    std::ifstream fp;
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        std::cerr << "This file READING for radialRead have some problem!!!\n";
    }
    getline(fp, line);
    double temp[N_RADIAL_GRID];
    for (int i = 0; i < N_radial; i++)
    {
        fp >> temp[i] >> Var[i];
        if (Var[i] < 0)
            Var[i] *= -1;
    }
    fp.close();
}

void TargetRead(string Path_l, string Path_r, double Var[])
{
    string line;
    std::ifstream fp;
    fp.open(Path_l, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        std::cerr << "This file READING for " + Path_l + " have some problem!!!\n";
    }
    getline(fp, line);
    double temp[N_RADIAL_GRID];
    for (int i = 0; i < N_radial; i++)
    {
        fp >> temp[i] >> Var[i];
        if (Var[i] < 0)
            Var[i] *= -1;
    }
    fp.close();
    fp.open(Path_r, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        std::cerr << "This file READING for " + Path_r + " have some problem!!!\n";
    }
    getline(fp, line);
    for (int i = 0; i < N_radial; i++)
    {
        fp >> temp[i] >> Var[i + N_radial];
        if (Var[i + N_radial] < 0)
            Var[i + N_radial] *= -1;
    }
    fp.close();
}

void PlasmaRead(string Path, std::vector<std::vector<double>> &Plasma, std::vector<double> &Plasma_Tri)
{
    string line;
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    getline(fp, line);
    for (int i = 0; i < gridCellCount(); i++)
    {
        fp >> a >> b >> c;
        fp >> Plasma[a + 1][b + 1];
    }
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        if (Grid4.if_in_plasmagrid(i))
        {
            Plasma_Tri[i] = Plasma[Grid4.tris_[i].neigh[9]][Grid4.tris_[i].neigh[10]];
        }
    }
    /*ofstream out_temp("doc/Plasma.txt");
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out_temp << i << "\t" << Grid4.tris_[i].neigh[9] << "\t";
        out_temp << Grid4.tris_[i].neigh[10] << "\t";
        if (Grid4.if_in_plasmagrid(i))
            out_temp << Plasma[Grid4.tris_[i].neigh[9]][Grid4.tris_[i].neigh[10]] << "\t";
        else
            out_temp << 0. << "\t";
        out_temp << Plasma_Tri[i] << endl;
    }
    out_temp.close();
    fp.close();*/
}

void V_0_read_Tri(string Path, std::vector<std::vector<double>> &Plasma)
{
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        fp >> Plasma[i][0] >> Plasma[i][1] >> Plasma[i][2];
    }
}

void T_0_read_Tri(string Path, std::vector<double> &Plasma)
{
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        fp >> Plasma[i];
    }
}

void PlasmaRead(string Path, std::vector<std::vector<double>> &Plasma)
{
    string line;
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    getline(fp, line);
    for (int i = 0; i < gridCellCount(); i++)
    {
        fp >> a >> b >> c;
        fp >> Plasma[a + 1][b + 1];
    }
}

void PlasmaRead(string Path, std::vector<double> &Plasma)
{
    string line;
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    double Plasma_temp[N_POLOIDAL_GRID][N_RADIAL_GRID];
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    getline(fp, line);
    for (int i = 0; i < gridCellCount(); i++)
    {
        fp >> a >> b >> c;
        fp >> Plasma_temp[a + 1][b + 1];
    }
    fp.close();
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        if (Grid4.if_in_plasmagrid(i))
        {
            Plasma[i] = Plasma_temp[Grid4.tris_[i].neigh[9]][Grid4.tris_[i].neigh[10]];
        }
    }
    /*ofstream out_temp("doc/Plasma.txt");
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out_temp << i << "\t" << Plasma[i] << endl;
    }
    out_temp.close();*/
}

void MatrixRead(int M, int N, string Path, double **Matrix)
{
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            fp >> Matrix[i][j];
    }
    fp.close();
}
