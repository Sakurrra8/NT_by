#include "Global.h"

void SOLPSPlasmaRead(string Path, std::vector<std::vector<double>> &Plasma);
void OutPec(string Path, string name, ADAS *PEC, int wave, std::vector<std::vector<double>> &n_i);
void MatrixRead(string Path, double **Matrix);

std::vector<std::vector<double>> n_e, T_e, n_Li_0, n_Li_1, n_Li_2, n_Li_3, n_D_0_pec, n_D2, n_W_1;
int N_poloidal_pec = 136;
int N_radial_pec = 40;
void CalPec()
{
    /// ADAS pec read
    PEC_Li1.readname("pec96#li_pjr#li1.dat");
    PEC_Li1.readpec(1, 1);
    PEC_Li2.readname("pec96#li_pju#li2.dat");
    PEC_Li2.readpec(1, 1);

    PEC_W_1.readname("pec40#w_ca#w1.dat");
    PEC_W_1.readpec(19, 0);
    PEC_W_8.readname("pec40#w_ls#w8.dat");
    PEC_W_8.readpec(19, 0);
    PEC_W_11.readname("pec40#w_ls#w11.dat");
    PEC_W_11.readpec(19, 0);
    PEC_W_15.readname("pec40#w_ic#w15.dat");
    PEC_W_15.readpec(19, 0);
    PEC_W_22.readname("pec40#w_ic#w22.dat");
    PEC_W_22.readpec(19, 0);
    PEC_W_24.readname("pec40#w_ic#w24.dat");
    PEC_W_24.readpec(19, 0);
    PEC_B_1.readname("pec93#b_llu#b1.dat");
    PEC_B_1.readpec(1, 1);

    // std::cout << PEC_W_1.PecCal(6.89E-02, 1e19, 0) << endl;
    n_e.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    T_e.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_Li_0.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_Li_1.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_Li_2.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_Li_3.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_D_0_pec.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0)); // defined in Global.cpp at first
    // n_D_1.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));
    n_W_1.assign(N_poloidal_pec, std::vector<double>(N_radial_pec, 0.0));

    /// read the SOLPS data
    string PecPath;
    PecPath = "/home/bianyu/nt/inputfile_Pec/case3/";
    /*SOLPSPlasmaRead(PecPath + "ne_2D.data", n_e);
    SOLPSPlasmaRead(PecPath + "te_2D.data", T_e);
    SOLPSPlasmaRead(PecPath + "nDatom_2D.data", n_D_0);
    // SOLPSPlasmaRead(PecPath + "nD+_2D.data", n_D_1);
    SOLPSPlasmaRead(PecPath + "nLi0_2D.data", n_Li_0);
    SOLPSPlasmaRead(PecPath + "nLi1_2D.data", n_Li_1);
    SOLPSPlasmaRead(PecPath + "nLi2_2D.data", n_Li_2);
    SOLPSPlasmaRead(PecPath + "nLi3_2D.data", n_Li_3);*/
    // MatrixRead(PecPath + "steady_distri_n_q1.txt", n_W_1);

    /// output the emission intensity
    string Path_pecout = "/data/leuven/379/vsc37950/nt/doc/pec/";
    ofstream pecout;
    pecout.open(Path_pecout + "Li1_14.txt");
    for (double i = 0.1; i < 1000; i++)
    {
        /*pecout << i << "\t" << PEC_B_1.PecCal(i, 1.2e19, 13) << "\t" << PEC_Li1.PecCal(i, 1.2e19, 5) << "\t";
        pecout << PEC_Li1.PecCal(i, 1.2e19, 67) << "\t";
        pecout << SCD12_H.cal(i, 1.2e19) << "\t" << SCD89_B.cal(i, 1.2e19, 1) << "\t" << ACD89_B.cal(i, 1.2e19, 2) << endl;*/
        pecout << i << "\t" << SCD96_Li.cal(i, 1.2e19, 1) << "\t" << ACD96_Li.cal(i, 1.2e19, 2) << endl;
    }
    /*for (int i = 0; i < PEC_Li1.num_waveLength(); i++)
    {
        if (i < 62)
            OutPec(Path_pecout, "pec_Li_1_" + std::to_string(i + 1), &PEC_Li1, i, n_Li_1);
        if (i >= 62)
            OutPec(Path_pecout, "pec_Li_1_" + std::to_string(i + 1), &PEC_Li1, i, n_Li_2);
    }
    for (int i = 0; i < PEC_Li2.num_waveLength(); i++)
    {
        if (i < (PEC_Li2.num_waveLength() + 1) / 2)
            OutPec(Path_pecout, "pec_Li_2_" + std::to_string(i + 1), &PEC_Li2, i, n_Li_2);
        if (i >= (PEC_Li2.num_waveLength() + 1) / 2)
            OutPec(Path_pecout, "pec_Li_2_" + std::to_string(i + 1), &PEC_Li2, i, n_Li_3);
    }*/
    /*for (int i = 0; i < N_poloidal_pec; i++)
    {
        for (int j = 0; j < N_radial_pec; j++)
        {
            std::cout << n_e[i][j] * n_Li_1[i][j] * PEC_Li1.PecCal(T_e[i][j], n_e[i][j], 5) + n_e[i][j] * n_Li_1[i][j] * PEC_Li1.PecCal(T_e[i][j], n_e[i][j], 36) + n_e[i][j] * n_Li_2[i][j] * PEC_Li1.PecCal(T_e[i][j], n_e[i][j], 67) << '\t';
        }
        std::cout << endl;
    }*/
    // OutPec(Path_pecout, "pec_W_1_" + std::to_string(11), &PEC_W_1, 10, n_W_1);
    /*ofstream pecout;
    pecout.open(Path_pecout + "cross_540_WLi.txt");
    for (double i = 0.1; i < 1000; i += 0.1)
    {
        pecout << i << '\t';
        pecout << PEC_Li1.PecCal(i, 1.2e19, 5) << '\t';
        pecout << PEC_Li1.PecCal(i, 1.2e19, 36) << '\t';
        pecout << PEC_Li1.PecCal(i, 1.2e19, 67) << '\t';

        pecout << PEC_W_1.PecCal(i, 1.2e19, 10) << '\t';
        pecout << PEC_W_8.PecCal(i, 1.2e19, 48) << '\t';
        pecout << PEC_W_11.PecCal(i, 1.2e19, 46) << '\t';
        pecout << PEC_W_15.PecCal(i, 1.2e19, 49) << '\t';
        pecout << PEC_W_22.PecCal(i, 1.2e19, 47) << '\t';
        pecout << PEC_W_22.PecCal(i, 1.2e19, 48) << '\t';
        pecout << PEC_W_24.PecCal(i, 1.2e19, 39) << '\t';

        for (int j = 0; j < PEC_Li1.num_waveLength(); j++)
        {
            pecout << PEC_Li1.PecCal(i, 1.2e19, j) << '\t';
        }
        pecout << endl;
    }*/
}

void SOLPSPlasmaRead(string Path, std::vector<std::vector<double>> &Plasma)
{
    int a, b, c;
    std::ifstream stream; // file stream
    string line;
    stream.open(Path, std::ios::in); // ios::in means read
    if (!stream.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    std::getline(stream, line);
    for (int i = 0; i < N_poloidal_pec * N_radial_pec; i++)
    {
        stream >> a >> b >> c;
        stream >> Plasma[a + 1][b + 1];
    }
    stream.close();
}

void OutPec(string Path, string name, ADAS *PEC, int wave, std::vector<std::vector<double>> &n_i)
{
    ofstream pecout;
    pecout.open(Path + name);
    for (int i = 0; i < N_poloidal_pec; i++)
    {
        for (int j = 0; j < N_radial_pec; j++)
        {
            pecout << n_e[i][j] * n_i[i][j] * PEC->PecCal(T_e[i][j], n_e[i][j], wave) << '\t';
        }
        pecout << endl;
    }
}

void MatrixRead(string Path, double **Matrix)
{
    std::ifstream fp;            // file stream
    fp.open(Path, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path + " have some problem!!!\n";
    }
    for (int i = 0; i < N_poloidal_pec; i++)
    {
        for (int j = 0; j < N_radial_pec; j++)
            fp >> Matrix[i][j];
    }
    fp.close();
}