#pragma once
#include "Global.h"
#define _CRT_SECURE_NO_WARNINGS

using namespace std;

void PlasmaRead(string Path, std::vector<std::vector<double>> &Plasma, std::vector<double> &Plasma_Tri);
void PlasmaRead(string Path, std::vector<std::vector<double>> &Plasma);
void PlasmaRead(string Path, std::vector<double> &Plasma);
ifstream &seekline(ifstream &in, int line);
void GridRead(string Path, int n);
void TargetRead(string Path_l, string Path_r, double Var[]);
void radial_Read(string Path, double Var[]);
void MatrixRead(int M, int N, string Path, double **Matrix);

std::ifstream fp; // file stream;
string Path_Grid,
    Path_Te, Path_Ti, Path_ne, Path_ni, Path_Vol, Path_ua, Path_NumPar_l, Path_NumPar_r,
    Path_n_D2_0, Path_n_D_0, Path_T_D_0, Path_D2_1, Path_T_D2_0,
    Path_n_T2_0, Path_n_T_0, Path_T_T_0, Path_T2_1, Path_T_T2_0,
    Path_Ei_Dion_l, Path_Ei_Dion_r, Path_bb, Path_bpol,
    Path_brad, Path_btor, Path_Epol, Path_Erad, Path_Dn;
int a, b, c;