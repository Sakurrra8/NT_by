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
void V_0_read_Tri(string Path, std::vector<std::vector<double>> &Plasma);
void T_0_read_Tri(string Path, std::vector<double> &Plasma);

