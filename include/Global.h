#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
// #include <string.h>
#include <string>
#include <cassert>
#include <malloc.h>
#include <unistd.h>
#include <iomanip>
#include <math.h>
#include <array>
#include "float.h"
#include "Reflect.h"
#include "EIRENE.h"
#include "ADAS.h"
#include "Grid.h"
#include "WallImpactTracker.h"
#include "WallEro.h"
#include "utils.h"
// #include <sys/stat.h>
// #include <sys/types.h>
// #include <dirent.h>
// #include <stdio.h>
// #include <stdint.h>
#include "/data/leuven/379/vsc37950/nt/opt/hdf5/include/H5Cpp.h"
#include "/data/leuven/379/vsc37950/nt/opt/hdf5/include/hdf5.h"
// #include "H5Cpp.h"
// #include <stdlib.h>

#define Hmass 1.674e-27
#define Dmass 3.34e-27
#define Tmass 5.010e-27
#define H2mass 3.348e-27
#define D2mass 6.695e-27
#define T2mass 1.002e-26
#define CD4mass 2.665e-26
#define CD3mass 2.498e-26
#define CD2mass 2.330e-26
#define CD1mass 2.163e-26
#define Cmass 1.99518e-26
#define Armass 6.6359e-26
#define pi 3.14159265
#define qe 1.602e-19
#define e 2.71828

using namespace std;
using namespace H5;
using namespace eirene;

void Initialize(int Input, char *settingfile[]);
void Read();
void CalPec();
void Prepare();
void Moncar();
void Output();
void Dump_2D_Global(double Var[98][38], string name);
void flight(int i, int j);
double Xrecycling(double rand_, int a, int b);
void read_extend_plasma();

////Plasma and device parameter
extern int N_poloidal, N_radial;
extern bool K_log, StepLog, K_H2_elastic, K_EcrossBDrift, backGridBoundry, K_H, K_D, K_T;
extern int K_CX_impurity, K_C, K_ReflectDirection, K_Reflect, K_Prob, K_back, K_dn, K_CX_DT, K_DT, K_flight, K_GRID, K_Pump, K_test1, K_test2, K_test3, K_Tn, K_Vi, K_Wallelement;
extern int numPar_flight, numPar_flight_Target, IfOut, K_Recyc, K_Rec, K_Maxwell, K_Ar, K_abnormal_transport, K_D2Flight, *K_Puff, K_database_Pra, K_mu, K_PartoPar;
extern int K_database, Num_Reflect, Num_Reflect_Core, numPar_flight_CD4, K_Ei, K_WallRefl, K_backScatter, K_MarColl, K_Methane, numPar_flight_D2, numPar_flight_T2;
extern int K_T_array, K_NNCs;
extern double dt, Te_core, ne_core, NumPar_now, Tn_core, T_N, T_wall, Num_D2_pump, Num_T2_pump, coefficient_D, coefficient_T;
extern double Num_CD4_pump, t_max, Ratio_T, Ratio_D_Coll, Ratio_T_Coll, Ratio_DT_Coll, coeff_puff;
extern double DTmass, D2Tmass, DT2mass, HH2mass, DD2mass, TT2mass, H2H2mass, D2D2mass, T2T2mass;
extern double Dn[38];

extern int MeshMode;

// 三维网格
extern std::vector<std::vector<std::vector<double>>> Grid;
extern std::vector<std::vector<std::vector<double>>> B;

// 真空区域参数
extern double Te_vacuum, ne_vacuum;

// 二维场量
extern std::vector<std::vector<double>> ne, Te, Ti, Te_coll, Ti_coll;

// 1 原子态密度
extern std::vector<std::vector<double>> n_H_1, n_D_1, n_T_1;

// H 相关
extern std::vector<double> n_H_0_Tri, T_H_0_Tri, T_H2_0_Tri, n_H2_0_Tri;
extern std::vector<std::vector<double>> ua_H_0_Tri, ua_H2_0_Tri;
extern std::vector<std::vector<double>> n_H_0, T_H_0, T_H2_0, n_H2_0;
extern std::vector<std::vector<double>> n_H2_1, ua_H_1;
extern std::array<double, 3> V_H_1_now, V_H_0_now, V_H2_0_now; // 原来是固定长度 3 的数组
extern std::vector<std::vector<double>> Vi_H, Ti_H_thermal;

// D 相关
extern std::vector<double> n_D_0_Tri, T_D_0_Tri, T_D2_0_Tri, n_D2_0_Tri;
extern std::vector<std::vector<double>> ua_D_0_Tri, ua_D2_0_Tri;
extern std::vector<std::vector<double>> n_D_0, T_D_0, T_D2_0, n_D2_0;
extern std::vector<std::vector<double>> n_D2_1, ua_D_1;
extern std::array<double, 3> V_D_1_now, V_D_0_now, V_D2_0_now; // 原来是固定长度 3 的数组
extern std::vector<std::vector<double>> Vi_D, Ti_D_thermal;

// T 相关
extern std::vector<double> n_T_0_Tri, T_T_0_Tri, T_T2_0_Tri, n_T2_0_Tri;
extern std::vector<std::vector<double>> ua_T_0_Tri, ua_T2_0_Tri;
extern std::vector<std::vector<double>> n_T_0, T_T_0, T_T2_0, n_T2_0;
extern std::vector<std::vector<double>> n_T2_1, ua_T_1;
extern std::array<double, 3> V_T_1_now, V_T_0_now, V_T2_0_now;
extern std::vector<std::vector<double>> Vi_T, Ti_T_thermal;

// 几何参数
extern std::vector<std::vector<double>> Volume, Area;

// 磁场
extern std::vector<std::vector<double>> Btot, Btor, Brad, Bpol;

// 电场与漂移速度
extern std::vector<std::vector<double>> Erad, Epol, Etor, vdEp, vdEr;
extern std::vector<std::vector<double>> epx, epy, erx, ery;

// 其它物理量分量
extern std::vector<std::vector<double>> crx, cry;
extern std::vector<std::vector<double>> dTip, dTir, dTix, dTiy, dTiz;

extern std::vector<double> V_1, V_2;

// extern double step_now;
// extern double cosang[410], sinang[410];
extern int num_CoreBoundry, num_GridBoundry, num_PFRBoundry;
extern double CoreBoundry[55][2], GridBoundry[100][2], PFRBoundry[55][2]; //, cosCore[60], sinCore[60];
extern std::vector<double> NumPar_H_recyc, Tn_H_recyc, NumPar_H2_recyc;
extern std::vector<double> NumPar_D_recyc, Tn_D_recyc, NumPar_D2_recyc;
extern std::vector<double> NumPar_T_recyc, Tn_T_recyc, NumPar_T2_recyc;
extern double Ei_Dion[80], InterscePoint[10][6], coeff_recyc_target, coeff_ercyc_wall;
// extern vector<double> S_target, S_wall;
extern double NumPar_wall_H[80], NumPar_wall_D[80], NumPar_wall_T[80];
//  extern double Xrecyc[80][2]; //, theta[80][2];
extern string name_Xlog, name_Vlog, Outputpath, Inputpath, Casepath;

extern std::ofstream Xlog;
extern std::ofstream Vlog;

/*extend_mesh*/
/*extern double dnx;
extern double dny;
extern int Nbc;
extern int sp_Ntotal;*/

extern vector<vector<double>> fnax_e;
extern vector<vector<double>> ne_e;
extern vector<vector<double>> ni_e;
extern vector<vector<double>> Te_e;
extern vector<vector<double>> Ti_e;
extern vector<vector<double>> dTip_e;
extern vector<vector<double>> dTir_e;
extern vector<vector<double>> ua_e;
extern vector<vector<double>> duap_e;
extern vector<vector<double>> Ex_e;
extern vector<vector<double>> Ey_e;
extern vector<vector<double>> Ez_e;
extern vector<vector<double>> Bp_e;
extern vector<vector<double>> Br_e;
extern vector<vector<double>> Bz_e;
extern vector<vector<double>> BB_e;
extern vector<vector<double>> dBp_e;
extern vector<vector<double>> eBx_e;
extern vector<vector<double>> eBy_e;
extern vector<vector<double>> eBz_e;

/// ADAS Reading preparation
extern ADAS SCD12_H; // Effective Ionisation Coefficients
extern ADAS ACD12_H; // Effective Recombination Coefficients
extern ADAS CCD96_D; // D-D Effective CX Coefficients
extern ADAS PLT12_H; // Line Emission from excitation
extern ADAS PRB12_H; // Recombination and Bremsstrahlung
extern ADAS PRC96_D; // Charge exchange
extern ADAS SCD96_T;
extern ADAS ACD96_T;
extern ADAS CCD96_T;
extern ADAS CCD96_H;
extern ADAS PLT96_T;
extern ADAS PRB96_T;
extern ADAS PRC96_T;
extern ADAS PRC96_H;
extern ADAS PLS96_T;
extern ADAS PLS89_H;

extern ADAS SCD96_He; // Effective Ionisation Coefficients
extern ADAS SCD96_He; // Effective Ionisation Coefficients
extern ADAS ACD96_He; // Effective Recombination Coefficients
extern ADAS CCD96_He; // D-D Effective CX Coefficients
extern ADAS PLT96_He; // Line Emission from excitation
extern ADAS PRB96_He; // Recombination and Bremsstrahlung Radiation
extern ADAS PRC96_He; // D-D Effective CX Radiation

extern ADAS SCD96_C;  // Effective Ionisation Coefficients
extern ADAS ACD96_C;  // Effective Recombination Coefficients
extern ADAS CCD96_C;  // D-D Effective CX Coefficients
extern ADAS PLT96_C;  // Line Emission from excitation
extern ADAS PRB96_C;  // Recombination and Bremsstrahlung Radiation
extern ADAS PRC96_C;  // D-D Effective CX Radiation
extern ADAS SCD89_Ar; // Effective Ionisation Coefficients
extern ADAS ACD89_Ar; // Effective Recombination Coefficients
extern ADAS CCD89_Ar; // D-D Effective CX Coefficients
extern ADAS PLT89_Ar; // Line Emission from excitation
extern ADAS PRB89_Ar; // Recombination and Bremsstrahlung Radiation
extern ADAS PRC89_Ar; // D-D Effective CX Radiation

extern ADAS SCD89_B;
extern ADAS ACD89_B;
extern ADAS SCD96_Li;
extern ADAS ACD96_Li;

extern ADAS PEC_Li1;
extern ADAS PEC_Li2;
extern ADAS PEC_B_1;

extern ADAS PEC_W_1;
extern ADAS PEC_W_8;
extern ADAS PEC_W_11;
extern ADAS PEC_W_15;
extern ADAS PEC_W_22;
extern ADAS PEC_W_24;

extern double H_W_n_a[4], H_W_E_a[4], H_W_n_b[6], H_W_E_b[6], H_W_epsilon_L;
extern double D_W_n_a[4], D_W_E_a[4], D_W_n_b[6], D_W_E_b[6], D_W_epsilon_L;
extern double T_W_n_a[4], T_W_E_a[4], T_W_n_b[6], T_W_E_b[6], T_W_epsilon_L;
extern double H_C_n_a[4], H_C_E_a[4], H_C_n_b[6], H_C_E_b[6], H_C_epsilon_L;
extern double D_C_n_a[4], D_C_E_a[4], D_C_n_b[6], D_C_E_b[6], D_C_epsilon_L;
extern double T_C_n_a[4], T_C_E_a[4], T_C_n_b[6], T_C_E_b[6], T_C_epsilon_L;
extern double H_Be_n_a[4], H_Be_E_a[4], H_Be_n_b[6], H_Be_E_b[6], H_Be_epsilon_L;
extern double D_Be_n_a[4], D_Be_E_a[4], D_Be_n_b[6], D_Be_E_b[6], D_Be_epsilon_L;
extern double T_Be_n_a[4], T_Be_E_a[4], T_Be_n_b[6], T_Be_E_b[6], T_Be_epsilon_L;
extern Reflect H_W;
extern Reflect D_W;
extern Reflect T_W;
extern Reflect H_C;
extern Reflect D_C;
extern Reflect T_C;
extern Reflect H_Be;
extern Reflect D_Be;
extern Reflect T_Be;
extern Reflect Ar_W;
extern Reflect He_W;

extern ofstream out;

enum EireneFit
{
    EIRENE_H0,
    EIRENE_H1,
    EIRENE_H2,
    EIRENE_H3,
    EIRENE_H4,
    EIRENE_H5,
    EIRENE_H6,
    EIRENE_H7,
    EIRENE_H8,
    EIRENE_H9,
    EIRENE_H10
};
/// EIRENE
extern EIRENE R3_2_3r_H4;  // MAR via H2+
extern EIRENE R3_2_3d_H4;  // MAD via H2+
extern EIRENE R3_2_3i_H4;  // MAI via H2+
extern EIRENE R2_2_17r_H4; // MAR via H-
extern EIRENE R2_2_17d_H4; // MAD via H-
extern EIRENE R2_2_5_H4;   //$ e + H_2 = e + H + H$
extern EIRENE R2_2_14_H2;
extern EIRENE R2_2_14_H4; //$ e + H_2^+ = H + H  $
extern EIRENE R2_2_14_H8; //$ e + H_2^+(v) = H(1s) + H^*(n)   (v=0-9, n \ge 2) $
extern EIRENE R2_2_9_H4;  //$ e + H_2 = 2e + H_2^+  $
extern EIRENE R2_1_5_H4;  //$ H + e = H^+ + 2e $
extern EIRENE R2_1_5_H10; //$ e + H = H^+ + 2e $
extern EIRENE R2_1_8_H4;  //$ H^+ + e = H(1s) $
extern EIRENE R2_1_8_H10; //$ H^+ + e = H(1s) $
extern EIRENE R2_2_5g_H4; //$ e + H_2 = e + H + H  $
extern EIRENE R2_2_10_H4; //$ e + H_2 = 2e + H + H^+  $
extern EIRENE R3_2_3_H2;  //$ p + H_2 = H + H_2^+  $
extern EIRENE R3_2_3_H3;  //$ p + H_2(v) = H + H_2^+  $
extern EIRENE R3_1_8_H3;  //$ p + H(1s) = H(1s) + p   $
extern EIRENE R0_3T_H3;   //$ p + H_2 = p + H_2 ,\ $  elastic
extern EIRENE R0_3D_H3;   //$ p + H_2 = p + H_2 ,\ $  elastic
extern EIRENE R2_2_12_H4; //$ e + H_2^+ = e + H + H^+  $
extern EIRENE R2_2_11_H4; //$ e + H_2^+ = 2e + H^+ + H^+ $
extern EIRENE R2_2_17_H2; //
extern EIRENE R7_2_3a_H4;
extern EIRENE R7_2_3b_H4;

extern EIRENE R5_3_1_H3; //$ He+ + He = He + He+

extern EIRENE RM_2_1;
extern EIRENE RM_2_2;
extern EIRENE RM_2_3;
extern EIRENE RM_2_4;
extern EIRENE RM_2_5;
extern EIRENE RM_2_8;
extern EIRENE RM_2_9;
extern EIRENE RM_2_10_1;
extern EIRENE RM_2_10_2;
extern EIRENE RM_2_11;
extern EIRENE RM_2_12;
extern EIRENE RM_2_13;
extern EIRENE RM_2_14;
extern EIRENE RM_2_23;
extern EIRENE RM_3_1_1;
extern EIRENE RM_3_1_2;
extern EIRENE RM_3_1_3;
extern EIRENE RM_3_1_4;
extern EIRENE RM_3_2;
extern EIRENE RM_2_15_1;
extern EIRENE RM_2_15_2;
extern EIRENE RM_2_16_1;
extern EIRENE RM_2_16_2;
extern EIRENE RM_2_17_1;
extern EIRENE RM_2_17_2;
extern EIRENE RM_2_18_1;
extern EIRENE RM_2_18_2;
extern EIRENE RM_2_19_1;
extern EIRENE RM_2_19_2;
extern EIRENE RM_2_20;
extern EIRENE RM_2_21;
extern EIRENE RM_2_22;

extern EIRENE R_H_H;
extern EIRENE R_H2_H2;
extern EIRENE R_H_H2;
extern EIRENE R_H2_H;

extern vector<WallEro> AlltoWall_;
extern vector<WallEro> CXtoWall_;
extern vector<WallEro> RectoWall_;
extern vector<WallEro> DisstoWall_;
extern vector<WallEro> AlltoTarget_;
extern vector<WallEro> CXtoTarget_;
extern vector<WallEro> RectoTarget_;
extern vector<WallEro> DisstoTarget_;
extern vector<WallEro> AlltoPlasmaBoundary_;
extern vector<WallEro> CXtoPlasmaBoundary_;
extern vector<WallEro> RectoPlasmaBoundary_;
extern vector<WallEro> DisstoPlasmaBoundary_;

extern WallImpactTracker D_W_sputtered;

extern GRID Grid1;
extern Grid_extern Grid3;
extern TriMesh Grid4;