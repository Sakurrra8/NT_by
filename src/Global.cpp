#include "Global.h"
#include "Particle.h"

// 1. 定义并初始化全局引擎
// 使用 random_device 产生种子
std::random_device rd;
std::mt19937 g_gen(rd());

/// path initialization
string Inputpath;            // path for Inputfiles
string Casepath;             //
string Outputpath;           // path for result Output
string name_Xlog, name_Vlog; // path for log Output

/// OP
bool K_log;               // OP: Particle way log
bool StepLog;             // OP: step log
bool K_H2_elastic;        // OP: molecule elastic collision
bool K_EcrossBDrift;      // E cross B Drift
int K_Recyc;              // OP: Plasma Recycling from Target
int K_Rec;                // OP: Plasma Recombination from Target
int K_Pump;               // OP: particle pump
int K_Methane;            // OP: CD4 pump near Target
int K_C;                  // OP: C pump near Target
int K_Ar;                 // OP: Ar pump near Target
int K_database;           // OP: 1: ADAS, 2: EIRENE
int K_database_Pra;       // OP: 1: ADAS, 2: EIRENE for D
int K_Maxwell;            // OP: 1: standard Maxwell, 2: fast Maxwell
int K_Ei;                 // OP: 1: use SOLPS-ITER EI, 2: use Te and Ti, 3: use Ti;
int K_WallRefl;           // OP: 1: Trim database, 2: Fixed parameters
int K_backScatter;        // OP: 1: Chance backing-scatter; 2: Real backing-scatter;
int K_MarColl;            // OP: 1: 2 Separate sub-reactions; 2: 1 Sum-reaction;
int K_Prob;               // OP: 1: add the cross-section; 2: add the collision probability
int K_CX_impurity;        // OP: CX collision recation
std::array<int, 5> K_Puff{};  // OP: Puff, 0 and 1 for different puff ports
int K_D2Flight;           // OP: D2 flight
int K_abnormal_transport; // OP: abnormal transport
int K_flight;             // OP: 1 for mean free path; 2 for fixed time step; 3 for
int K_Reflect;            // OP: 1 for empirical formula; 2 for Trim database; 3 for calculate from Trim
int K_ReflectDirection;   // OP: 1 for cosing distribution; 2 for forward-reflect distribution; 3 for mirror reflection
bool backGridBoundry;     // OP: make Particles go back to Grid when Particle flight
bool K_H;                 // OP: H calculation
bool K_D;                 // OP: D calculation
bool K_T;                 // OP: T calculation
int K_back;               // OP: plasma backgroung. 1 for H, 2 for D, 3 for T
int K_dn;                 // OP: 1 for fixed 0.3; 2 file for reading from SOLPS
double Ratio_T;           // OP: n_T / n_D
int K_CX_DT;              // OP: charge exchange collision between D and T
int K_DT;                 // OP: 1 for CFEDR D-T simulation;
int K_Wallelement;        // OP: 1 for W, 2 for C
int K_T_array = 1;        // OP: energy hist in grid

/// simulate parameter
double dt;                // Time step
double t_max;             // Max time for Particle
int numPar_flight;        // number of test flight each Grid
int numPar_flight_D2;     //
int numPar_flight_T2;     //
int numPar_flight_CD4;    // number of test flight of Methane
int numPar_flight_Target; // sum of test flight from Wall recycling
double Num_CD4_pump;      //
double Num_D2_pump;       //
double Num_T2_pump;       // number of test flight of pumping particle

int K_GRID = 1;
int K_test1 = 0;
int K_test2 = 0;
int K_test3 = 0;
int K_CD4;
int K_Tn = 1;       // 1 for relative speed energy for collision, 2 for 1.5*Tn for speed;
int K_Vi = 1;       // 1 for Ti maxwell and calculate the relative v; 2 for don't calculate the relative v;
int K_mu = 3;       // 3 is the default option;
int K_PartoPar = 0; // OP: if the PartoPar is working
int K_NNCs = 0;
int K_Roulette = 0;
int K_Splitting = 0;
int K_H5Output = 1;
double W_RouletteMin = 0.05;
double W_RouletteTarget = 0.2;
double W_SplitMax = 5.0;
double W_SplitTarget = 1.0;
int MaxSplit = 4;

int MeshMode; // Option 1: Orthogonal grid; 2: Full grid; 3: Triangular grid

double Dn[N_RADIAL_GRID];

/// @brief Plasma and device parameter
double InterscePoint[10][6] = {0}; // 0: distance of particle and interscePoint; 1,2: x,y index; 3: number of boundary; 4:types of intersection boundaries
double NumPar_now = 0;
double T_N = 0.02;    // Temperature of injected impurities
double T_wall = 0.02; // Temperture of Target(eV)

// 真空区域参数
double Te_vacuum = 0.3, ne_vacuum = 1e13;
// 核心区域参数
double Te_core = 1200, ne_core = 9e19, Tn_core = 1200;

// double cosang[410], sinang[410];
// double cosCore[60], sinCore[60];
double coefficient_D = 2., coefficient_T = 3.;

/// @brief Variable about wall
int num_CoreBoundry = 48;
int num_GridBoundry = N_poloidal - 2;
int num_PFRBoundry = 48;
double CoreBoundry[55][2], GridBoundry[100][2], PFRBoundry[55][2];

/// @brief Recycling coefficient at wall and target
double coeff_recyc_target, coeff_ercyc_wall, coeff_puff;

/// @brief 2D Grid Value
int N_poloidal = N_POLOIDAL_GRID;
int N_radial = N_RADIAL_GRID;

// 三维网格
std::vector<std::vector<std::vector<double>>> Grid;
std::vector<std::vector<std::vector<double>>> B;

// 二维场量
std::vector<std::vector<double>> ne, Te, Ti, Te_coll, Ti_coll;

// 1 原子态密度
std::vector<std::vector<double>> n_H_1, n_D_1, n_T_1;
std::vector<double> n_H_1_Tri, n_D_1_Tri, n_T_1_Tri;

// H 相关
std::vector<double> n_H_0_Tri, T_H_0_Tri, T_H2_0_Tri, n_H2_0_Tri;
std::vector<std::vector<double>> ua_H_0_Tri, ua_H2_0_Tri;
std::vector<std::vector<double>> n_H_0, T_H_0, T_H2_0, n_H2_0;
std::vector<std::vector<double>> n_H2_1, ua_H_1;
std::array<double, 3> V_H_1_now, V_H_0_now, V_H2_0_now;
std::vector<std::vector<double>> Vi_H, Ti_H_thermal;

// D 相关
std::vector<double> n_D_0_Tri, T_D_0_Tri, T_D2_0_Tri, n_D2_0_Tri;
std::vector<std::vector<double>> ua_D_0_Tri, ua_D2_0_Tri;
std::vector<std::vector<double>> n_D_0, T_D_0, T_D2_0, n_D2_0;
std::vector<std::vector<double>> n_D2_1, ua_D_1;
std::array<double, 3> V_D_1_now, V_D_0_now, V_D2_0_now;
std::vector<std::vector<double>> Vi_D, Ti_D_thermal;

// T 相关
std::vector<double> n_T_0_Tri, T_T_0_Tri, T_T2_0_Tri, n_T2_0_Tri;
std::vector<std::vector<double>> ua_T_0_Tri, ua_T2_0_Tri;
std::vector<std::vector<double>> n_T_0, T_T_0, T_T2_0, n_T2_0;
std::vector<std::vector<double>> n_T2_1, ua_T_1;
std::array<double, 3> V_T_1_now, V_T_0_now, V_T2_0_now;
std::vector<std::vector<double>> Vi_T, Ti_T_thermal;

// 几何参数
std::vector<std::vector<double>> Volume, Area;

// 磁场
std::vector<std::vector<double>> Btot, Btor, Brad, Bpol;

// 电场与漂移速度
std::vector<std::vector<double>> Erad, Epol, Etor, vdEp, vdEr;
std::vector<std::vector<double>> epx, epy, erx, ery;

// 其它物理量分量
std::vector<std::vector<double>> crx, cry;
std::vector<std::vector<double>> dTip, dTir, dTix, dTiy, dTiz;

/// @brief Target Value
std::vector<double> NumPar_H_recyc, Tn_H_recyc, NumPar_H2_recyc;
std::vector<double> NumPar_D_recyc, Tn_D_recyc, NumPar_D2_recyc;
std::vector<double> NumPar_T_recyc, Tn_T_recyc, NumPar_T2_recyc;
double NumPar_wall_H[80], NumPar_wall_D[80], NumPar_wall_T[80];
double Xrecyc[80][2], theta[80][2], Ei_Dion[80];

/// @brief variable for cross-CX
double Ratio_D_Coll, Ratio_T_Coll, Ratio_DT_Coll;
double HHmass, DTmass, D2Tmass, DT2mass, HH2mass, DD2mass, TT2mass, H2H2mass, D2D2mass, T2T2mass;

/// variable for debug
int IfOut = 0;        // variable for judging if patricle escape
int Num_Reflect;      // Number of particle reflections
int Num_Reflect_Core; // Number of Charge particle reflections from Core

/// @brief backscattering Coefficient fitting when D on W from Trim database
double H_W_n_a[4] = {0.4116e0, -0.9122e-1, 0.6956e0, 0.1357e1};                      //
double H_W_E_a[4] = {0.2039e0, -0.1501e0, 0.1064e1, 0.1359e1};                       //
double H_W_n_b[6] = {0.2393e0, 0.7644e-1, 0.8777e0, 0.1444e1, 0.7384e1, 0.1545e0};   //
double H_W_E_b[6] = {0.1923e0, -0.1605e0, 0.1187e1, 0.1582e1, 0.1424e0, 0.1044e1};   //
double H_W_epsilon_L = 9.86986e3;                                                    //
double D_W_n_a[4] = {0.4453e0, -0.8079e-1, 0.5964e0, 0.1344e1};                      //
double D_W_E_a[4] = {0.2547e0, -0.1146e0, 0.1195e1, 0.1259e1};                       //
double D_W_n_b[6] = {0.3246e0, -0.8766e-1, 0.5956e0, 0.1356e1, 0.3579e1, 0.6475e-1}; //
double D_W_E_b[6] = {0.2040e0, -0.1114e0, 0.1204e1, 0.1420e1, 0.4466e0, 0.3570e0};   //
double D_W_epsilon_L = 9.92326e3;                                                    //
double T_W_n_a[4] = {0.4632e0, -0.6965e-1, 0.5252e0, 0.1322e1};                      //
double T_W_E_a[4] = {0.2547e0, -0.1113e0, 0.8170e0, 0.1354e1};                       //
double T_W_n_b[6] = {0.2408e0, 0.4954e-1, 0.5626e0, 0.1402e1, 0.9493e1, -0.3282e-1}; //
double T_W_E_b[6] = {0.1253e0, -0.1479e0, 0.8593e0, 0.1429e1, 0.9369e1, 0.2049e0};   //
double T_W_epsilon_L = 9.97718e3;                                                    //

double H_C_n_a[4] = {0.2119e0, -0.2217e0, 0.2461e0, 0.1316e1};                        //
double H_C_E_a[4] = {0.8723e-1, -0.2768e0, 0.3988e0, 0.1282e1};                       //
double H_C_n_b[6] = {0.2016e0, -0.2404e0, 0.1888e0, 0.1392e1, 0.2136e-2, 0.1082e1};   //
double H_C_E_b[6] = {0.7983e-1, -0.3088e0, 0.2742e0, 0.1368e1, 0.8600e-3, 0.1261e1};  //
double H_C_epsilon_L = 4.14659e2;                                                     //
double D_C_n_a[4] = {0.1526e0, -0.2304e0, 0.2113e0, 0.1287e1};                        //
double D_C_E_a[4] = {0.5142e1, -0.2714e0, 0.2668e0, 0.1316e1};                        //
double D_C_n_b[6] = {0.1375e0, -0.2648e0, 0.1494e0, 0.1467e1, 0.6183e-1, 0.7212e0};   //
double D_C_E_b[6] = {0.4806e-1, -0.2963e0, 0.1866e0, 0.1410e1, 0.1778e-2, 0.1109e1};  //
double D_C_epsilon_L = 4.46507e2;                                                     //
double T_C_n_a[4] = {0.1044e0, -0.2618e0, 0.2052e0, 0.1228e1};                        //
double T_C_E_a[4] = {0.2889e-1, -0.2868e0, 0.2312e0, 0.1293e1};                       //
double T_C_n_b[6] = {0.4854e-1, -0.1796e0, 0.2203e0, 0.1280e1, 0.4869e2, -0.3260e-1}; //
double T_C_E_b[6] = {0.1738e-1, -0.1534e0, 0.2513e0, 0.1446e1, 0.2318e1, 0.9630e-1};  //
double T_C_epsilon_L = 4.78673e2;                                                     //

double H_Be_n_a[4] = {0.1878e0, -0.2469e0, 0.2448e0, 0.1271e1};                       //
double H_Be_E_a[4] = {0.7524e-1, -0.2871e0, 0.4301e0, 0.1233e+1};                     //
double H_Be_n_b[6] = {0.1652e0, -0.2991e0, 0.1196e0, 0.1374e1, 0.4694e-4, 0.1476e1};  //
double H_Be_E_b[6] = {0.6656e-1, -0.3327e0, 0.2542e0, 0.1342e1, 0.1402e-2, 0.1149e1}; //
double H_Be_epsilon_L = 2.56510e2;                                                    //
double D_Be_n_a[4] = {0.1169e0, -0.2841e0, 0.2138e0, 0.1205e1};                       //
double D_Be_E_a[4] = {0.3521e-1, -0.3095e0, 0.2794e0, 0.1243e+1};                     //
double D_Be_n_b[6] = {0.9725e-1, -0.2597e0, 0.1623e0, 0.1446e1, 0.2505e0, 0.4941e0};  //
double D_Be_E_b[6] = {0.3305e-1, -0.3314e0, 0.1982e0, 0.1304e1, 0.4603e-3, 0.1178e1}; //
double D_Be_epsilon_L = 2.82110e2;                                                    //
double T_Be_n_a[4] = {0.7312e-1, -0.3123e0, 0.2143e0, 0.1150e1};                      //
double T_Be_E_a[4] = {0.1725e-1, -0.2993e0, 0.2554e0, 0.1250e+1};                     //
double T_Be_n_b[6] = {0.5283e-1, -0.3846e0, 0.1825e0, 0.1171e1, 0.6160e2, 0.3601e-1}; //
double T_Be_E_b[6] = {0.7229e-2, -0.2502e0, 0.3048e0, 0.1415e1, 0.3736e1, 0.1321e0};  //
double T_Be_epsilon_L = 3.07966e2;                                                    //

/// @brief other Value
std::vector<double> V_1, V_2;

/// @brief area of the target and wall
vector<double> S_target; // area of the target
vector<double> S_wall;   // area of the wall

/// Particle Change
PartoPar PartoPar_Temp;

/*extend_mesh*/
vector<vector<double>> fnax_e;
vector<vector<double>> ne_e;
vector<vector<double>> ni_e;
vector<vector<double>> Te_e;
vector<vector<double>> Ti_e;
vector<vector<double>> dTip_e;
vector<vector<double>> dTir_e;
vector<vector<double>> ua_e;
vector<vector<double>> duap_e;
vector<vector<double>> Ex_e;
vector<vector<double>> Ey_e;
vector<vector<double>> Ez_e;
vector<vector<double>> Bp_e;
vector<vector<double>> Br_e;
vector<vector<double>> Bz_e;
vector<vector<double>> BB_e;
vector<vector<double>> dBp_e;
vector<vector<double>> eBx_e;
vector<vector<double>> eBy_e;
vector<vector<double>> eBz_e;

/// ADAS class
ADAS SCD12_H("scd12_h.dat"); // Effective Ionisation Coefficients
ADAS ACD12_H("acd12_h.dat"); // Effective Recombination Coefficients
ADAS CCD96_D("ccd96_d.dat"); // D-D Effective CX Coefficients
ADAS PLT12_H("plt12_h.dat"); // Line Emission from excitation
ADAS PRB12_H("prb12_h.dat"); // Recombination and Bremsstrahlung Radiation
ADAS PRC96_D("prc96_d.dat"); // D-D Effective CX RadiationADAS SCD96_He("scd96_he.dat"); // Effective Ionisation Coefficients
ADAS SCD96_T("scd96_t.dat");
ADAS ACD96_T("acd96_t.dat");
ADAS CCD96_T("ccd96_t.dat");
ADAS CCD96_H("ccd96_h.dat");
ADAS PLT96_T("plt96_t.dat");
ADAS PRB96_T("prb96_t.dat");
ADAS PRC96_T("prc96_t.dat");
ADAS PRC96_H("prc96_h.dat");
ADAS PLS96_T("pls96_t.dat");
ADAS PLS89_H("pls89_h.dat");

ADAS SCD96_He("scd96_he.dat"); // Effective Ionisation Coefficients
ADAS ACD96_He("acd96_he.dat"); // Effective Recombination Coefficients
ADAS CCD96_He("ccd96_he.dat"); // D-D Effective CX Coefficients
ADAS PLT96_He("plt96_he.dat"); // Line Emission from excitation
ADAS PRB96_He("prb96_he.dat"); // Recombination and Bremsstrahlung Radiation
ADAS PRC96_He("prc96_he.dat"); // D-D Effective CX Radiation

ADAS SCD96_C("scd96_c.dat");
ADAS ACD96_C("acd96_c.dat");
ADAS CCD96_C("ccd96_c.dat");
ADAS PLT96_C("plt96_c.dat"); // Line Emission from excitation
ADAS PRB96_C("prb96_c.dat"); // Recombination and Bremsstrahlung Radiation
ADAS PRC96_C("prc96_c.dat"); // D-D Effective CX Radiation

ADAS SCD89_Ar("scd89_ar.dat");
ADAS ACD89_Ar("acd89_ar.dat");
ADAS CCD89_Ar("ccd89_ar.dat");
ADAS PLT89_Ar("plt89_ar.dat"); // Line Emission from excitation
ADAS PRB89_Ar("prb89_ar.dat"); // Recombination and Bremsstrahlung Radiation
ADAS PRC89_Ar("prc89_ar.dat"); // D-D Effective CX Radiation

ADAS SCD89_B("scd89_b.dat");
ADAS ACD89_B("acd89_b.dat");

ADAS SCD96_Li("scd96_li.dat");
ADAS ACD96_Li("acd96_li.dat");

ADAS PEC_Li1;
ADAS PEC_Li2;
ADAS PEC_B_1;

ADAS PEC_W_1;
ADAS PEC_W_8;
ADAS PEC_W_11;
ADAS PEC_W_15;
ADAS PEC_W_22;
ADAS PEC_W_24;

Reflect H_W;
Reflect D_W;
Reflect T_W;
Reflect H_C;
Reflect D_C;
Reflect T_C;
Reflect H_Be;
Reflect D_Be;
Reflect T_Be;
Reflect Ar_W;
Reflect He_W;

std::ofstream Xlog;
std::ofstream Vlog;

std::ofstream out;

// Particle
Particle D("D", Dmass, 1, 3);
Particle D2("D2", D2mass, 1, 4);
Particle CD4("CD4", CD4mass, 1, 5);
Particle CD3("CD3", CD3mass, 1, 6);
Particle CD2("CD2", CD2mass, 1, 7);
Particle CD1("CD1", CD1mass, 1, 8);
Particle C("C", Cmass, 6, 8);
Particle Ar("Ar", Armass, 18, 10);
Particle H("H", Hmass, 1, 13);
Particle H2("H2", H2mass, 1, 14);
Particle T("T", Tmass, 1, 15);
Particle T2("T2", T2mass, 1, 16);
Particle *P;

/// EIRENE
EIRENE R3_2_3r_H4(EIRENE_H4, 6710);  // Mar via H2+
EIRENE R3_2_3d_H4(EIRENE_H4, 6767);  // Mar via H2+
EIRENE R3_2_3i_H4(EIRENE_H4, 6827);  // Mar via H2+
EIRENE R2_2_17r_H4(EIRENE_H4, 7007); // Mar via H-
EIRENE R2_2_17d_H4(EIRENE_H4, 7061); // Mar via H-
EIRENE R2_2_5_H4(EIRENE_H4, 4511);   //$ e + H_2 = e + H + H$
EIRENE R2_2_14_H2(EIRENE_H2, 2157);  //$ e + H_2^+ = H + H  $
EIRENE R2_2_14_H4(EIRENE_H4, 4844);  //$ e + H_2^+ = H + H  $
EIRENE R2_2_14_H8(EIRENE_H8, 7147);  //$ e + H_2^+(v) = H(1s) + H^*(n)   (v=0-9, n \ge 2) $
EIRENE R2_2_9_H4(EIRENE_H4, 4628);   //$ e + H_2 = 2e + H_2^+  $
EIRENE R2_1_5_H4(EIRENE_H4, 4156);   //$ H + e = H^+ + 2e $
EIRENE R2_1_5_H10(EIRENE_H10, 7843); //$ e + H = H^+ + 2e $
EIRENE R2_1_8_H4(EIRENE_H4, 4329);   //$ H^+ + e = H(EIRENE_1s) $
EIRENE R2_1_8_H10(EIRENE_H10, 7962); //$ H^+ + e = H(EIRENE_1s) $
EIRENE R2_2_5g_H4(EIRENE_H4, 4571);  //$ e + H_2 = e + H + H  $
EIRENE R2_2_10_H4(EIRENE_H4, 4679);  //$ e + H_2 = 2e + H + H^+  $
EIRENE R3_2_3_H2(EIRENE_H2, 2253);   //$ p + H_2 = H + H_2^+  $

EIRENE R3_2_3_H3(EIRENE_H3, 3841); //$ p + H_2(EIRENE_v) = H + H_2^+  $
EIRENE R3_1_8_H3(EIRENE_H3, 3683); //$ p + H(EIRENE_1s) = H(EIRENE_1s) + p   $
EIRENE R0_3T_H3(EIRENE_H3, 2972);  //$ p + H_2 = p + H_2 ,\ $  elastic
EIRENE R0_3D_H3(EIRENE_H3, 3026);  //$ p + H_2 = p + H_2 ,\ $  elastic

EIRENE R2_2_12_H4(EIRENE_H4, 4788); //$ e + H_2^+ = e + H + H^+  $
EIRENE R2_2_11_H4(EIRENE_H4, 4732); //$ e + H_2^+ = 2e + H^+ + H^+ $
EIRENE R2_2_17_H2(EIRENE_H2, 2188); // $ e + H_2^+ = H + H^- $
EIRENE R7_2_3a_H4(EIRENE_H4, 6888);
EIRENE R7_2_3b_H4(EIRENE_H4, 6948);

EIRENE R5_3_1_H3(EIRENE_H3, 5533, "hydhel"); //$ He+ + He = He + He+

EIRENE RM_2_1(EIRENE_H2, 1020, "methane");
EIRENE RM_2_2(EIRENE_H2, 1031, "methane");
EIRENE RM_2_3(EIRENE_H2, 1048, "methane");
EIRENE RM_2_4(EIRENE_H2, 1059, "methane");
EIRENE RM_2_5(EIRENE_H2, 1076, "methane");
EIRENE RM_2_8(EIRENE_H2, 1093, "methane");
EIRENE RM_2_9(EIRENE_H2, 1104, "methane");
EIRENE RM_2_10_1(EIRENE_H2, 1126, "methane");
EIRENE RM_2_10_2(EIRENE_H2, 1137, "methane");
EIRENE RM_2_11(EIRENE_H2, 1154, "methane");
EIRENE RM_2_12(EIRENE_H2, 1165, "methane");
EIRENE RM_2_13(EIRENE_H2, 1182, "methane");
EIRENE RM_2_14(EIRENE_H2, 1193, "methane");
EIRENE RM_2_15_1(EIRENE_H2, 1214, "methane");
EIRENE RM_2_15_2(EIRENE_H2, 1225, "methane");
EIRENE RM_2_16_1(EIRENE_H2, 1246, "methane");
EIRENE RM_2_16_2(EIRENE_H2, 1257, "methane");
EIRENE RM_2_17_1(EIRENE_H2, 1278, "methane");
EIRENE RM_2_17_2(EIRENE_H2, 1289, "methane");
EIRENE RM_2_18_1(EIRENE_H2, 1310, "methane");
EIRENE RM_2_18_2(EIRENE_H2, 1321, "methane");
EIRENE RM_2_19_1(EIRENE_H2, 1342, "methane");
EIRENE RM_2_19_2(EIRENE_H2, 1353, "methane");
EIRENE RM_2_20(EIRENE_H2, 1370, "methane");
EIRENE RM_2_21(EIRENE_H2, 1381, "methane");
EIRENE RM_2_22(EIRENE_H2, 1393, "methane");
EIRENE RM_2_23(EIRENE_H2, 1410, "methane");
EIRENE RM_3_1_1(EIRENE_H3, 1437, "methane");
EIRENE RM_3_1_2(EIRENE_H3, 1489, "methane");
EIRENE RM_3_1_3(EIRENE_H3, 1541, "methane");
EIRENE RM_3_1_4(EIRENE_H3, 1593, "methane");
EIRENE RM_3_2(EIRENE_H3, 1645, "methane");

EIRENE R_H_H(EIRENE_H2, 1026, "ammonx_elsa"); //
EIRENE R_H2_H2(EIRENE_H2, 1033, "ammonx_elsa");
EIRENE R_H_H2(EIRENE_H2, 1040, "ammonx_elsa");
EIRENE R_H2_H(EIRENE_H2, 1050, "ammonx_elsa");

GRID Grid1(N_poloidal, N_radial);
Grid_extern Grid3;
TriMesh Grid4;

WallImpactTracker D_W_sputtered;
