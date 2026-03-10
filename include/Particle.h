#pragma once
#include "Global.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include "ADAS.h"
#include "EIRENE.h"
#include "Reflect.h"
#include "stack"
#include "queue"
#include <cassert>

using namespace std;
using namespace H5;

double Maxwell(double E, double Mass);
double Random();
void Edge(int Start[2], int End[2], double a[], double b[], int *k);
bool IfinEdge(double a[], double b[], int k, double x, double y);
int Xfind(int Start[], int End[], double x, double y);
int Yfind(int Start[], int End[], double x, double y);
void Vinit_cosine(std::vector<double> &Vt, double Rcos, double Rsin);
int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
						  double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y);

int ZonefromXY(int X, int Y);
int PointandLine(double line_x1, double line_y1, double line_x2, double line_y2, double point_x, double point_y);
string fatename(int i);
string sourcename(int i);

class ParCollCar
{
private:
	int K_;
	int e_or_i_;
	int cc;

	double Sn_[98][38];
	double Smu_[98][38];
	double Smu1_[98][38];
	double Smu2_[98][38];
	double SE_[98][38];
	double Num_n_[98][38];
	double Num_mu_[98][38];
	double Num_mu1_[98][38];
	double Num_mu2_[98][38];
	double Num_E_[98][38];
	double Pra_[98][38];
	double crossSection[98][38];
	double cs_[98][38];
	double SE_n_[98][38];
	double Num_E_n_[98][38];

	double Cor_cs_;
	double cs_now_;
	vector<double> V_relative_;
	double V_2_relative_;
	double V_relative_2_;
	double num_trimesh_;

	vector<double> Tri_cs_;
	vector<double> Tri_Sn_;
	vector<double> Tri_Smu_;
	vector<double> Tri_Smu_0_;
	vector<double> Tri_Smu_1_;
	vector<double> Tri_Smu_2_;
	vector<double> Tri_Smu1_;
	vector<double> Tri_Smu2_;
	vector<double> Tri_SE_;
	vector<double> Tri_Num_n_;
	vector<double> Tri_Num_mu_;
	vector<double> Tri_Num_mu_0_;
	vector<double> Tri_Num_mu_1_;
	vector<double> Tri_Num_mu_2_;
	vector<double> Tri_Num_mu1_;
	vector<double> Tri_Num_mu2_;
	vector<double> Tri_Num_E_;
	vector<double> Tri_Pra_;
	vector<double> Tri_crossSection;
	vector<double> Tri_SE_n_;
	vector<double> Tri_Num_E_n_;

public:
	ParCollCar() {}
	void Clear();
	int K();
	void Set_K(int i);
	void Stat(int Charge, double *n[98][38], ADAS *PraADAS, const std::vector<std::vector<double>> &ni = {});
	void Stat1(int Charge, double *n[98][38], ADAS *PraADAS, const std::vector<std::vector<double>> &ni = {});
	void StatEIRENE(int Charge, double *n[98][38], EIRENE *PraEIRENE, double E, double **ni = NULL);
	void RecStat(int Charge, double *n[98][38], ADAS *PraADAS);

	void Stat_Tri(int Charge, const std::vector<std::vector<double>> &n, ADAS *PraADAS, const std::vector<std::vector<double>> &ni = {});
	void StatEIRENE_Tri(int Charge, std::vector<std::vector<double>> &n, EIRENE *PraEIRENE, double E, const std::vector<std::vector<double>> &ni = {});
	void RecStat_Tri(int Charge, const std::vector<std::vector<double>> &n, ADAS *PraADAS);

	void OutSn(std::ofstream &in);
	void OutSmu(std::ofstream &in);
	void OutSE(std::ofstream &in);
	void OutPra(std::ofstream &in);
	void OutSmu1(std::ofstream &in);
	void OutSmu2(std::ofstream &in);
	void OutSE_n(std::ofstream &in);
	double Sn(int X, int Y);
	double Smu(int X, int Y);
	double SE(int X, int Y);
	void Setcs(int i, double cs);
	void Setcs(int i, int j, double cs);
	double cs(int XY[3], double Zone);
	double cs(int i, int j);
	double cs(int i);
	void SetCor_cs(double Cor_cs);
	double Cor_cs();
	double cs_now();

	void n_Add(int XY[], double Var);
	void Mu_Add(int XY[], double Var);
	void Mu_Add(int XY, double Var0, double Var1, double Var2);
	void Mu_Add1(int XY[], double Var);
	void Mu_Add2(int XY[], double Var);
	void E_Add(int XY[], double Var);
	void Pra_Add(int XY[], double Var);
	void E_n_Add(int XY[], double Var);
	void SetPra(int i, int j, double Var);

	void n_Add(int XY, double Var);
	void Mu_Add(int XY, double Var);
	void E_Add(int XY, double Var);
	void E_n_Add(int XY, double Var);
	void Pra_Add(int XY, double Var);
	void SetPra_Tri(int i, double Var);
	void OutSn_Tri(std::ofstream &in);
	void OutSmu_Tri(std::ofstream &in);
	void OutSE_Tri(std::ofstream &in);
	void OutPra_Tri(std::ofstream &in);

	void initialize(int e_or_i, int num_trimesh);
	void ADASInput(ADAS *ParColl, int Charge, int Par = 2, int cross_CX = 0);
	void EIRENEInput(EIRENE *ParColl);

	void Setcs_now(double cs);

	void SetTrics_cs(std::vector<std::vector<int>> &Tri_B2);
	void Set_V_relative(double Vi0, double Vi1, double Vi2, double Vn0, double Vn1, double Vn2);
	double V_relative(int i);
	double V_2_relative();
	double V_relative_2();
};

// 简单的3D向量结构体
struct Vector3
{
	double x, y, z;

	Vector3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

	// 加法
	Vector3 operator+(const Vector3 &other) const
	{
		return Vector3(x + other.x, y + other.y, z + other.z);
	}

	// [新增] 减法 (计算向量差 p1 - p0 需要)
	Vector3 operator-(const Vector3 &other) const
	{
		return Vector3(x - other.x, y - other.y, z - other.z);
	}

	// 数乘
	Vector3 operator*(double scalar) const
	{
		return Vector3(x * scalar, y * scalar, z * scalar);
	}

	// [新增] 点积 (计算夹角和反射需要)
	double dot(const Vector3 &other) const
	{
		return x * other.x + y * other.y + z * other.z;
	}

	// 模长
	double norm() const
	{
		return std::sqrt(x * x + y * y + z * z);
	}

	// 叉乘
	Vector3 cross(const Vector3 &other) const
	{
		return Vector3(
			y * other.z - z * other.y,
			z * other.x - x * other.z,
			x * other.y - y * other.x);
	}

	// 归一化
	Vector3 normalize() const
	{
		double n = norm();
		if (n == 0)
			return Vector3(0, 0, 0);
		return Vector3(x / n, y / n, z / n);
	}
};

Vector3 calculate_dissociation_velocity(Vector3 v_mol, double E_diss, double m_h = 1.673e-27);

/**
 * @brief  记录跨区域输运：粒子数、动量、能量及能谱
 * @note   判区逻辑由调用方提供；当检测到区域 from→to 时，调用 accumulate()
 */
class FluxTracker
{
public:
	FluxTracker();
	FluxTracker(int nRegion,
				int nBin,
				double eMin,
				double eMax,
				bool trackMom = true,
				bool trackPow = true);
	bool FluxTraInit(int nRegion, int nBin,
					 double eMin, double eMax,
					 int nx, int ny,
					 bool trackMom, bool trackPow);

	/// 记录一次跨区事件
	/// @param from,to  区域索引
	/// @param E        粒子能量 (eV)
	/// @param vx,vy,vz 速度分量 (m/s)；若不统计动量可传 0
	/// @param w        权重
	void accumulate(int from, int to,
					double E,
					double vx, double vy, double vz,
					double w);

	/// 在(i,j)网格位置统计能谱、动量和能量
	/// @param i,j 网格坐标
	/// @param E 粒子能量 (eV)
	/// @param vx,vy,vz 速度分量 (m/s)
	/// @param w 粒子权重
	void accumulateGrid(int i, int j,
						double E,
						double vx, double vy, double vz,
						double w);
	void normalizeByVolume();

	/// 输出到 HDF5：/flux, /hist, (/momFlux), (/powFlux)
	void write_H5(const std::string &fileName,
				  const std::string &group = "/flux_stat") const;

private:
	int NREG_, NBIN_;
	double E_MIN_, E_MAX_;
	bool trackMom_, trackPow_; // 网格维度
	int NX_{0}, NY_{0};
	bool gridEnabled_{false}; // 是否开启网格统计// 网格能谱：NX × NY × NBIN
	std::vector<std::vector<std::vector<double>>> gridHist_;

	// 网格动量：NX × NY × 3
	std::vector<std::vector<std::array<double, 3>>> gridMom_;

	// 网格能量：NX × NY
	std::vector<std::vector<double>> gridPow_;

	std::vector<std::vector<double>> flux_;					  // NREG×NREG
	std::vector<std::vector<std::vector<double>>> hist_;	  // NREG×NREG×NBIN
	std::vector<std::vector<std::array<double, 3>>> momFlux_; // (可选)
	std::vector<std::vector<double>> powFlux_;				  // (可选)
	bool inited_{false};									  // 标记是否已初始化
};

class Particle
{
private:
	/// @brief fixed property of this particle
	string name_;
	double mass_;
	double MaxCharge_;
	vector<double> Pump_;
	int Index_;

	int num_trimesh_;
	int Tri_Index_;
	std::vector<std::vector<int>> Tri_B2_;
	std::vector<std::vector<double>> Tri_Totalcs_;
	std::vector<std::vector<double>> Tri_CollProb_;
	std::vector<std::vector<double>> Tri_lambda_;

	/// @brief Properties of particles in flight
	int Charge_;
	int ChargeTag_;
	double X_[3];			// position
	std::vector<double> V_; // velocity
	double V_Charge_[8];	// velocity when Particle Charge != 0
	double X_new_[3];		// position at next time step
	double dt_;				// the dt calculated from lambda_min
	double dt_trace_;		// the dt calculated from trace
	int Zone_;				// Particle Zone in Tokamak
	double Tn_;				// Particle Temperture
	double Weight_;			// Particle weight
	int XY_[3];				// Particle position Grid Index in Tokamak
	int GridIndex_;			// the mesh Index

	/// @brief Temporary attribute during flight
	int boundary_start_;   // Particle start from which boundry of this grid
	int IfColl_;		   // If the particle Coll in here
	int IfFlightOut_;	   //
	int fate_[2];		   // fate[0] is source, fate[1] is sink
	int sourcePar_[2];	   //
	int sourceGrid_[2];	   // the grid that the particle from
	double sourceWall_[4]; //
	double Xold_[3];	   // position at previous time step
	int Zone_old_;		   // Particle Zone in Tokamak(X)
	int GridIndex_old_;	   // GridIndex = XY_[0] * 38 + XY_[1]
	double *lambda_min_;
	double lambda_now_;
	double d_flight_;
	double Rand_flight_;
	double *lambda_[98][38];

	/// @brief Initial launch parameters
	int Num_Target_[80];		 // for Particle from Recycling
	double NumPar_Target_[80];	 // for Particle from Recycling
	double Weight_Target_[80];	 // for Particle from Recycling
	double NumPar_sum_Target_;	 // for Particle from Recycling
	double NumPar_Grid_[98][38]; // for Particle from Grid Collision
	double Weight_Grid_[98][38]; // for Particle from Grid Collision
	double NumPar_sum_Grid_;	 // for Particle from Grid Collision

	/// @brief Statistical parameters
	double *n_[98][38];
	double *T_[98][38];
	double *E_[98][38];
	double *V_Grid_[98][38][4];
	double *V_D_1_[98][38][4];
	double *Num_V_D_1_[98][38];
	double *Sum_n_[98][38];
	double *Sum_E_[98][38];
	double *Sum_V_[98][38][4];
	double *Sum_V_D_1_[98][38][4];

	/// @brief Triangle Statistical parameters
	std::vector<std::vector<double>> Tri_n_;
	std::vector<std::vector<double>> Tri_T_;
	std::vector<std::vector<double>> Tri_E_;
	std::vector<std::vector<std::vector<double>>> Tri_V_Grid_;
	std::vector<std::vector<std::vector<double>>> Tri_V_D_1_;
	std::vector<std::vector<double>> Tri_Num_V_D_1_;
	std::vector<std::vector<double>> Tri_Sum_n_;
	std::vector<std::vector<double>> Tri_Sum_E_;
	std::vector<std::vector<std::vector<double>>> Tri_Sum_V_;
	std::vector<std::vector<std::vector<double>>> Tri_Sum_V_D_1_;
	std::vector<double> Tri_NumPar_Grid_;

	double *V_Grid_CX_Ion_Be_[98][38][4];
	double *V_Grid_CX_Ion_Af_[98][38][4];
	double *V_Grid_CX_Neu_Be_[98][38][4];
	double *V_Grid_CX_Neu_Af_[98][38][4];
	double *Sum_V_CX_Ion_Be_[98][38][4];
	double *Sum_V_CX_Ion_Af_[98][38][4];
	double *Sum_V_CX_Neu_Be_[98][38][4];
	double *Sum_V_CX_Neu_Af_[98][38][4];

	vector<vector<double>> Sum_n_array_[98][38];
	vector<double> T_array_;
	int Num_array_;

	/*double Flux_GridBoundry[97][2]; // Inner Flux when second index is 0, Out Flux when second index is 1
	double Flux_CoreBoundry[49][2];
	double Flux_PFRBoundry[49][2];
	double Flux_Target_l[36][2];
	double Flux_Target_r[36][2];
	double EneFlux_GridBoundry[97][2]; // Inner Flux when second index is 0, Out Flux when second index is 1
	double EneFlux_CoreBoundry[49][2];
	double EneFlux_PFRBoundry[49][2];
	double EneFlux_Target_l[36][2];
	double EneFlux_Target_r[36][2];*/

	/// @brief collision parameters
	vector<double> CoreTotalcs_;
	vector<double> CoreCollProb_;
	vector<double> CollProb_[98][38];
	vector<double> Totalcs_[98][38];

	vector<vector<vector<double>>> n_Flux_Grid_;
	vector<vector<vector<double>>> T_Flux_Grid_;
	vector<vector<double>> First_CX_;
	FluxTracker FT_;

public:
	Particle();
	Particle(string particle_name, double mass, int MaxCharge, int Index);
	~Particle();
	void ParInit(vector<int> &K_collision, std::vector<std::vector<int>> Tri_B2);

	/// @brief fixed property of this particle
	vector<ParCollCar> Ion_;
	vector<ParCollCar> Rec_;
	vector<ParCollCar> CX_;
	vector<ParCollCar> CX_DT_;
	vector<ParCollCar> Ela_;
	vector<ParCollCar> Ela_DT_;
	vector<ParCollCar> MAR_;
	vector<ParCollCar> Diss1_;
	vector<ParCollCar> Diss2_;
	vector<ParCollCar> Diss3_;
	vector<vector<ParCollCar>> DS_;
	vector<ParCollCar> R_with_H_;
	vector<ParCollCar> R_with_H2_;

	double Ion_factor;
	double CX_factor;
	double CX_DT_factor;
	double Rec_factor;
	double MAR_factor;
	double Diss1_factor;
	double Diss2_factor;
	double Diss3_factor;
	double Ela_factor;
	double Ela_DT_factor;
	double DS_factor[4];
	double R_with_H_factor;
	double R_with_H2_factor;

	vector<WallEro> AlltoWall_;
	vector<WallEro> CXtoWall_;
	vector<WallEro> RectoWall_;
	vector<WallEro> DisstoWall_;
	vector<WallEro> AlltoTarget_;
	vector<WallEro> CXtoTarget_;
	vector<WallEro> RectoTarget_;
	vector<WallEro> DisstoTarget_;
	vector<WallEro> AlltoPlasmaBoundary_;
	vector<WallEro> CXtoPlasmaBoundary_;
	vector<WallEro> RectoPlasmaBoundary_;
	vector<WallEro> DisstoPlasmaBoundary_;

	/*double MAR_cs_[98][38], MAR_Cor_cs_;	 // for D2
	double Diss1_cs_[98][38], Diss1_Cor_cs_;	 // for D2
	double Diss2_cs_[98][38], Diss2_Cor_cs_; // for D2
	double Diss3_cs_[98][38], Diss3_Cor_cs_; // for D2
	double DS_cs_[98][38][4], DS_Cor_cs_[4];*/

	void Particlefrom(Particle *A, double K = 1, int Charge = -100);
	void CalWeight1(int num);
	void CalWeight2(double NumPar[76], int num);
	void Init(int k, int z = 0);
	void Vinit_neutralReflect(std::vector<double> &Vt, double Rcos, double Rsin);
	void Vinit_ionReflect(std::vector<double> &Vt, double Rcos, double Rsin, int z);
	void VtoVcharge();
	void VchargetoV();
	void Vchargefix(); // when the charged particle flight to next grid, V_charge[0,1,2] should be fixed
	void track();
	void CalLambda();
	void Caltrace();
	void Caltrace_Tri();
	void NumParStat();
	void CalProb();
	double CollProb();
	void Coll();
	int H2PCollCal();
	void regression();
	void Setdt(double dt_);
	double FlightTime();
	void SetX(int i, double X);
	void SetV(int i, double V);
	void SetTn(double T);
	void SetZone(int zone);
	void SetsourceWall(double *sourceWall);
	void divimp_Ti();
	void divimp_F();
	void divimp_E();
	void divimp_drift_E();
	void divimp_anomalous_diffusion();
	void SetX_new(int i, double X);
	void SetWeight(double w);
	void Pump_add(int i, double Sum);
	double Pump(int i);
	void Stat(int n);
	void Stat_Tri(int n);
	// void FluxCal_Grid();
	void FluxCal_Target();
	void FluxOutput();
	void Find(); // Find the index and Zone of [ X_new ];
	int Index(int i);
	double X(int i);
	double X_new(int i);
	double V(int i);
	double V_Charge(int i);
	double V_Grid(int i, int j, int k, int m);
	double V_D_1(int i, int j, int k, int m);
	double V_Grid_CX_Ion_Af(int i, int j, int k, int m);
	double V_Grid_CX_Ion_Be(int i, int j, int k, int m);
	double V_Grid_CX_Neu_Af(int i, int j, int k, int m);
	double V_Grid_CX_Neu_Be(int i, int j, int k, int m);
	double Vreal(int i);
	double sourceWall(int i);
	int XY(int i);
	int Tri_Index();
	// int XY_new(int i);
	int boundary_start();
	void SetXY(int i, int j);
	double Tn();
	int IfColl();
	int IfParticleOut(int i);
	int IfFlightOut();
	void SetIfFlightOut(int i);
	double lambda(int i, int j, int charge);
	double lambda_min(int charge);
	double lambda_now();
	void addCollProb(int zhonglei, ADAS *Coll1 = NULL, EIRENE *Coll2 = NULL);
	void CalTn();
	double n(int i, int j, int k);
	double T_n(int i, int j, int k);
	double Gamma();
	double Weight_Target(int i);
	double NumPar_Target(int i);
	void Dump_Flux();
	void Dump_array_2D(string type, string coll, int Charge);
	void Dump_2D(string type, string coll = "", int Charge = 0);
	void Dump_Tri(string type, string coll, int Charge = 0);
	void Clear(int n = 0);
	int NowZone();
	int Num_Target(int i);
	int Zone();
	// int Zone_new();
	int Charge();
	int ChargeTag();
	void SetChargeTag(int i);
	double Weight();
	void setname(string particle_name);
	void setPar(string particle_name, double mass, int Charge);
	void setmass(double m);
	void setfate(int i, int FATE, int sourcePar);
	void SetCharge(int i);
	double mass();
	string name();
	int MaxCharge();
	int fate(int i);
	int sourcePar(int i);
	double CollProb(int i, int j, int k);
	double CoreTotalcs(int i);
	void calIonRateADAS(ADAS *ParColl);
	void calRecRateADAS(ADAS *ParColl);
	void calCXRateADAS(ADAS *ParColl, int CX_cross = 0);

	void AddWallEro(int num_Ero_wall);
	void OutWallEro(int fate);
	void AddTargetEro(int num_Ero_wall);
	void OutTargetEro(int fate);
	void AddPlasmaBoundaryEro(int num_Ero_wall);
	void OutPlasmaBoundaryEro(int fate);

	double CalAngle(int num_wall);
	double CalAngle(double Sin, double Cos);

	void emit_particle(const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, int mode, Vector3 incidentDir = Vector3(0, 0, 0));
	/*void Out_n(string path);
	void Out_Sn_ion(string path);
	void Out_Sn_rec(string path);
	void Out_Sn_CX(string path);
	void Out_SE_ion(string path);
	void Out_SE_rec(string path);
	void Out_SE_CX(string path);
	void Out_Smu_ion(string path);
	void Out_Smu_rec(string path);
	void Out_Smu_CX(string path);
	void Out_Pra_ion(string path);
	void Out_Pra_rec(string path);
	void Out_Pra_CX(string path);
	void Out_Sn_MAR(string path);
	void Out_Sn_Diss(string path);
	void Out_Sn_Diss2(string path);*/
};

class PartoPar
{
private:
	bool ifstore_;
	bool ifChange_;
	queue<Particle *> Q_;
	queue<Particle> PQ_;
	queue<std::vector<double>> toV_;
	queue<int> toSource_;
	queue<int> toCollPar_;

public:
	PartoPar();
	bool ifstore();
	bool ifChange();
	void store(Particle *PP, double factor, Particle *QQ, std::vector<double> &V, int SOURCE, int toCollPar, int Charge);
	void ParChange();
};

void Push(Particle *PP);
void Push_Tri(Particle *PP);
void WallReflect();

// Particle
extern Particle H;
extern Particle H2;
extern Particle D;
extern Particle D2;
extern Particle T;
extern Particle T2;
extern Particle CD4;
extern Particle CD3;
extern Particle CD2;
extern Particle CD1;
extern Particle C;
extern Particle Ar;
extern Particle *P;

extern PartoPar PartoPar_Temp;
