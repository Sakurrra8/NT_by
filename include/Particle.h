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

void Edge(int Start[2], int End[2], double a[], double b[], int *k);
bool IfinEdge(double a[], double b[], int k, double x, double y);
int Xfind(int Start[], int End[], double x, double y);
int Yfind(int Start[], int End[], double x, double y);
int ZonefromXY(int X, int Y);
string fatename(int i);
string sourcename(int i);

enum class EireneDensitySource
{
	None = 0,
	Electron,
	HIon,
	DIon,
	TIon,
	HNeutralTri,
	H2NeutralTri,
	DNeutralTri,
	D2NeutralTri,
	TNeutralTri,
	T2NeutralTri
};

enum class EireneArgumentMode
{
	ElectronDensityTemperature = 0,
	SameDensityTemperature
};

enum class SourceStratum
{
	Unknown = 0,
	IT,
	OT,
	MCW,
	Recombination,
	Puff,
	Methane,
	Carbon,
	Argon,
	PlasmaBoundary,
	Core,
	Count
};

string sourceStratumName(SourceStratum source);

class ParCollCar
{
private:
	int K_;
	int e_or_i_;
	int cc;

	double Sn_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Smu_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Smu1_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Smu2_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double SE_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_n_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_mu_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_mu1_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_mu2_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_E_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Pra_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double crossSection[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double cs_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double SE_n_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double Num_E_n_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	// The launch summary is not Sn_by_stratum. A future source-resolved
	// diagnostic must add B2 and triangle counters keyed by SourceStratum
	// without replacing the existing aggregate Sn_ / Tri_Sn_ values. It must
	// separately cover Ion_D, Diss2_D2, DS1_D2p, DS2_D2p, DS3_D2p, Ion_D2,
	// and CX_D2 before comparing with EIRENE IT/OT/MCW source terms.

	double Cor_cs_;
	double cs_now_;
	bool defer_stats_{false};
	double deferred_stat_scale_{1.0};
	vector<double> pendingNumN_, pendingNumMu_, pendingNumMu1_, pendingNumMu2_, pendingNumE_, pendingNumEN_;
	vector<double> pendingTriNumN_, pendingTriNumMu_, pendingTriNumMu0_, pendingTriNumMu1_, pendingTriNumMu2_, pendingTriNumE_, pendingTriNumEN_;
	std::size_t statGridIndex(int i, int j) const;
	void clearDeferredStats();
	vector<double> V_relative_;
	double V_2_relative_;
	double V_relative_2_;
	int num_trimesh_;
	EIRENE *eirene_rate_{nullptr};
	EireneDensitySource eirene_density_source_{EireneDensitySource::None};
	EireneArgumentMode eirene_argument_mode_{EireneArgumentMode::ElectronDensityTemperature};
	double eirene_scale_{1.0};
	double eireneDensity(int i, int j) const;
	double eireneDensityTri(int tri) const;
	double eireneTemperatureForDensityTri(int tri) const;
	double eireneRate(int i, int j) const;
	double eireneRateTri(int tri) const;

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
	void BeginDeferredStats(double scale);
	void EndDeferredStats();
	int K();
	void Set_K(int i);
	void Stat(int Charge, double *n[N_POLOIDAL_GRID][N_RADIAL_GRID], ADAS *PraADAS, const std::vector<std::vector<double>> &ni = {});
	void Stat1(int Charge, double *n[N_POLOIDAL_GRID][N_RADIAL_GRID], ADAS *PraADAS, const std::vector<std::vector<double>> &ni = {});
	void StatEIRENE(int Charge, double *n[N_POLOIDAL_GRID][N_RADIAL_GRID], EIRENE *PraEIRENE, double E, double **ni = NULL);
	void RecStat(int Charge, double *n[N_POLOIDAL_GRID][N_RADIAL_GRID], ADAS *PraADAS);

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
	double Sn(int XY);
	double Smu(int XY);
	double SE(int XY);
	void Setcs(int i, double cs);
	void Setcs(int i, int j, double cs);
	double cs(int XY[3], double Zone);
	double cs(int i, int j);
	double cs(int i);
	void SetCor_cs(double Cor_cs);
	double Cor_cs();
	double cs_now();
	double EireneRate(double density, double arg1, double arg2) const;

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
	void EIRENEInput(EIRENE *ParColl,
					 EireneDensitySource density_source = EireneDensitySource::Electron,
					 EireneArgumentMode argument_mode = EireneArgumentMode::ElectronDensityTemperature,
					 double scale = 1.0);

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
	int NREG_{0}, NBIN_{0};
	double E_MIN_{0.0}, E_MAX_{0.0};
	bool trackMom_{false}, trackPow_{false};
	int NX_{0}, NY_{0};
	bool gridEnabled_{false};

	std::vector<double> gridHist_;				 // NX x NY x NBIN
	std::vector<std::array<double, 3>> gridMom_; // NX x NY x 3
	std::vector<double> gridPow_;				 // NX x NY
	std::vector<double> flux_;					 // NREG x NREG
	std::vector<double> hist_;					 // NREG x NREG x NBIN
	std::vector<std::array<double, 3>> momFlux_; // NREG x NREG x 3
	std::vector<double> powFlux_;				 // NREG x NREG
	bool inited_{false};

	std::size_t regionIndex(int from, int to) const;
	std::size_t histIndex(int from, int to, int bin) const;
	std::size_t gridIndex(int i, int j) const;
	std::size_t gridHistIndex(int i, int j, int bin) const;
};

class Particle
{
private:
	/// @brief fixed property of this particle
	string name_;
	double mass_;
	int MaxCharge_;
	bool owns_storage_{false};
	void allocateStorage();
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
	double X_old_[3];		// position at last time step
	double dt_;				// the dt calculated from lambda_min
	double dt_trace_;		// the dt calculated from trace
	int Zone_;				// Particle Zone in Tokamak
	double Tn_;				// Particle Temperture
	double Weight_;			// Particle weight
	int split_depth_{0};	// Number of variance-reduction split generations
	int importance_region_{-1}; // Last independently defined physical importance region
	SourceStratum source_stratum_{SourceStratum::Unknown};
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
	double *lambda_[N_POLOIDAL_GRID][N_RADIAL_GRID];

	/// @brief Initial launch parameters
	double NumPar_Target_;		 // for Particle from Recycling
	double Weight_Target_;		 // for Particle from Recycling
	double NumPar_sum_Target_;	 // for Particle from Recycling
	double NumPar_Grid_[N_POLOIDAL_GRID][N_RADIAL_GRID]; // for Particle from Grid Collision
	double Weight_Grid_[N_POLOIDAL_GRID][N_RADIAL_GRID]; // for Particle from Grid Collision
	double NumPar_sum_Grid_;	 // for Particle from Grid Collision

	std::vector<double> Recycled_counts_;
	std::vector<double> Recombin_counts_;
	std::vector<double> Recycled_cdf_;
	std::vector<double> Recombin_cdf_;
	double Recycled_source_sum_{0.};
	double Recombin_source_sum_{0.};
	std::vector<double> T_Init_;

	/// @brief Contiguous storage backing the pointer views below
	std::vector<double> nStorage_, TStorage_, EStorage_;
	std::vector<double> sumNStorage_, sumEStorage_, lambdaStorage_;
	std::vector<double> numVD1Storage_, lambdaMinStorage_;
	std::vector<double> vGridStorage_, vD1Storage_, sumVStorage_, sumVD1Storage_;
	std::vector<double> vCxIonBeforeStorage_, vCxIonAfterStorage_;
	std::vector<double> vCxNeutralBeforeStorage_, vCxNeutralAfterStorage_;
	std::vector<double> sumVCxIonBeforeStorage_, sumVCxIonAfterStorage_;
	std::vector<double> sumVCxNeutralBeforeStorage_, sumVCxNeutralAfterStorage_;
	bool defer_flight_stats_{false};
	double deferred_flight_stat_scale_{1.0};
	std::vector<double> pendingGridN_, pendingGridE_, pendingGridV_;
	std::vector<double> pendingTriN_, pendingTriE_, pendingTriV_;
	std::vector<double> pendingCxIonBefore_, pendingCxIonAfter_;
	std::vector<double> pendingCxNeutralBefore_, pendingCxNeutralAfter_;
	std::vector<std::array<double, static_cast<std::size_t>(SourceStratum::Count)>> launchedWeightByStratum_;
	std::vector<std::array<unsigned long long, static_cast<std::size_t>(SourceStratum::Count)>> launchedEventsByStratum_;
	std::size_t gridScalarIndex(int i, int j, int charge) const;
	std::size_t gridVectorIndex(int i, int j, int component, int charge) const;
	std::size_t triScalarIndex(int tri, int charge) const;
	std::size_t triVectorIndex(int tri, int component, int charge) const;
	void addDensityStatGrid(int i, int j, int charge, double n);
	void addFlightStatGrid(int i, int j, int charge, double n, double energy, const std::vector<double> &v);
	void addFlightStatTri(int tri, int charge, double n, double energy, const std::vector<double> &v);
	void addCxVelocityBeforeGrid(int i, int j, int component, int charge, double ion, double neutral, double n);
	void addCxVelocityAfterGrid(int i, int j, int component, int charge, double ion, double neutral, double n);
	double collisionStatWeight() const;
	double diagnosticEventWeight() const;
	void recordSourceLaunch();
	void beginDeferredCollisionStats(double scale);
	void endDeferredCollisionStats();
	void SampleIonVelocity(int isotope);

	/// @brief Statistical pointer views
	double *n_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *T_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *E_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *V_Grid_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *V_D_1_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Num_V_D_1_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *Sum_n_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *Sum_E_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	double *Sum_V_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Sum_V_D_1_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];

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
	std::vector<double> Tri_D2p_track_time_;
	std::vector<std::array<double, 3>> Tri_D2p_DS_weight_;
	std::vector<double> Tri_D2p_boundary_loss_weight_;
	unsigned long long D2p_created_by_ion_{0};
	unsigned long long D2p_created_by_cx_{0};
	unsigned long long D2p_track_steps_{0};
	unsigned long long D2p_DS_events_[3]{0, 0, 0};
	unsigned long long D2p_secondary_D_events_[2]{0, 0};
	unsigned long long D2p_boundary_loss_{0};
	unsigned long long D2p_max_steps_loss_{0};
	unsigned long long D2p_current_flight_steps_{0};
	double D2p_created_by_ion_weight_{0.};
	double D2p_created_by_cx_weight_{0.};
	double D2p_DS_weight_[3]{0., 0., 0.};
	double D2p_secondary_D_weight_[2]{0., 0.};
	double D2p_boundary_loss_weight_{0.};
	double D2p_max_steps_loss_weight_{0.};
	double D2p_sum_weight_segment_dt_{0.};
	double D2p_sum_weight_segment_length_{0.};
	double D2p_sum_weight_charge_speed_{0.};
	double D2p_sum_weight_neutral_speed_{0.};
	double D2p_sum_weight_segment_speed_{0.};
	unsigned long long D2p_fallback_to_neutral_velocity_count_{0};
	double D2p_fallback_to_neutral_velocity_weight_{0.};
	double D2p_sum_segment_dt_{0.};
	double D2p_sum_segment_length_{0.};
	double D2p_min_charge_speed_{0.};
	double D2p_max_charge_speed_{0.};
	unsigned long long D2p_initial_velocity_events_[2]{0, 0};
	double D2p_initial_velocity_weight_[2]{0., 0.};
	double D2p_initial_sum_weight_neutral_speed_[2]{0., 0.};
	double D2p_initial_sum_weight_charge_speed_[2]{0., 0.};
	double D2p_initial_sum_weight_abs_v_parallel_[2]{0., 0.};
	double D2p_initial_sum_weight_v_perp_[2]{0., 0.};
	double D2p_initial_sum_weight_neutral_energy_[2]{0., 0.};
	double D2p_initial_sum_weight_charge_energy_[2]{0., 0.};
	double D2p_initial_min_charge_speed_[2]{0., 0.};
	double D2p_initial_max_charge_speed_[2]{0., 0.};
	double D2p_sum_weight_length_over_neutral_speed_{0.};
	double D2p_sum_weight_length_over_charge_speed_{0.};
	double D2p_sum_weight_length_over_ion_thermal_speed_{0.};
	double D2p_CX_sum_weight_length_over_neutral_speed_{0.};
	double D2p_CX_sum_weight_length_over_ion_thermal_speed_{0.};
	unsigned long long D2p_low_charge_speed_count_[4]{0, 0, 0, 0};
	double D2p_low_charge_speed_weight_dt_[4]{0., 0., 0., 0.};
	bool D2p_current_created_by_cx_{false};
	static constexpr unsigned long long MaxD2pFlightSteps = 100000;
	static constexpr unsigned long long MaxNeutralTriSteps = 1000;

	double *V_Grid_CX_Ion_Be_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *V_Grid_CX_Ion_Af_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *V_Grid_CX_Neu_Be_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *V_Grid_CX_Neu_Af_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Sum_V_CX_Ion_Be_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Sum_V_CX_Ion_Af_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Sum_V_CX_Neu_Be_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];
	double *Sum_V_CX_Neu_Af_[N_POLOIDAL_GRID][N_RADIAL_GRID][4];

	vector<array<double, 51>> Sum_n_array_[N_POLOIDAL_GRID][N_RADIAL_GRID];
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
	vector<double> CollProb_[N_POLOIDAL_GRID][N_RADIAL_GRID];
	vector<double> Totalcs_[N_POLOIDAL_GRID][N_RADIAL_GRID];

	vector<double> n_Flux_Grid_; // [N_POLOIDAL_GRID][N_RADIAL_GRID][4], row-major
	vector<double> T_Flux_Grid_; // [N_POLOIDAL_GRID][N_RADIAL_GRID][4], row-major
	std::size_t fluxGridIndex(int i, int j, int direction) const;
	FluxTracker FT_;

public:
	struct State
	{
		int Tri_Index;
		int Charge;
		int ChargeTag;
		int Zone;
		int XY[3];
		int GridIndex;
		int boundary_start;
		int IfColl;
		int IfFlightOut;
		int fate[2];
		int sourcePar[2];
		int sourceGrid[2];
		int Zone_old;
		int GridIndex_old;
		double X[3];
		double X_new[3];
		double X_old[3];
		double Xold[3];
		double V_Charge[8];
		double dt;
		double dt_trace;
		double Tn;
		double Weight;
		int splitDepth;
		int importanceRegion;
		SourceStratum sourceStratum;
		double sourceWall[4];
		double lambda_now;
		double d_flight;
		double Rand_flight;
		unsigned long long D2p_current_flight_steps;
		bool D2p_current_created_by_cx;
		std::vector<double> V;
		double intersection[10][6];
	};

	Particle();
	Particle(string particle_name, double mass, int MaxCharge, int Index);
	Particle(const Particle &) = delete;
	Particle &operator=(const Particle &) = delete;
	~Particle();
	void ParInit(vector<int> &K_collision, const std::vector<std::vector<int>> &Tri_B2);

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
	vector<WallEro> AlltoTarget_;
	vector<WallEro> CXtoTarget_;
	vector<WallEro> RectoTarget_;
	vector<WallEro> AlltoPlasmaBoundary_;
	vector<WallEro> CXtoPlasmaBoundary_;
	vector<WallEro> RectoPlasmaBoundary_;

	/*double MAR_cs_[N_POLOIDAL_GRID][N_RADIAL_GRID], MAR_Cor_cs_;	 // for D2
	double Diss1_cs_[N_POLOIDAL_GRID][N_RADIAL_GRID], Diss1_Cor_cs_;	 // for D2
	double Diss2_cs_[N_POLOIDAL_GRID][N_RADIAL_GRID], Diss2_Cor_cs_; // for D2
	double Diss3_cs_[N_POLOIDAL_GRID][N_RADIAL_GRID], Diss3_Cor_cs_; // for D2
	double DS_cs_[N_POLOIDAL_GRID][N_RADIAL_GRID][4], DS_Cor_cs_[4];*/

	void Particlefrom(Particle *A, double K = 1, int Charge = -100);
	int sampleTargetPlate(const std::vector<double> &Counts);
	int sampleRecycledTarget();
	double RecycledSourceSum();
	double RecombinSourceSum();
	int sampleRecombinCell();
	void RecycledCal(std::vector<double> &NumPar_wall);
	void RecycledCal(std::vector<double> &NumPar_wall, std::vector<double> &T_Init);
	void RecombinCal(std::vector<double> &NumPar_wall, std::vector<double> &T_Init);
	void CalWeight1(int num);
	void CalWeight2(std::vector<double> &NumPar, int num);
	void Init(int k, int z = 0);
	void VtoVcharge();
	void VchargetoV();
	void Vchargefix(); // when the charged particle flight to next grid, V_charge[0,1,2] should be fixed
	void advanceChargedMinimalTrace(double requested_dt);
	bool isHydrogenMoleculeIon() const;
	void markD2pJustCreated(bool created_by_cx);
	void track();
	void CalLambda();
	void Caltrace();
	void Caltrace_Tri();
	void NumParStat();
	void BeginDeferredFlightStats(double scale);
	void EndDeferredFlightStats();
	void ApplyRussianRoulette();
	State SaveState() const;
	void RestoreState(const State &state);
	void ApplyRegionalImportance(std::queue<State> &pending_states);
	int PhysicalImportanceRegion() const;
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
	void divimp_Ti(double step_dt);
	void divimp_F(double step_dt);
	void divimp_E(double step_dt);
	void divimp_drift_E(double step_dt);
	void divimp_anomalous_diffusion();
	void SetX_new(int i, double X);
	void SetWeight(double w);
	void Pump_add(int i, double Sum);
	double Pump(int i);
	void Stat(int n);
	void Stat_Tri(int n);
	void DumpD2pBalance_B2();
	void DumpD2pBalance_Tri();
	void DumpD2pPhysicsDecomposition_B2();
	void DumpD2pTrackLengthTri();
	void UseD2pTransportDensityForOutput();
	void AppendSourceStratumSummary(std::ostream &out) const;
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
	double V_Grid_Tri(int i, int k, int m);
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
	double Tn_Tri(int i, int k);
	int IfColl();
	int IfParticleOut(int i);
	int IfFlightOut();
	void SetIfFlightOut(int i);
	double lambda(int i, int j, int charge);
	double lambda_min(int charge);
	double lambda_now();
	void CalTn();
	double n(int i, int j, int k);
	double n_Tri(int i, int k);
	double T_n(int i, int j, int k);
	double Gamma();
	double Weight_Target();
	double NumPar_Target();
	double NumPar_Grid();
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
	SourceStratum sourceStratum() const;
	void setSourceStratum(SourceStratum source);
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
	void NeutralFluxStatistics(int k, int oritation);
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
