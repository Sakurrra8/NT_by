#ifndef _REFLECT_H_
#define _REFLECT_H_
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>

using namespace std;

int BinarySearch(const vector<double> &array, int low, int high, double key);
class Reflect
{
private:
	double n_a_[4];
	double n_b_[6];
	double E_a_[4];
	double E_b_[6];
	double epsilon_;
	double epsilon_L_;

	int num_ne_;
	int num_na_;
	vector<double> Rn_ne_;
	vector<double> Rn_na_;
	vector<vector<double>> Rn_;
	vector<double> RE_ne_;
	vector<double> RE_na_;
	vector<vector<double>> RE_;

public:
	void SetReflect(double n_a[4], double E_a[4], double n_b[6], double E_b[6], double epsilon_L);
	void ReadTrimData(int K_Reflect, int num_ne, int num_na, string file_Rn, string file_RE);
	double n_RefCoeff(int K_Reflect, double Energy_0, double angle = 30);
	double E_RefCoeff(int K_Reflect, double Energy_0, double angle = 30);
};

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

// 存储单个能级和角度网格点的 TRIM 数据块
struct TrimBlock
{
	double m1 = 0.0;					  // 入射粒子质量
	double z2 = 0.0;					  // 壁面原子序数
	double m2 = 0.0;					  // 壁面原子质量
	double E0 = 0.0;					  // 网格点入射能 (eV)
	double theta0 = 0.0;				  // 网格点入射角 (度)
	double RN = 0.0;					  // 粒子反射系数
	std::vector<double> energy_quantiles; // 5个能量分布分位数
};

struct DWReflectionSample
{
	double energy_eV = 0.0;
	double cos_polar = 1.0;
	double cos_azimuth = 1.0;
};

class DWTrimReflection
{
private:
	struct Block
	{
		double incident_energy_eV = 0.0;
		double incident_angle_deg = 0.0;
		double particle_reflection = 0.0;
		double energy_reflection = 0.0;
		std::array<double, 7> energy{};
		std::array<std::array<double, 7>, 5> cos_polar{};
		std::array<std::array<std::array<double, 7>, 5>, 5> cos_azimuth{};
	};

	std::array<double, 12> energies_{};
	std::array<double, 7> angles_{};
	std::array<std::array<Block, 7>, 12> blocks_{};
	bool loaded_ = false;

	static double SampleQuantile(const std::array<double, 7> &values, double xi);
	static int FindInterval(const double *grid, int size, double value);
	const Block &BlockAt(int energy_index, int angle_index) const;

public:
	bool Load(const std::string &filename);
	bool IsLoaded() const;
	double ReflectionProbability(double incident_energy_eV, double incident_angle_deg) const;
	double MeanReflectedEnergy(double incident_energy_eV, double incident_angle_deg) const;
	DWReflectionSample Sample(double incident_energy_eV, double incident_angle_deg,
							  double xi_energy, double xi_polar, double xi_azimuth,
							  double minimum_energy_eV = 0.0) const;
};

// EIRENE ILREF=2 default surface model for D incident on Fe.
class EireneDFeReflection
{
public:
	static constexpr double ThermalTemperatureEV = 0.0388;

	static double ReflectionProbability(double incident_energy_eV,
										double incident_angle_deg);
	static double MeanReflectedEnergy(double incident_energy_eV,
								  double incident_angle_deg);
	static double SampleReflectedEnergy(double incident_energy_eV,
									double incident_angle_deg,
									double xi);
};

#endif
