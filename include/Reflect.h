#ifndef _REFLECT_H_
#define _REFLECT_H_
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

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

#endif