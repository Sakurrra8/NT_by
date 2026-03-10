#ifndef _REFLECT_H_
#define _REFLECT_H_
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int BinarySearch(vector<double> &array, int low, int high, double key);
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
#endif