#ifndef _ADAS_H_
#define _ADAS_H_
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;
int BinarySearch(double *array, int low, int high, double key);

class ADAS
{
private:
	string name_;
	int num_Z_ = 0;          // nuclear charge
	int num_Dens_ = 0;       // number of densities
	int num_Temp_ = 0;       // number of temperature
	int num_Zmin_ = 0;       // lowest ion charge + 1
	int num_Zmax_ = 0;       // highest ion charge
	int num_waveLength_ = 0; // number of the Photon Emissivity Coefficients
	double *waveLength_ = nullptr; // wavelength of the Photon Emissivity Coefficients
	double **Te_pec_ = nullptr;
	double **ne_pec_ = nullptr;
	double *Te_ = nullptr;
	double *ne_ = nullptr;
	double **cross_section_ = nullptr;
	void release();

public:
	ADAS();
	~ADAS();
	ADAS(const ADAS &) = delete;
	ADAS &operator=(const ADAS &) = delete;
	ADAS(string name);
	void read();
	void readpec(int num_skip, int K_unit);
	void readname(string name);
	void readbase(ifstream &in);
	void readne(ifstream &in);
	void readTe(ifstream &in);
	void readinfo(ifstream &in);

	string name();
	double waveLength(int i);
	double *Te(int i);
	double *ne(int i);
	double num_waveLength();
	double *cross_section(int i, int Charge = 0);
	double cal(double Te, double ne, int Charge = 0);
	double PecCal(double Te, double ne, int wave);
};
#endif