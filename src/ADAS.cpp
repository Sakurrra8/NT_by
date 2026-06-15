#include "ADAS.h"

ADAS::ADAS() {}
ADAS::ADAS(string name)
{
	name_ = name;
	read();
}

ADAS::~ADAS()
{
	release();
}

void ADAS::release()
{
	if (cross_section_ != nullptr)
	{
		const int rows = (num_waveLength_ > 0) ? num_waveLength_ : num_Zmax_;
		for (int i = 0; i < rows; ++i)
			free(cross_section_[i]);
		free(cross_section_);
		cross_section_ = nullptr;
	}
	if (ne_pec_ != nullptr)
	{
		for (int i = 0; i < num_waveLength_; ++i)
			free(ne_pec_[i]);
		free(ne_pec_);
		ne_pec_ = nullptr;
	}
	if (Te_pec_ != nullptr)
	{
		for (int i = 0; i < num_waveLength_; ++i)
			free(Te_pec_[i]);
		free(Te_pec_);
		Te_pec_ = nullptr;
	}
	free(waveLength_);
	free(Te_);
	free(ne_);
	waveLength_ = nullptr;
	Te_ = nullptr;
	ne_ = nullptr;
	num_Z_ = 0;
	num_Dens_ = 0;
	num_Temp_ = 0;
	num_Zmin_ = 0;
	num_Zmax_ = 0;
	num_waveLength_ = 0;
}

void ADAS::readname(string name)
{
	ADAS::name_ = name;
}

void ADAS::readbase(ifstream &in)
{
	in >> ADAS::num_Z_;	   // nuclear charge
	in >> ADAS::num_Dens_; // number of densities
	in >> ADAS::num_Temp_; // number of temperature
	in >> ADAS::num_Zmin_; // lowest ion charge + 1
	in >> ADAS::num_Zmax_; // highest ion charge
}

void ADAS::readne(ifstream &in)
{
	for (int i = 0; i < ADAS::num_Dens_; i++)
	{
		in >> ADAS::ne_[i];
		ADAS::ne_[i] += 6;
	}
}

void ADAS::readTe(ifstream &in)
{
	for (int i = 0; i < ADAS::num_Temp_; i++)
	{
		in >> ADAS::Te_[i];
	}
}

void ADAS::readinfo(ifstream &in)
{
	string line;
	for (int k = 0; k < num_Zmax_; k++)
	{
		std::getline(in, line);
		std::getline(in, line);
		for (int i = 0; i < num_Temp_; i++)
			for (int j = 0; j < num_Dens_; j++)
			{
				in >> cross_section_[k][i * num_Dens_ + j];
				cross_section_[k][i * num_Dens_ + j] -= 6;
			}
	}
}

string ADAS::name()
{
	return name_;
}

double *ADAS::Te(int i)
{
	return &ADAS::Te_[i];
}

double *ADAS::ne(int i)
{
	return &ADAS::ne_[i];
}

double ADAS::num_waveLength()
{
	return num_waveLength_;
}

/// @brief Cross Section
/// @param i for i = num_Dens_ * b + a
/// @param Charge Ionize to which Charge, or Recombination from which Charge
/// @return
double *ADAS::cross_section(int i, int Charge)
{
	return &ADAS::cross_section_[Charge - 1][i];
}

double ADAS::cal(double Te, double ne, int Charge)
{
	if (num_Dens_ < 2 || num_Temp_ < 2 || ne_ == nullptr || Te_ == nullptr ||
		cross_section_ == nullptr || Charge < 0 || Charge >= num_Zmax_)
	{
		std::cerr << "ADAS.cal() called with unavailable data: " << name_ << std::endl;
		return 0.0;
	}
	if (ne < 0)
	{
		std::cerr << "T in ADAS.cal() less than 0" << endl;
		ne = 0.1;
	}
	if (Te < 0)
	{
		std::cerr << "T in ADAS.cal() less than 0" << endl;
		Te = 0.1;
	}
	Te = log10(Te);
	ne = log10(ne);
	int a, b, c, d, k;
	double x2x1, y2y1, x2x, y2y, yy1, xx1;

	a = BinarySearch(ne_, 0, num_Dens_ - 1, ne);
	b = BinarySearch(Te_, 0, num_Temp_ - 1, Te);

	if (a == -1 || b == -1)
	{
		if (ne < ne_[0])
		{
			a = 0;
			ne = ne_[0];
		}
		else if (ne > ne_[num_Dens_ - 1])
		{
			a = num_Dens_ - 2;
			ne = ne_[num_Dens_ - 1];
		}
		else if (a == -1)
			cerr << "BinartSearch in ADAS has some problem." << endl;
		if (Te < Te_[0])
		{
			b = 0;
			Te = Te_[0];
		}
		else if (Te > Te_[num_Temp_ - 1])
		{
			b = num_Temp_ - 2;
			Te = Te_[num_Temp_ - 1];
		}
		else if (b == -1)
			cerr << "BinartSearch in ADAS has some problem." << endl;
	}

	x2x1 = ne_[a + 1] - ne_[a];
	y2y1 = Te_[b + 1] - Te_[b];
	x2x = ne_[a + 1] - ne;
	y2y = Te_[b + 1] - Te;
	xx1 = ne - ne_[a];
	yy1 = Te - Te_[b];

	return pow(10., 1.0 / (x2x1 * y2y1) * (cross_section_[Charge][num_Dens_ * b + a] * x2x * y2y + cross_section_[Charge][num_Dens_ * b + (a + 1)] * xx1 * y2y + cross_section_[Charge][num_Dens_ * (b + 1) + a] * x2x * yy1 + cross_section_[Charge][num_Dens_ * (b + 1) + (a + 1)] * xx1 * yy1));
}

double ADAS::waveLength(int i)
{
	return waveLength_[i];
}

int BinarySearch(double *array, int low, int high, double key)
{
	if (array == nullptr || low < 0 || high <= low)
		return -1;
	if (key < array[low] || key > array[high])
		return -1;
	if (key == array[high])
		return high - 1;

	while (high - low > 1)
	{
		const int mid = low + (high - low) / 2;
		if (key < array[mid])
			high = mid;
		else
			low = mid;
	}
	return low;
}

void ADAS::read()
{
	// int resolved;
	// double x, y;
		release();
	string path_ADAS, input, line;
	// resolved = 0;
	// x = y = 0;

	/// path
	input = "Inputfile/database/";
	path_ADAS = input + name_;
	std::ifstream fp;
	fp.open(path_ADAS, std::ios::in);
	// tell you if
	if (!fp.is_open())
	{
		std::cerr << "Cannot open ADAS file: " << path_ADAS << std::endl;
		return;
	}
	readbase(fp);
	if (!fp || num_Dens_ < 2 || num_Temp_ < 2 || num_Zmax_ <= 0)
	{
		std::cerr << "Invalid ADAS dimensions in: " << path_ADAS << std::endl;
		release();
		return;
	}
	ne_ = (double *)malloc(sizeof(double) * num_Dens_);
	Te_ = (double *)malloc(sizeof(double) * num_Temp_);
	cross_section_ = (double **)malloc(sizeof(double *) * (num_Zmax_));
	for (int i = 0; i < num_Zmax_; i++)
		cross_section_[i] = (double *)malloc(sizeof(double) * (num_Dens_ * num_Temp_));
	std::getline(fp, line);
	std::getline(fp, line);
	readne(fp);
	readTe(fp);
	readinfo(fp);
	fp.close();
}

/// @brief read the PEC data from ADAS
/// @param num_skip The number of rows before the data
/// @param K_unit Whether there are units after the wavelength
void ADAS::readpec(int num_skip, int K_unit)
{
	// int resolved;
	// double x, y;
		release();
	string path_ADAS, input, line;
	// resolved = 0;
	// x = y = 0;

	/// path
	input = "Inputfile/database/";
	path_ADAS = input + name_;
	std::ifstream fp;
	fp.open(path_ADAS, std::ios::in);
	// tell you if
	if (!fp.is_open())
	{
		std::cerr << "Cannot open ADAS PEC file: " << path_ADAS << std::endl;
		return;
	}
	fp >> num_waveLength_; // nuclear charge
	if (!fp || num_waveLength_ <= 0)
	{
		std::cerr << "Invalid ADAS PEC header in: " << path_ADAS << std::endl;
		release();
		return;
	}
	for (int i = 0; i < num_skip; i++)
		std::getline(fp, line);
	waveLength_ = (double *)malloc(sizeof(double) * num_waveLength_);
	ne_pec_ = (double **)malloc(sizeof(double *) * num_waveLength_);
	Te_pec_ = (double **)malloc(sizeof(double *) * num_waveLength_);
	cross_section_ = (double **)malloc(sizeof(double *) * (num_waveLength_));
	for (int i = 0; i < num_waveLength_; i++)
	{
		fp >> waveLength_[i];
		// std::cout << waveLength_[i] << endl;
		if (K_unit)
			fp >> line;
		// std::cout << line << endl;
		fp >> num_Dens_; // number of densities
		// std::cout << num_Dens_ << endl;
		fp >> num_Temp_; // number of temperature
		// std::cout << num_Temp_ << endl;
		std::getline(fp, line);

		ne_pec_[i] = (double *)malloc(sizeof(double) * num_Dens_);
		Te_pec_[i] = (double *)malloc(sizeof(double) * num_Temp_);
		cross_section_[i] = (double *)malloc(sizeof(double) * (num_Dens_ * num_Temp_));
		for (int j = 0; j < num_Dens_; j++)
		{
			fp >> ne_pec_[i][j];
			// std::cout << ne_pec_[i][j] << '\t';
			ne_pec_[i][j] = log10(ne_pec_[i][j] * 1e6);
		}
		for (int j = 0; j < num_Temp_; j++)
		{
			fp >> Te_pec_[i][j];
			// std::cout << Te_pec_[i][j] << '\t';
			Te_pec_[i][j] = log10(Te_pec_[i][j]);
		}
		for (int j = 0; j < num_Dens_ * num_Temp_; j++)
		{
			fp >> cross_section_[i][j];
			// std::cout << cross_section_[i][j] << '\t';
			cross_section_[i][j] = log10(cross_section_[i][j]);
		}
		// if (i == 10)
		// for (int k = 0; k < 7; k++)
		// std::cout << pow(10, cross_section_[i][3 + k * 14]) << endl;
	}
	fp.close();
}

double ADAS::PecCal(double Te, double ne, int wave)
{
	if (num_Dens_ < 2 || num_Temp_ < 2 || wave < 0 || wave >= num_waveLength_ ||
		ne_pec_ == nullptr || Te_pec_ == nullptr || cross_section_ == nullptr ||
		ne_pec_[wave] == nullptr || Te_pec_[wave] == nullptr || cross_section_[wave] == nullptr)
	{
		std::cerr << "ADAS.PecCal() called with unavailable data: " << name_ << std::endl;
		return 0.0;
	}
	if (ne < 0)
	{
		std::cerr << "T in ADAS.cal() less than 0" << endl;
		ne = 0.1;
	}
	if (Te < 0)
	{
		std::cerr << "T in ADAS.cal() less than 0" << endl;
		Te = 0.1;
	}
	Te = log10(Te);
	ne = log10(ne);
	int a, b, c, d, k;
	double x2x1, y2y1, x2x, y2y, yy1, xx1;

	a = BinarySearch(ne_pec_[wave], 0, num_Dens_ - 1, ne);
	b = BinarySearch(Te_pec_[wave], 0, num_Temp_ - 1, Te);
	// std::cout << a << '\t' << b << endl;
	if (a == -1 || b == -1)
	{
		if (ne < ne_pec_[wave][0])
		{
			a = 0;
			ne = ne_pec_[wave][0];
		}
		else if (ne > ne_pec_[wave][num_Dens_ - 1])
		{
			a = num_Dens_ - 2;
			ne = ne_pec_[wave][num_Dens_ - 1];
		}
		else if (a == -1)
			cerr << "BinartSearch in ADAS has some problem." << endl;
		if (Te < Te_pec_[wave][0])
		{
			b = 0;
			Te = Te_pec_[wave][0];
		}
		else if (Te > Te_pec_[wave][num_Temp_ - 1])
		{
			b = num_Temp_ - 2;
			Te = Te_pec_[wave][num_Temp_ - 1];
		}
		else if (b == -1)
			cerr << "BinartSearch in ADAS has some problem." << endl;
	}

	x2x1 = ne_pec_[wave][a + 1] - ne_pec_[wave][a];
	y2y1 = Te_pec_[wave][b + 1] - Te_pec_[wave][b];
	x2x = ne_pec_[wave][a + 1] - ne;
	y2y = Te_pec_[wave][b + 1] - Te;
	xx1 = ne - ne_pec_[wave][a];
	yy1 = Te - Te_pec_[wave][b];

	return pow(10., 1.0 / (x2x1 * y2y1) * (cross_section_[wave][num_Temp_ * a + b] * x2x * y2y + cross_section_[wave][num_Temp_ * a + (b + 1)] * x2x * yy1 + cross_section_[wave][num_Temp_ * (a + 1) + b] * xx1 * y2y + cross_section_[wave][num_Temp_ * (a + 1) + (b + 1)] * xx1 * yy1));
}