#include "Reflect.h"
#include <stdexcept>

void Reflect::SetReflect(double n_a[4], double E_a[4], double n_b[6], double E_b[6], double epsilon_L)
{
	epsilon_L_ = epsilon_L;
	for (int i = 0; i < 4; i++)
	{
		n_a_[i] = n_a[i];
		E_a_[i] = E_a[i];
	}
	for (int i = 0; i < 6; i++)
	{
		n_b_[i] = n_b[i];
		E_b_[i] = E_b[i];
	}
}

void Reflect::ReadTrimData(int K_Reflect, int num_ne, int num_na, string file_Rn, string file_RE)
{
	if (K_Reflect == 2)
	{
		num_ne_ = num_ne;
		num_na_ = num_na;
		Rn_na_.resize(num_na);
		Rn_ne_.resize(num_ne);
		Rn_.resize(num_ne);
		for (int i = 0; i < num_ne; i++)
		{
			Rn_[i].resize(num_na);
		}
		RE_na_.resize(num_na);
		RE_ne_.resize(num_ne);
		RE_.resize(num_ne);
		for (int i = 0; i < num_ne; i++)
		{
			RE_[i].resize(num_na);
		}
		string line;
		std::ifstream fp;
		fp.open(file_Rn, std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for Rn_ have some problem!!!\n";
		}
		for (int i = 0; i < num_na; i++)
		{
			fp >> Rn_na_[i];
			// std::cout << Rn_na_[i] << endl;
		}
		for (int i = 0; i < num_ne; i++)
		{
			fp >> Rn_ne_[i];
			Rn_ne_[i] = log(Rn_ne_[i]);
			// std::cout << Rn_ne_[i] << endl;
			for (int j = 0; j < num_na; j++)
			{
				fp >> Rn_[i][j];
				Rn_[i][j] = log(Rn_[i][j]);
				// std::cout << Rn_[i][j] << '\t';
			}
			// std::cout << endl;
		}
		fp.close();

		fp.open(file_RE, std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for RE_ have some problem!!!\n";
		}
		for (int i = 0; i < num_na; i++)
		{
			fp >> RE_na_[i];
		}
		for (int i = 0; i < num_ne; i++)
		{
			fp >> RE_ne_[i];
			RE_ne_[i] = log(RE_ne_[i]);
			for (int j = 0; j < num_na; j++)
			{
				fp >> RE_[i][j];
				RE_[i][j] = log(RE_[i][j]);
			}
		}
		fp.close();
	}
	if (K_Reflect == 3)
	{
		num_ne_ = num_ne;
		num_na_ = num_na;
		Rn_na_.resize(num_na);
		Rn_ne_.resize(num_ne);
		Rn_.resize(num_ne);
		for (int i = 0; i < num_ne; i++)
		{
			Rn_[i].resize(num_na);
		}
		RE_na_.resize(num_na);
		RE_ne_.resize(num_ne);
		RE_.resize(num_ne);
		for (int i = 0; i < num_ne; i++)
		{
			RE_[i].resize(num_na);
		}

		string line;
		std::ifstream fp;
		fp.open("Inputfile/database/E_input.txt", std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for E_input have some problem!!!\n";
		}
		for (int i = 0; i < num_ne; i++)
		{
			fp >> Rn_ne_[i];
			Rn_ne_[i] = log(Rn_ne_[i]);
			RE_ne_[i] = Rn_ne_[i];
		}
		fp.close();
		fp.open("Inputfile/database/angle.txt", std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for angle have some problem!!!\n";
		}
		for (int i = 0; i < num_na; i++)
		{
			fp >> Rn_na_[i];
			RE_na_[i] = Rn_na_[i];
		}
		fp.close();
		fp.open(file_Rn, std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for Rn_ have some problem!!!\n";
		}
		for (int i = 0; i < num_ne; i++)
		{
			for (int j = 0; j < num_na; j++)
			{
				fp >> Rn_[i][j];
				Rn_[i][j] = log(Rn_[i][j]);
				// std::cout << Rn_[i][j] << '\t';
			}
			// std::cout << endl;
		}
		fp.close();
		fp.open(file_RE, std::ios::in); // ios::in means read
		if (!fp.is_open())
		{
			std::cerr << "This file READING for RE_ have some problem!!!\n";
		}
		for (int i = 0; i < num_ne; i++)
		{
			for (int j = 0; j < num_na; j++)
			{
				fp >> RE_[i][j];
				RE_[i][j] = log(RE_[i][j]);
				// std::cout << Rn_[i][j] << '\t';
			}
			// std::cout << endl;
		}
		fp.close();
	}

	/*for (int i = 0; i < num_ne; i++)
		std::cout << exp(Rn_ne_[i]) << endl;
	for (int i = 0; i < num_na; i++)
		std::cout << Rn_na_[i] << endl;
	for (int i = 0; i < num_ne; i++)
	{
		for (int j = 0; j < num_na; j++)
		{
			std::cout << exp(Rn_[i][j]) << '\t';
		}
		std::cout << endl;
	}*/
}

/// @brief K_Reflect: 1 for empirical formula; 2 for Trim database
double Reflect::n_RefCoeff(int K_Reflect, double Energy_0, double angle)
{
	if (K_Reflect == 1)
	{
		epsilon_ = Energy_0 / epsilon_L_;
		if (Energy_0 <= 1.e5)
			return std::min(1.0, n_a_[0] * pow(epsilon_, n_a_[1]) / (1 + n_a_[2] * pow(epsilon_, n_a_[3])));
		else
			return std::min(1.0, n_b_[0] * pow(epsilon_, n_b_[1]) / (1 + n_b_[2] * pow(epsilon_, n_b_[3])) * (1 + n_b_[4] * pow(epsilon_, n_b_[5])));
	}
	if (K_Reflect == 2 || K_Reflect == 3)
	{
		if (Energy_0 < 0.99)
			return 0;
		Energy_0 = log(Energy_0);
		int a, b, c, d, k;
		double x2x1, y2y1, x2x, y2y, yy1, xx1;
		a = BinarySearch(Rn_ne_, 0, num_ne_ - 1, Energy_0);
		b = BinarySearch(Rn_na_, 0, num_na_ - 1, angle);
		// std::cerr << a << '\t' << b << endl;
		if (a == -1)
		{
			if (Energy_0 <= Rn_ne_[0])
			{
				a = 0;
			}
			else if (Energy_0 >= Rn_ne_[num_ne_ - 1])
			{
				a = num_ne_ - 2;
				Energy_0 = Rn_ne_[num_ne_ - 1];
			}
			else
			{
				std::cerr << Energy_0 << endl;
				cerr << "BinartSearch for ne of n in Reflect has some problem." << endl;
			}
		}
		if (b == -1)
		{
			if (angle > Rn_na_[num_na_ - 1] && angle <= 90.1)
			{
				b = num_na_ - 2;
				angle = Rn_na_[num_na_ - 1];
			}
			else
			{
				std::cerr << "angle and Rn_na_[num_na_ - 1]: " << angle << "\t" << Rn_na_[num_na_ - 1] << endl;
				std::cerr << "BinartSearch of angle in n_Reflect has some problem." << endl;
			}
			if (angle < Rn_na_[0])
			{
				b = 0;
				angle = Rn_na_[0];
			}
			else if (angle > Rn_na_[num_na_ - 1])
			{
				b = num_na_ - 2;
				angle = Rn_na_[num_na_ - 1];
			}
			else if (b == -1)
				cerr << "BinartSearch in ADAS has some problem." << endl;
		}

		x2x1 = Rn_ne_[a + 1] - Rn_ne_[a];
		y2y1 = Rn_na_[b + 1] - Rn_na_[b];
		x2x = Rn_ne_[a + 1] - Energy_0;
		y2y = Rn_na_[b + 1] - angle;
		xx1 = Energy_0 - Rn_ne_[a];
		yy1 = angle - Rn_na_[b];

		return std::min(1.0, exp(1.0 / (x2x1 * y2y1) * (Rn_[a][b] * x2x * y2y + Rn_[a + 1][b] * xx1 * y2y + Rn_[a][b + 1] * x2x * yy1 + Rn_[a + 1][b + 1] * xx1 * yy1)));
	}
	throw std::invalid_argument("n_RefCoeff(): invalid K_Reflect value (must be 1, 2, or 3)");
}

/// @brief K_Reflect: 1 for empirical formula; 2 for Trim database
double Reflect::E_RefCoeff(int K_Reflect, double Energy_0, double angle)
{
	if (K_Reflect == 1)
	{
		epsilon_ = Energy_0 / epsilon_L_;
		if (Energy_0 <= 1.e5)
			return std::min(1.0, E_a_[0] * pow(epsilon_, E_a_[1]) / (1 + E_a_[2] * pow(epsilon_, E_a_[3])));
		else
			return std::min(1.0, E_b_[0] * pow(epsilon_, E_b_[1]) / (1 + E_b_[2] * pow(epsilon_, E_b_[3])) * (1 + E_b_[4] * pow(epsilon_, E_b_[5])));
	}
	else if (K_Reflect == 2 || K_Reflect == 3)
	{
		Energy_0 = log(Energy_0);
		int a, b;
		double x2x1, y2y1, x2x, y2y, yy1, xx1;
		a = BinarySearch(RE_ne_, 0, num_ne_ - 1, Energy_0);
		b = BinarySearch(RE_na_, 0, num_na_ - 1, angle);

		if (a == -1)
		{
			if (Energy_0 < RE_ne_[0])
			{
				a = 0;
				// std::cout << Energy_0 << endl;
				// c = 1;
			}
			else if (Energy_0 > RE_ne_[num_ne_ - 1])
			{
				a = num_ne_ - 2;
				Energy_0 = RE_ne_[num_ne_ - 1];
			}
			else if (a == -1)
				cerr << "BinartSearch for ne of E in Reflect has some problem." << endl;
		}
		if (b == -1)
		{
			if (angle > Rn_na_[num_na_ - 1] && angle <= 90.1)
			{
				b = num_na_ - 2;
				angle = Rn_na_[num_na_ - 1];
			}
			else
			{
				std::cerr << "angle and Rn_na_[num_na_ - 1]: " << angle << "\t" << Rn_na_[num_na_ - 1] << endl;
				cerr << "BinartSearch of angle in E_Reflect has some problem." << endl;
			}
			if (angle < RE_na_[0])
			{
				b = 0;
				angle = RE_na_[0];
			}
			else if (angle > RE_na_[num_na_ - 1])
			{
				b = num_na_ - 2;
				angle = RE_na_[num_na_ - 1];
			}
			else if (b == -1)
				cerr << "BinartSearch in ADAS has some problem." << endl;
		}

		x2x1 = RE_ne_[a + 1] - RE_ne_[a];
		y2y1 = RE_na_[b + 1] - RE_na_[b];
		x2x = RE_ne_[a + 1] - Energy_0;
		y2y = RE_na_[b + 1] - angle;
		xx1 = Energy_0 - RE_ne_[a];
		yy1 = angle - RE_na_[b];
		return std::min(1.0, exp(1.0 / (x2x1 * y2y1) * (RE_[a][b] * x2x * y2y + RE_[a + 1][b] * xx1 * y2y + RE_[a][b + 1] * x2x * yy1 + RE_[a + 1][b + 1] * xx1 * yy1)));
	}
	throw std::invalid_argument("n_RefCoeff(): invalid K_Reflect value (must be 1, 2, or 3)");
}

int BinarySearch(const vector<double> &array, int low, int high, double key)
{
	if (low < 0 || high <= low || high >= static_cast<int>(array.size()))
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

namespace
{
constexpr std::array<double, 7> kQuantileGrid{{0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0}};

double clampValue(double value, double low, double high)
{
	return std::max(low, std::min(high, value));
}

double linear(double a, double b, double weight)
{
	return a + weight * (b - a);
}

double bilinear(double q00, double q10, double q01, double q11,
				double energy_weight, double angle_weight)
{
	return linear(linear(q00, q10, energy_weight),
				  linear(q01, q11, energy_weight), angle_weight);
}

std::vector<double> numbersInLine(const std::string &line)
{
	std::vector<double> values;
	std::istringstream input(line);
	double value = 0.0;
	while (input >> value)
		values.push_back(value);
	return values;
}
}

double DWTrimReflection::SampleQuantile(const std::array<double, 7> &values, double xi)
{
	xi = clampValue(xi, 0.0, 1.0);
	const int interval = FindInterval(kQuantileGrid.data(), static_cast<int>(kQuantileGrid.size()), xi);
	const double width = kQuantileGrid[interval + 1] - kQuantileGrid[interval];
	const double weight = width > 0.0 ? (xi - kQuantileGrid[interval]) / width : 0.0;
	return linear(values[interval], values[interval + 1], weight);
}

int DWTrimReflection::FindInterval(const double *grid, int size, double value)
{
	if (size < 2)
		throw std::invalid_argument("DWTrimReflection grid must contain at least two points");
	if (value <= grid[0])
		return 0;
	if (value >= grid[size - 1])
		return size - 2;
	const double *upper = std::upper_bound(grid, grid + size, value);
	return static_cast<int>(upper - grid) - 1;
}

const DWTrimReflection::Block &DWTrimReflection::BlockAt(int energy_index, int angle_index) const
{
	return blocks_[energy_index][angle_index];
}

bool DWTrimReflection::Load(const std::string &filename)
{
	std::ifstream input(filename);
	if (!input)
	{
		std::cerr << "Cannot open D-on-W TRIM reflection database: " << filename << '\n';
		return false;
	}

	std::vector<Block> parsed;
	std::string line;
	while (std::getline(input, line))
	{
		const std::vector<double> header = numbersInLine(line);
		if (header.size() < 10 || std::abs(header[0] - 1.0) > 1e-8 ||
			std::abs(header[1] - 2.01) > 1e-6 || std::abs(header[2] - 74.0) > 1e-8)
			continue;

		Block block;
		block.incident_energy_eV = header[4];
		block.incident_angle_deg = header[5];
		block.particle_reflection = header[6];
		block.energy_reflection = header[7];

		auto readDataRow = [&input, &line]() {
			while (std::getline(input, line))
			{
				const std::vector<double> values = numbersInLine(line);
				if (!values.empty())
					return values;
			}
			return std::vector<double>{};
		};

		std::vector<double> values = readDataRow();
		if (values.size() != 7)
			throw std::runtime_error("Invalid energy row in D-on-W TRIM database");
		for (int q = 0; q < 5; ++q)
			block.energy[q + 1] = values[q];
		block.energy[0] = values[5];
		block.energy[6] = values[6];

		for (int energy_q = 0; energy_q < 5; ++energy_q)
		{
			values = readDataRow();
			if (values.size() != 7)
				throw std::runtime_error("Invalid polar-angle row in D-on-W TRIM database");
			for (int q = 0; q < 5; ++q)
				block.cos_polar[energy_q][q + 1] = values[q];
			block.cos_polar[energy_q][0] = values[5];
			block.cos_polar[energy_q][6] = values[6];
		}

		for (int energy_q = 0; energy_q < 5; ++energy_q)
			for (int polar_q = 0; polar_q < 5; ++polar_q)
			{
				values = readDataRow();
				if (values.size() != 7)
					throw std::runtime_error("Invalid azimuth row in D-on-W TRIM database");
				for (int q = 0; q < 5; ++q)
					block.cos_azimuth[energy_q][polar_q][q + 1] = values[q];
				block.cos_azimuth[energy_q][polar_q][0] = values[5];
				block.cos_azimuth[energy_q][polar_q][6] = values[6];
			}
		parsed.push_back(block);
	}

	if (parsed.size() != 84)
	{
		std::cerr << "D-on-W TRIM database contains " << parsed.size()
				  << " blocks; expected 84\n";
		return false;
	}

	for (int energy_index = 0; energy_index < 12; ++energy_index)
	{
		energies_[energy_index] = parsed[energy_index * 7].incident_energy_eV;
		for (int angle_index = 0; angle_index < 7; ++angle_index)
		{
			const Block &block = parsed[energy_index * 7 + angle_index];
			if (energy_index == 0)
				angles_[angle_index] = block.incident_angle_deg;
			if (std::abs(block.incident_energy_eV - energies_[energy_index]) > 1e-8 ||
				std::abs(block.incident_angle_deg - angles_[angle_index]) > 1e-8)
				throw std::runtime_error("Unexpected D-on-W TRIM energy/angle block ordering");
			blocks_[energy_index][angle_index] = block;
		}
	}

	loaded_ = true;
	return true;
}

bool DWTrimReflection::IsLoaded() const
{
	return loaded_;
}

double DWTrimReflection::ReflectionProbability(double incident_energy_eV,
												double incident_angle_deg) const
{
	if (!loaded_)
		throw std::runtime_error("D-on-W TRIM reflection database is not loaded");

	const double energy = clampValue(incident_energy_eV, energies_.front(), energies_.back());
	const double angle = clampValue(incident_angle_deg, angles_.front(), angles_.back());
	const int ei = FindInterval(energies_.data(), static_cast<int>(energies_.size()), energy);
	const int ai = FindInterval(angles_.data(), static_cast<int>(angles_.size()), angle);
	const double log_energy = std::log(energy);
	const double energy_weight =
		(log_energy - std::log(energies_[ei])) /
		(std::log(energies_[ei + 1]) - std::log(energies_[ei]));
	const double angle_weight = (angle - angles_[ai]) / (angles_[ai + 1] - angles_[ai]);

	return clampValue(
		bilinear(BlockAt(ei, ai).particle_reflection,
				 BlockAt(ei + 1, ai).particle_reflection,
				 BlockAt(ei, ai + 1).particle_reflection,
				 BlockAt(ei + 1, ai + 1).particle_reflection,
				 energy_weight, angle_weight),
		0.0, 1.0);
}

double DWTrimReflection::MeanReflectedEnergy(double incident_energy_eV,
											  double incident_angle_deg) const
{
	if (!loaded_)
		throw std::runtime_error("D-on-W TRIM reflection database is not loaded");

	const double energy = clampValue(incident_energy_eV, energies_.front(), energies_.back());
	const double angle = clampValue(incident_angle_deg, angles_.front(), angles_.back());
	const int ei = FindInterval(energies_.data(), static_cast<int>(energies_.size()), energy);
	const int ai = FindInterval(angles_.data(), static_cast<int>(angles_.size()), angle);
	const double energy_weight =
		(std::log(energy) - std::log(energies_[ei])) /
		(std::log(energies_[ei + 1]) - std::log(energies_[ei]));
	const double angle_weight = (angle - angles_[ai]) / (angles_[ai + 1] - angles_[ai]);

	auto blockMean = [](const Block &block) {
		double sum = 0.0;
		for (int q = 1; q <= 5; ++q)
			sum += block.energy[q];
		return sum / 5.0;
	};
	return clampValue(
		bilinear(blockMean(BlockAt(ei, ai)), blockMean(BlockAt(ei + 1, ai)),
				 blockMean(BlockAt(ei, ai + 1)), blockMean(BlockAt(ei + 1, ai + 1)),
				 energy_weight, angle_weight),
		0.0, incident_energy_eV);
}

DWReflectionSample DWTrimReflection::Sample(double incident_energy_eV,
											 double incident_angle_deg,
											 double xi_energy,
											 double xi_polar,
											 double xi_azimuth) const
{
	if (!loaded_)
		throw std::runtime_error("D-on-W TRIM reflection database is not loaded");

	const double energy = clampValue(incident_energy_eV, energies_.front(), energies_.back());
	const double angle = clampValue(incident_angle_deg, angles_.front(), angles_.back());
	const int ei = FindInterval(energies_.data(), static_cast<int>(energies_.size()), energy);
	const int ai = FindInterval(angles_.data(), static_cast<int>(angles_.size()), angle);
	const double energy_weight =
		(std::log(energy) - std::log(energies_[ei])) /
		(std::log(energies_[ei + 1]) - std::log(energies_[ei]));
	const double angle_weight = (angle - angles_[ai]) / (angles_[ai + 1] - angles_[ai]);

	auto sampleBlock = [xi_energy, xi_polar, xi_azimuth](const Block &block) {
		DWReflectionSample sample;
		sample.energy_eV = SampleQuantile(block.energy, xi_energy);

		constexpr std::array<double, 5> conditional_grid{{0.1, 0.3, 0.5, 0.7, 0.9}};
		const double energy_xi = clampValue(xi_energy, conditional_grid.front(), conditional_grid.back());
		const int energy_q = FindInterval(conditional_grid.data(), 5, energy_xi);
		const double energy_q_weight =
			(energy_xi - conditional_grid[energy_q]) /
			(conditional_grid[energy_q + 1] - conditional_grid[energy_q]);
		sample.cos_polar = linear(
			SampleQuantile(block.cos_polar[energy_q], xi_polar),
			SampleQuantile(block.cos_polar[energy_q + 1], xi_polar),
			energy_q_weight);

		const double polar_xi = clampValue(xi_polar, conditional_grid.front(), conditional_grid.back());
		const int polar_q = FindInterval(conditional_grid.data(), 5, polar_xi);
		const double polar_q_weight =
			(polar_xi - conditional_grid[polar_q]) /
			(conditional_grid[polar_q + 1] - conditional_grid[polar_q]);
		const double azimuth_low = linear(
			SampleQuantile(block.cos_azimuth[energy_q][polar_q], xi_azimuth),
			SampleQuantile(block.cos_azimuth[energy_q][polar_q + 1], xi_azimuth),
			polar_q_weight);
		const double azimuth_high = linear(
			SampleQuantile(block.cos_azimuth[energy_q + 1][polar_q], xi_azimuth),
			SampleQuantile(block.cos_azimuth[energy_q + 1][polar_q + 1], xi_azimuth),
			polar_q_weight);
		sample.cos_azimuth = linear(azimuth_low, azimuth_high, energy_q_weight);
		return sample;
	};

	const DWReflectionSample q00 = sampleBlock(BlockAt(ei, ai));
	const DWReflectionSample q10 = sampleBlock(BlockAt(ei + 1, ai));
	const DWReflectionSample q01 = sampleBlock(BlockAt(ei, ai + 1));
	const DWReflectionSample q11 = sampleBlock(BlockAt(ei + 1, ai + 1));

	DWReflectionSample result;
	result.energy_eV = clampValue(
		bilinear(q00.energy_eV, q10.energy_eV, q01.energy_eV, q11.energy_eV,
				 energy_weight, angle_weight),
		0.0, incident_energy_eV);
	result.cos_polar = clampValue(
		bilinear(q00.cos_polar, q10.cos_polar, q01.cos_polar, q11.cos_polar,
				 energy_weight, angle_weight),
		0.0, 1.0);
	result.cos_azimuth = clampValue(
		bilinear(q00.cos_azimuth, q10.cos_azimuth, q01.cos_azimuth, q11.cos_azimuth,
				 energy_weight, angle_weight),
		-1.0, 1.0);
	return result;
}
