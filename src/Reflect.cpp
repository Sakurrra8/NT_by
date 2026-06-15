#include "Reflect.h"

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