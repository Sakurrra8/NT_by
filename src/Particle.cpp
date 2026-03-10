#include "Particle.h"

Particle::Particle()
{
	V_.resize(3, 0.0);
}

Particle::Particle(string particle_name, double mass, int MaxCharge, int Index)
{
	name_ = particle_name;
	mass_ = mass;
	Index_ = Index;
	Zone_ = 0;
	MaxCharge_ = MaxCharge;
	ChargeTag_ = 0;
	Num_array_ = 0;
	V_.resize(3, 0.0);
	/*for (double i = 0.1; i <= 0.9; i = i + 0.1)
	{
		T_array_.push_back(i);
		Num_array_ += 1;
		// std::cout << i << '\t' << Num_array_ << endl;
	}*/
	// std::cout << endl;
	for (int i = 1; i <= 50; i++)
	{
		T_array_.push_back((double)i);
		Num_array_ += 1;
		// std::cout << i << '\t' << Num_array_ << endl;
	}
	for (int i = 0; i < 98; i++)
	{
		n_Flux_Grid_.push_back(vector<vector<double>>());
		T_Flux_Grid_.push_back(vector<vector<double>>());
		First_CX_.push_back(vector<double>());
		for (int j = 0; j < 38; j++)
		{
			n_Flux_Grid_[i].push_back(vector<double>());
			T_Flux_Grid_[i].push_back(vector<double>());
			First_CX_[i].push_back(double());
			for (int k = 0; k < 4; k++)
			{
				n_Flux_Grid_[i][j].push_back(double());
				T_Flux_Grid_[i][j].push_back(double());
			}
			n_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			T_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			E_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			Sum_n_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			Sum_E_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			lambda_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			for (int k = 0; k < MaxCharge_ + 1; k++)
			{
				CollProb_[i][j].push_back(0);
				Totalcs_[i][j].push_back(0);
				Sum_n_array_[i][j].push_back(vector<double>());
				for (int m = 0; m < Num_array_ + 1; m++)
				{
					Sum_n_array_[i][j][k].push_back(0);
				}
			}

			// CollProb_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			// Totalcs_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
			for (int k = 0; k < 3; k++)
			{
				V_Grid_CX_Ion_Be_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				V_Grid_CX_Ion_Af_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				V_Grid_CX_Neu_Be_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				V_Grid_CX_Neu_Af_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				V_Grid_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_CX_Ion_Be_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_CX_Ion_Af_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_CX_Neu_Be_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_CX_Neu_Af_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));

				V_D_1_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Sum_V_D_1_[i][j][k] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				Num_V_D_1_[i][j] = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));
				for (int m = 0; m < MaxCharge_ + 1; m++)
				{
					Sum_V_[i][j][k][m] = 0.;
					Sum_V_CX_Ion_Be_[i][j][k][m] = 0;
					Sum_V_CX_Ion_Af_[i][j][k][m] = 0;
					Sum_V_CX_Neu_Be_[i][j][k][m] = 0;
					Sum_V_CX_Neu_Af_[i][j][k][m] = 0;
					V_Grid_[i][j][k][m] = 0.;
					V_Grid_CX_Ion_Be_[i][j][k][m] = 0;
					V_Grid_CX_Ion_Af_[i][j][k][m] = 0;
					V_Grid_CX_Neu_Be_[i][j][k][m] = 0;
					V_Grid_CX_Neu_Af_[i][j][k][m] = 0;
					V_D_1_[i][j][k][m] = 0;
					Sum_V_D_1_[i][j][k][m] = 0;
					Num_V_D_1_[i][j][m] = 0;
				}
			}
			NumPar_Grid_[i][j] = 0;
		}
	}

	lambda_min_ = (double *)malloc(sizeof(double) * (MaxCharge_ + 1));

	Pump_.push_back(0);
	Pump_.push_back(0);
	Pump_.push_back(0);

	FT_.FluxTraInit(10, 200, 0.0, 1500, 98, 38, true, true); // 初始化 nRegion, nBin, eMin, eMax
}

Particle::~Particle()
{
	for (int i = 0; i < 98; i++)
		for (int j = 0; j < 38; j++)
		{
			free(n_[i][j]);
			free(T_[i][j]);
			free(Sum_n_[i][j]);
			free(lambda_[i][j]);
			for (int k = 0; k < 3; k++)
			{
				free(V_Grid_[i][j][k]);
				free(Sum_V_[i][j][k]);
			}
		}
	free(lambda_min_);
}

/// @brief Particle *this get Particle *A's All
/// @param A old Particle *
/// @param K *this Weight_/*A's Weight_, Default is 1
/// @param Charge *this Charge, Default is *A's
void Particle::Particlefrom(Particle *A, double K, int Charge)
{
	for (int i = 0; i < 3; i++)
	{
		X_[i] = A->X(i);
		X_new_[i] = A->X_new(i);
		V_[i] = A->V(i);
		// V_[i] = 0;
	}
	for (int i = 0; i < 8; i++)
	{
		V_Charge_[i] = A->V_Charge(i);
		// V_Charge_[i] = 0;
	}
	Zone_ = A->Zone();
	// Zone_new_ = A->Zone_new();
	XY_[0] = A->XY(0);
	XY_[1] = A->XY(1);
	Tri_Index_ = A->Tri_Index();

	Tn_ = A->Tn();

	// XY_new_[0] = A->XY_new(0);
	// XY_new_[1] = A->XY_new(1);
	boundary_start_ = A->boundary_start();
	//
	Weight_ = K * A->Weight();
	if (Charge == -100)
		Charge_ = A->Charge();
	else
		Charge_ = Charge;
}

void Particle::ParInit(vector<int> &K_collision, std::vector<std::vector<int>> Tri_B2)
{
	num_trimesh_ = Tri_B2.size();
	Tri_Totalcs_.resize(num_trimesh_);
	Tri_CollProb_.resize(num_trimesh_);
	Tri_lambda_.resize(num_trimesh_);

	Tri_n_.resize(num_trimesh_);
	Tri_T_.resize(num_trimesh_);
	Tri_E_.resize(num_trimesh_);
	Tri_V_Grid_.resize(num_trimesh_);
	Tri_V_D_1_.resize(num_trimesh_);
	Tri_Num_V_D_1_.resize(num_trimesh_);
	Tri_Sum_n_.resize(num_trimesh_);
	Tri_Sum_E_.resize(num_trimesh_);
	Tri_Sum_V_.resize(num_trimesh_);
	Tri_Sum_V_D_1_.resize(num_trimesh_);
	Tri_NumPar_Grid_.resize(num_trimesh_);

	for (int i = 0; i < num_trimesh_; i++)
	{
		Tri_n_[i].resize(MaxCharge_ + 1);
		Tri_T_[i].resize(MaxCharge_ + 1);
		Tri_E_[i].resize(MaxCharge_ + 1);
		Tri_Sum_n_[i].resize(MaxCharge_ + 1);
		Tri_Sum_E_[i].resize(MaxCharge_ + 1);
		Tri_Num_V_D_1_[i].resize(MaxCharge_ + 1);
		Tri_V_Grid_[i].resize(3);
		Tri_V_D_1_[i].resize(3);
		Tri_Sum_V_[i].resize(3);
		Tri_Sum_V_D_1_[i].resize(3);

		for (int j = 0; j < 3; j++)
		{
			Tri_V_Grid_[i][j].resize(MaxCharge_ + 1);
			Tri_V_D_1_[i][j].resize(MaxCharge_ + 1);
			Tri_Sum_V_[i][j].resize(MaxCharge_ + 1);
			Tri_Sum_V_D_1_[i][j].resize(MaxCharge_ + 1);
		}

		Tri_Totalcs_[i].resize(MaxCharge_ + 1);
		Tri_CollProb_[i].resize(MaxCharge_ + 1);
		Tri_lambda_[i].resize(MaxCharge_ + 1);

		Tri_B2_.push_back(std::vector<int>());
		if (Tri_B2[i][9] == -2)
		{
			Tri_B2_[i].push_back(-1);
			Tri_B2_[i].push_back(-1);
		}
		else
		{
			Tri_B2_[i].push_back(Tri_B2[i][9]);
			Tri_B2_[i].push_back(Tri_B2[i][10]);
		}
	}
	/*ofstream out_temp("doc/neigh.txt");
	for (int i = 0; i < num_trimesh_; i++)
	{
		out_temp << i << "\t" << Tri_B2_[i][0] << "\t" << Tri_B2_[i][1] << endl;
	}
	out_temp.close();*/

	for (int i = 0; i < MaxCharge_ + 1; i++)
	{
		CoreCollProb_.push_back(0);
		CoreTotalcs_.push_back(0);

		Ion_.push_back(ParCollCar());
		Rec_.push_back(ParCollCar());
		CX_.push_back(ParCollCar());
		Ion_[i].initialize(1, num_trimesh_);
		Rec_[i].initialize(1, num_trimesh_);
		CX_[i].initialize(2, num_trimesh_);
	}

	if (this == &H2 || this == &D2 || this == &T2 || this == &H || this == &D || this == &T)
	{
		R_with_H_.push_back(ParCollCar());
		R_with_H2_.push_back(ParCollCar());
		CX_DT_.push_back(ParCollCar());

		CX_DT_[0].initialize(2, num_trimesh_);
		R_with_H_[0].initialize(3, num_trimesh_);
		R_with_H2_[0].initialize(3, num_trimesh_);
	}

	if (this == &H2 || this == &D2 || this == &T2)
	{
		Ela_.push_back(ParCollCar());
		MAR_.push_back(ParCollCar());
		Diss1_.push_back(ParCollCar());
		Diss2_.push_back(ParCollCar());
		Diss3_.push_back(ParCollCar());

		Ela_[0].initialize(2, num_trimesh_);
		MAR_[0].initialize(1, num_trimesh_);
		Diss1_[0].initialize(1, num_trimesh_);
		Diss2_[0].initialize(1, num_trimesh_);
		Diss3_[0].initialize(1, num_trimesh_);
		Ela_DT_.push_back(ParCollCar());
		Ela_DT_[0].initialize(2, num_trimesh_);

		DS_.resize(2);
		DS_[1].push_back(ParCollCar());
		DS_[1].push_back(ParCollCar());
		DS_[1].push_back(ParCollCar());
		DS_[1].push_back(ParCollCar());
		DS_[1][0].initialize(1, num_trimesh_);
		DS_[1][1].initialize(1, num_trimesh_);
		DS_[1][2].initialize(1, num_trimesh_);
		DS_[1][3].initialize(1, num_trimesh_);
	}
	if (K_Methane)
	{
		DS_.resize(2);
		DS_[1].resize(4);
		DS_[1][0].initialize(1, num_trimesh_);
		DS_[1][1].initialize(1, num_trimesh_);
		DS_[1][2].initialize(1, num_trimesh_);
		DS_[1][3].initialize(1, num_trimesh_);
	}

	AlltoWall_.push_back(WallEro());
	CXtoWall_.push_back(WallEro());
	RectoWall_.push_back(WallEro());
	DisstoWall_.push_back(WallEro());
	AlltoWall_.begin()->WallEroInit(Grid1.Wall_num(), N_poloidal, N_radial);
	CXtoWall_.begin()->WallEroInit(Grid1.Wall_num(), N_poloidal, N_radial);
	RectoWall_.begin()->WallEroInit(Grid1.Wall_num(), N_poloidal, N_radial);
	DisstoWall_.begin()->WallEroInit(Grid1.Wall_num(), N_poloidal, N_radial);

	AlltoTarget_.push_back(WallEro());
	CXtoTarget_.push_back(WallEro());
	RectoTarget_.push_back(WallEro());
	DisstoTarget_.push_back(WallEro());
	AlltoTarget_.begin()->WallEroInit(N_radial * 2, N_poloidal, N_radial);
	CXtoTarget_.begin()->WallEroInit(N_radial * 2, N_poloidal, N_radial);
	RectoTarget_.begin()->WallEroInit(N_radial * 2, N_poloidal, N_radial);
	DisstoTarget_.begin()->WallEroInit(N_radial * 2, N_poloidal, N_radial);

	AlltoPlasmaBoundary_.push_back(WallEro());
	CXtoPlasmaBoundary_.push_back(WallEro());
	RectoPlasmaBoundary_.push_back(WallEro());
	DisstoPlasmaBoundary_.push_back(WallEro());
	AlltoPlasmaBoundary_.begin()->WallEroInit(Grid1.PLasma_Grid_Boundry_num(), N_poloidal, N_radial);
	CXtoPlasmaBoundary_.begin()->WallEroInit(Grid1.PLasma_Grid_Boundry_num(), N_poloidal, N_radial);
	RectoPlasmaBoundary_.begin()->WallEroInit(Grid1.PLasma_Grid_Boundry_num(), N_poloidal, N_radial);
	DisstoPlasmaBoundary_.begin()->WallEroInit(Grid1.PLasma_Grid_Boundry_num(), N_poloidal, N_radial);
}
void Particle::Vinit_ionReflect(std::vector<double> &Vt, double Rcos, double Rsin, int z)
{
	double ionIncidentTheta = pi / 3;
	if (K_ReflectDirection == 1)
		do
		{
			double ztheta = 6.283185307 * Random();
			double zcost = cos(ztheta);
			double zsint = sin(ztheta);
			double zpsi = Random() * 3.1415927;
			double zcosp = abs(cos(zpsi));
			double zsinp = sin(zpsi);
			double Vx = zsint * zsinp;
			double Vy = zcosp;
			double Vz = zsinp * zcost;

			Vt[0] = Rcos * Vx - Rsin * Vy;
			Vt[1] = Rsin * Vx + Rcos * Vy;
			Vt[2] = Vz;
		} while (Vt[2] == 1 || Vt[0] < -3e8);
	else if (K_ReflectDirection == 2)
	{
		double normalLine[3] = {-Rsin, Rcos, 0};
		double V_incident[3], projector[3], B_target[3];
		if (z < N_radial) // inner target
		{
			B_target[0] = B[1][(int)z][0] / sqrt(pow(B[1][(int)z][0], 2) + pow(B[1][(int)z][1], 2));
			B_target[1] = B[1][(int)z][0] / sqrt(pow(B[1][(int)z][0], 2) + pow(B[1][(int)z][1], 2));
		}
		else // outer target
		{
			B_target[0] = B[96][(int)z - N_radial][0] / sqrt(pow(B[96][(int)z - N_radial][0], 2) + pow(B[96][(int)z - N_radial][1], 2));
			B_target[1] = B[96][(int)z - N_radial][0] / sqrt(pow(B[96][(int)z - N_radial][0], 2) + pow(B[96][(int)z - N_radial][1], 2));
		}

		z = (B_target[0] * normalLine[0] + B_target[1] * normalLine[1]) * (B_target[1] * normalLine[0] - B_target[0] * normalLine[1]); // Determine which quadrant B is in the target plate
		if (z < 0)
		{
			double R = Random() * pi / 2.;
			V_incident[0] = -(normalLine[0] * cos(ionIncidentTheta) - normalLine[1] * sin(ionIncidentTheta)) * cos(R);
			V_incident[1] = -(normalLine[1] * cos(ionIncidentTheta) + normalLine[0] * sin(ionIncidentTheta)) * cos(R);
			V_incident[2] = sin(R);
		}
		else
		{
			double R = Random() * pi / 2.;
			V_incident[0] = (normalLine[0] * cos(ionIncidentTheta) - normalLine[1] * sin(ionIncidentTheta)) * cos(R);
			V_incident[1] = (normalLine[1] * cos(ionIncidentTheta) + normalLine[0] * sin(ionIncidentTheta)) * cos(R);
			V_incident[2] = sin(R);
		}

		do
		{
			double ztheta = 2 * pi * Random();
			double zcost = cos(ztheta);
			double zsint = sin(ztheta);
			double zpsi = Random() * pi;
			double zcosp = abs(cos(zpsi));
			double zsinp = sin(zpsi);
			double Vx = zsint * zsinp;
			double Vy = zcosp;
			double Vz = zsinp * zcost;

			Vt[0] = Rcos * Vx - Rsin * Vy;
			Vt[1] = Rsin * Vx + Rcos * Vy;
			Vt[2] = Vz;

			z = (V_[0] * Rcos + V_[1] * Rsin) * (Rcos * Vt[0] + Rcos * Vt[1]);
		} while (Vt[2] == 1 || Vt[0] < -3e8 || z < 0);
	}
	else if (K_ReflectDirection == 3)
	{
		double normalLine[3] = {-Rsin, Rcos, 0};
		double V_incident[3], projector[3], B_target[3];
		if (z < N_radial) // inner target
		{
			B_target[0] = B[1][(int)z][0] / sqrt(pow(B[1][(int)z][0], 2) + pow(B[1][(int)z][1], 2));
			B_target[1] = B[1][(int)z][1] / sqrt(pow(B[1][(int)z][0], 2) + pow(B[1][(int)z][1], 2));
		}
		else // outer target
		{
			B_target[0] = B[96][(int)z - N_radial][0] / sqrt(pow(B[96][(int)z - N_radial][0], 2) + pow(B[96][(int)z - N_radial][1], 2));
			B_target[1] = B[96][(int)z - N_radial][1] / sqrt(pow(B[96][(int)z - N_radial][0], 2) + pow(B[96][(int)z - N_radial][1], 2));
		}

		double z_temp = (B_target[0] * normalLine[0] + B_target[1] * normalLine[1]) * (B_target[1] * normalLine[0] - B_target[0] * normalLine[1]); // Determine which quadrant B is in the target plate
		if (z_temp < 0)
		{
			double R = Random() * pi / 2.;
			V_incident[0] = -(normalLine[0] * cos(ionIncidentTheta) - normalLine[1] * sin(ionIncidentTheta)) * cos(R);
			V_incident[1] = -(normalLine[1] * cos(ionIncidentTheta) + normalLine[0] * sin(ionIncidentTheta)) * cos(R);
			V_incident[2] = -sin(R);
		}
		else
		{
			double R = Random() * pi / 2.;
			V_incident[0] = (normalLine[0] * cos(ionIncidentTheta) - normalLine[1] * sin(ionIncidentTheta)) * cos(R);
			V_incident[1] = (normalLine[1] * cos(ionIncidentTheta) + normalLine[0] * sin(ionIncidentTheta)) * cos(R);
			V_incident[2] = -sin(R);
		}

		projector[0] = (V_incident[0] * Rcos + V_incident[1] * Rsin) * Rcos;
		projector[1] = (V_incident[0] * Rcos + V_incident[1] * Rsin) * Rsin;
		for (int i = 0; i < 2; i++)
		{
			Vt[i] = 2 * projector[i] - V_incident[i];
		}
		Vt[2] = V_incident[2];
	}
}

void Particle::Vinit_neutralReflect(std::vector<double> &Vt, double Rcos, double Rsin)
{
	if (K_ReflectDirection == 1)
		do
		{
			double ztheta = 6.283185307 * Random();
			double zcost = cos(ztheta);
			double zsint = sin(ztheta);
			double zpsi = Random() * 3.1415927;
			double zcosp = abs(cos(zpsi));
			double zsinp = sin(zpsi);
			double Vx = zsint * zsinp;
			double Vy = zcosp;
			double Vz = zsinp * zcost;

			Vt[0] = Rcos * Vx - Rsin * Vy;
			Vt[1] = Rsin * Vx + Rcos * Vy;
			Vt[2] = Vz;
		} while (Vt[2] == 1 || Vt[0] < -3e8);
	if (K_ReflectDirection == 2)
	{
		double z, projector[3], normalLine[3] = {-Rsin, Rcos, 0};
		do
		{
			double ztheta = 2 * pi * Random();
			double zcost = cos(ztheta);
			double zsint = sin(ztheta);
			double zpsi = Random() * pi;
			double zcosp = abs(cos(zpsi));
			double zsinp = sin(zpsi);
			double Vx = zsint * zsinp;
			double Vy = zcosp;
			double Vz = zsinp * zcost;

			Vt[0] = Rcos * Vx - Rsin * Vy;
			Vt[1] = Rsin * Vx + Rcos * Vy;
			Vt[2] = Vz;

			z = (V_[0] * Rcos + V_[1] * Rsin) * (Rcos * Vt[0] + Rcos * Vt[1]);
		} while (Vt[2] == 1 || Vt[0] < -3e8 || z < 0);
	}
	else if (K_ReflectDirection == 3)
	{
		double V_incident[3], projector[3];
		double normalLine[3] = {-Rsin, Rcos, 0};
		for (int i = 0; i < 3; i++)
		{
			V_incident[i] = V_[i] / sqrt(V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]);
		}
		// std::cout << name_ + " X: " << X_[0] << ", " << X_[1] << ", " << X_[2] << endl;
		// std::cout << name_ + " before reflect: " << V_incident[0] << ", " << V_incident[1] << ", " << V_incident[2] << endl;
		projector[0] = (V_incident[0] * Rcos + V_incident[1] * Rsin) * Rcos;
		projector[1] = (V_incident[0] * Rcos + V_incident[1] * Rsin) * Rsin;
		for (int i = 0; i < 2; i++)
		{
			Vt[i] = 2 * projector[i] - V_incident[i];
		}
		Vt[2] = V_incident[2];
		// std::cout << name_ + " after reflect: " << Vt[0] << ", " << Vt[1] << ", " << Vt[2] << endl;
	}
}

void Particle::Init(int k, int z)
{
	double vel = 0.;
	double deltaX, deltaY, deltaL, cosCX, sinCX;
	if (k == 1) // neutral particle from recycling
	{
		fate_[0] = 2;
		sourcePar_[0] = 1;
		Charge_ = 0;
		ChargeTag_ = 0;
		boundary_start_ = 0;
		X_[0] = Xrecyc[z][0];
		X_[1] = Xrecyc[z][1];
		// double rand_temp = Random();
		// X_[0] = Xrecycling(rand_temp, z, 0);
		// X_[1] = Xrecycling(rand_temp, z, 1);
		X_[2] = 0;

		Vector3 InciDir, P1, P2, P3;
		double Bt_temp, B2_temp;
		if (z < 38)
		{
			XY_[0] = 1;
			XY_[1] = z;
			Tri_Index_ = Grid4.targetIndex(z, 0);
			B2_temp = sqrt(B[1][z][0] * B[1][z][0] + B[1][z][1] * B[1][z][1]);
			Bt_temp = B[1][z][2];
			InciDir = {-0.354 * B2_temp, -0.612 * B2_temp, Bt_temp};
			P1 = {Grid[1][z][3], Grid[1][z][7], 0};
			P2 = {Grid[1][z][0], Grid[1][z][4], 0};
			P3 = {Grid[1][z][3], Grid[1][z][7], -1};
		}
		else
		{
			XY_[0] = 96;
			XY_[1] = z - 38;
			Tri_Index_ = Grid4.targetIndex(z, 0);
			B2_temp = sqrt(B[96][XY_[1]][0] * B[96][XY_[1]][0] + B[96][XY_[1]][1] * B[96][XY_[1]][1]);
			Bt_temp = B[96][XY_[1]][2];
			InciDir = {0.612 * B2_temp, 0.354 * B2_temp, Bt_temp};
			P1 = {Grid[96][z - 38][1], Grid[96][z - 38][5], 0};
			P2 = {Grid[96][z - 38][2], Grid[96][z - 38][6], 0};
			P3 = {Grid[96][z - 38][1], Grid[96][z - 38][5], -1};
		}
		if (XY_[1] >= 19 && XY_[1] <= 36)
			Zone_ = 3;
		else if (XY_[0] <= 24)
			Zone_ = 4;
		else if (XY_[0] <= 72)
			Zone_ = 2;
		else
			Zone_ = 5;
		GridIndex_ = XY_[0] * 38 + XY_[1];
		// std::cout << "Tri_Index: " << Tri_Index_ << endl;
		/*if (z < 38 && X_[0] < 1.361)
			X_[0] = 1.3611;*/
		/*Vector3 P1 = {Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z, 0};
		Vector3 P2 = {Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z, 0};
		Vector3 P3 = {Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].z, 0};*/

		if (this == &H)
		{
			Tn_ = Tn_H_recyc[z];
			// Vinit_ionReflect(V_, theta[z][0], theta[z][1], z);
			emit_particle(P1, P2, P3, 2, InciDir);
		}
		if (this == &D)
		{
			Tn_ = Tn_D_recyc[z];
			// Vinit_ionReflect(V_, -theta[z][0], -theta[z][1], z);
			emit_particle(P1, P2, P3, 2, InciDir);
			// std::cout << "first X: " << X_[0] << "\t" << X_[1] << "\t" << X_[2] << endl;
			// std::cout << "first V: " << V_[0] << "\t" << V_[1] << "\t" << V_[2] << endl;
		}
		if (this == &T)
		{
			Tn_ = Tn_T_recyc[z];
			Vinit_ionReflect(V_, theta[z][0], theta[z][1], z);
		}
		else if (this == &H2 || this == &D2 || this == &T2)
		{
			Tn_ = T_wall;
			// Vinit_cosine(V_, theta[z][0], theta[z][1]);
			emit_particle(P1, P2, P3, 1, InciDir);
		}
		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);

		if (K_test1)
		{
			std::cout << Weight_Target_[z] << endl;
		}
		Weight_ = Weight_Target_[z];
		NumPar_now = NumPar_Target_[z];
	}
	else if (k == 2) // neutral particle after CX for D or T
	{
		/*m_A for Ion, m_B for Neutral, R is the same
		unit vector with a random direction			*/
		double V_ion[3], V_temp[3], m_A, m_B, Modules;
		if (this == &H)
			m_A = Hmass;
		else if (this == &D || this == &D2)
		{
			m_A = Dmass;
			if (K_Vi == 1)
			{
				V_ion[0] = V_D_1_now[0];
				V_ion[1] = V_D_1_now[0];
				V_ion[2] = V_D_1_now[0];
			}
		}
		else if (this == &T || this == &T2)
		{
			m_A = Tmass;
			if (K_Vi == 1)
			{
				V_ion[0] = V_T_1_now[0];
				V_ion[1] = V_T_1_now[0];
				V_ion[2] = V_T_1_now[0];
			}
		}
		m_B = mass_;
		if (Charge_ == 0)
			for (int i = 0; i < 3; i++)
				V_temp[i] = V_[i];
		else // if the calculate is needful
		{
			pause();
			double Rand_temp1, Rand_temp2;
			Rand_temp1 = 2 * pi * Random();
			Rand_temp1 = pi * Random();
			double V_temp1[3], V_temp2[3];
			V_temp1[0] = sin(Rand_temp1) * sin(Rand_temp2);
			V_temp1[1] = sin(Rand_temp1) * cos(Rand_temp2);
			V_temp1[2] = cos(Rand_temp1);
			V_temp2[0] = V_temp1[1] * B[XY_[0]][XY_[1]][2] - V_temp1[2] * B[XY_[0]][XY_[1]][1];
			V_temp2[1] = V_temp1[2] * B[XY_[0]][XY_[1]][0] - V_temp1[0] * B[XY_[0]][XY_[1]][2];
			V_temp2[2] = V_temp1[0] * B[XY_[0]][XY_[1]][1] - V_temp1[1] * B[XY_[0]][XY_[1]][0];
			double Module = sqrt(pow(V_temp2[0], 2) + pow(V_temp2[1], 2) + pow(V_temp2[2], 2));
			V_temp[0] = V_Charge_[0] + V_Charge_[7] * V_temp2[0] / Module; // * erx[XY_[0]][XY_[1]];
			V_temp[1] = V_Charge_[1] + V_Charge_[7] * V_temp2[1] / Module; // * ery[XY_[0]][XY_[1]];
			V_temp[2] = V_Charge_[2] + V_Charge_[7] * V_temp2[2] / Module;
		}
		double Vtheta, Vphi;
		if (K_Vi == 2)
		{
			if (K_Maxwell == 1)
				vel = Maxwell(Ti[XY_[0]][XY_[1]], m_A);
			else if (K_Maxwell == 2)
				vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / m_A);
			Vtheta = 2 * pi * Random();
			Vphi = pi * Random();
			V_ion[0] = vel * sin(Vphi) * cos(Vtheta);
			V_ion[1] = vel * sin(Vphi) * sin(Vtheta);
			V_ion[2] = vel * cos(Vphi);
		}
		Vtheta = 2 * pi * Random();
		Vphi = pi * Random();
		Modules = sqrt(pow(V_temp[0] - V_ion[0], 2) + pow(V_temp[1] - V_ion[1], 2) + pow(V_temp[2] - V_ion[2], 2));
		// z = 1 for particle B, z = 2 for particle A
		V_1[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] - m_A * Modules * sin(Vphi) * cos(Vtheta));
		V_1[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] - m_A * Modules * sin(Vphi) * sin(Vtheta));
		V_1[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] - m_A * Modules * cos(Vphi));
		V_2[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] + m_B * Modules * sin(Vphi) * cos(Vtheta));
		V_2[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] + m_B * Modules * sin(Vphi) * sin(Vtheta));
		V_2[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] + m_B * Modules * cos(Vphi));
		if (K_test3)
		{
			std::cout << "2V of neu: " << V_1[0] << ", " << V_1[1] << ", " << V_1[2] << " T:"
					  << 1. / 2. * m_B * (pow(V_1[0], 2) + pow(V_1[1], 2) + pow(V_1[2], 2)) / qe << endl;
			std::cout << "2V of ion: " << V_2[0] << ", " << V_2[1] << ", " << V_2[2] << " T:"
					  << 1. / 2. * m_A * (pow(V_2[0], 2) + pow(V_2[1], 2) + pow(V_2[2], 2)) / qe << endl;
		}
		if (z == 1)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_1[i];
		}
		if (z == 2)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_2[i];
		}
		// setfate(0, "rec", "D");
		vel = 1.;
	}
	else if (k == 3) // neutral particle from Wall, z = 1 for
	{
		if (MeshMode == 1)
		{
			if (K_Wallelement == 1)
			{
				if (this == &D)
					Tn_ = Tn_ * D_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / D_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if (this == &H)
					Tn_ = Tn_ * H_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / H_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if (this == &T)
					Tn_ = Tn_ * T_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / T_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if ((this == &H2 || this == &D2 || this == &T2) && z == 1)
					Tn_ = T_wall;
			}
			else if (K_Wallelement == 2)
			{
				if (this == &D)
					Tn_ = Tn_ * D_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / D_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if (this == &H)
					Tn_ = Tn_ * H_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / H_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if (this == &T)
					Tn_ = Tn_ * T_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3])) / T_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle((int)InterscePoint[0][3]));
				else if ((this == &H2 || this == &D2 || this == &T2) && z == 1)
					Tn_ = T_wall;
			}
			X_[0] = InterscePoint[0][1];
			X_[1] = InterscePoint[0][2];
			X_[2] = 0;

			if (InterscePoint[0][4] == 11)
			{
				XY_[0] = -1;
				XY_[1] = -1;
				Zone_ = 6;
			}

			if (z == 0) // reflect of neutral particle
			{
				if (InterscePoint[0][4] == 11)
				{
					Vinit_neutralReflect(V_, cosang[(int)InterscePoint[0][3]], sinang[(int)InterscePoint[0][3]]);
				}
				else if (InterscePoint[0][4] == 1)
				{
					Vinit_neutralReflect(V_, Grid1.Cos_target((int)InterscePoint[0][3]), Grid1.Sin_target((int)InterscePoint[0][3]));
				}
				else
				{
					std::cout << "where is the reflect point???" << endl;
					pause();
				}
			}
			if (z == 1) // thermalrelease of neutral particle
			{
				if (InterscePoint[0][4] == 11)
				{
					Vinit_cosine(V_, cosang[(int)InterscePoint[0][3]], sinang[(int)InterscePoint[0][3]]);
				}
				else if (InterscePoint[0][4] == 1)
				{
					Vinit_cosine(V_, Grid1.Cos_target((int)InterscePoint[0][3]), Grid1.Sin_target((int)InterscePoint[0][3]));
				}
				else
				{
					std::cout << "where is the reflect point???" << endl;
					pause();
				}
				//  std::cout << Charge_ << '\t' << ChargeTag_ << endl;
				ChargeTag_ = 0;
				// std::cout << V_[0] << "\t" << V_[1] << "\t" << V_[2] << endl;
			}
			if (z == 2 || Charge_ > 0)
			{
				// std::cout << theta[(int)InterscePoint[0][3]][0] << '\t' << theta[(int)InterscePoint[0][3]][1] << endl;
				Vinit_ionReflect(V_, theta[(int)InterscePoint[0][3]][0], theta[(int)InterscePoint[0][3]][1], InterscePoint[0][3]);
				// std::cout << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
				ChargeTag_ = 0;
			}
		}
		else if (MeshMode == 3)
		{
			double Cos_temp = -Grid4.Cos_trimesh(Tri_Index_, (int)InterscePoint[0][3]);
			double Sin_temp = -Grid4.Sin_trimesh(Tri_Index_, (int)InterscePoint[0][3]);
			/*std::cout << "index: " << Tri_Index_ << endl;
			std::cout << "Inter[0][3]: " << (int)InterscePoint[0][3] << endl;
			std::cout << "jisuan: " << Grid4.Cos_trimesh(262, 2) << endl;
			std::cout << "sanjiao: " << Cos_temp << "\t" << Sin_temp << endl;*/
			if (K_Wallelement == 1)
			{
				if (this == &D)
					Tn_ = Tn_ * D_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / D_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if (this == &H)
					Tn_ = Tn_ * H_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / H_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if (this == &T)
					Tn_ = Tn_ * T_W.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / T_W.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if ((this == &H2 || this == &D2 || this == &T2) && z == 1)
					Tn_ = T_wall;
			}
			else if (K_Wallelement == 2)
			{
				if (this == &D)
					Tn_ = Tn_ * D_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / D_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if (this == &H)
					Tn_ = Tn_ * H_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / H_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if (this == &T)
					Tn_ = Tn_ * T_C.E_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp)) / T_C.n_RefCoeff(K_Reflect, 1.5 * Tn_, CalAngle(Sin_temp, Cos_temp));
				else if ((this == &H2 || this == &D2 || this == &T2) && z == 1)
					Tn_ = T_wall;
			}
			X_[0] = InterscePoint[0][1];
			X_[1] = InterscePoint[0][2];
			X_[2] = 0;

			if (InterscePoint[0][4] == 11)
				Zone_ = 6;
			if (InterscePoint[0][4] == 1)
				Zone_ = 4;

			if (z == 0) // reflect of neutral particle
			{
				if (InterscePoint[0][4] == 11 || InterscePoint[0][4] == 1)
				{
					Vinit_neutralReflect(V_, Cos_temp, Sin_temp);
				}
				else
				{
					std::cout << "where is the reflect point???" << endl;
					pause();
				}
			}
			if (z == 1) // thermalrelease of neutral particle
			{
				if (InterscePoint[0][4] == 11 || InterscePoint[0][4] == 1)
				{
					Vinit_cosine(V_, Cos_temp, Sin_temp);
				}
				else
				{
					std::cout << "where is the reflect point???" << endl;
					pause();
				}
				//  std::cout << Charge_ << '\t' << ChargeTag_ << endl;
				ChargeTag_ = 0;
				// std::cout << V_[0] << "\t" << V_[1] << "\t" << V_[2] << endl;
			}
			if (z == 2 || Charge_ > 0)
			{
				std::cout << "error way 1" << endl;
				std::cout << theta[(int)InterscePoint[0][3]][0] << '\t' << theta[(int)InterscePoint[0][3]][1] << endl;
				Vinit_ionReflect(V_, theta[(int)InterscePoint[0][3]][0], theta[(int)InterscePoint[0][3]][1], InterscePoint[0][3]);
				std::cout << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
				ChargeTag_ = 0;
				exit(0);
			}
		}

		Charge_ = 0;
		InterscePoint[0][0] = -1;
		InterscePoint[0][1] = -100;
		InterscePoint[0][2] = -100;
		InterscePoint[0][3] = -1;
		InterscePoint[0][4] = -1;

		/*std::cout << "After reflect:" << endl;
		std::cout << name() << ',' << Charge() << ',' << fatename(fate(0)) << ',' << sourcename(0) << endl;
		std::cout << "X," << X(0) << ',' << X(1) << ';' << X_new(0) << ',' << X_new(1) << endl;
		std::cout << "V," << V(0) << ',' << V(1) << endl;
		std::cout << endl;*/

		// std::cout << z << '\t' << X_[0] << '\t' << X_[1] << '\t' << X_[2] <<
		// '\t'; std::cout << '\t' << name_ << '\t' << Charge_ << '\t' << V_[0] <<
		// '\t' << V_[1] << '\t' << V_[2] << endl;
		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);
		// std::cout << name_ + ": " << Tn_ << " " << V_[0] << "," << V_[1] << "," << V_[2] << endl;
		// SetsourceWall(InterscePoint[0]);
		// std::cout << Tn_ << endl;
	}
	else if (k == 4) // neutral particle from Recombination
	{
		Charge_ = 0;
		fate_[0] = 4;
		sourcePar_[0] = 12;

		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);
		double Vtheta = 2 * pi * Random();
		double Vphi = pi * Random();
		/*VtoVcharge();
		VchargetoV();*/
		if (MeshMode == 3)
		{
			int a = Tri_B2_[z][0];
			int b = Tri_B2_[z][1];
			Tn_ = Ti[a][b];
			Tri_Index_ = z;
			Weight_ = Weight_Grid_[a][b];
			NumPar_now = Tri_NumPar_Grid_[z];
			// std::cout << "Weight: " << Weight_ << endl;
			X_[0] = (Grid4.nodes_[Grid4.tris_[z].v[0]].r + Grid4.nodes_[Grid4.tris_[z].v[1]].r + Grid4.nodes_[Grid4.tris_[z].v[2]].r) / 3.;
			X_[1] = (Grid4.nodes_[Grid4.tris_[z].v[0]].z + Grid4.nodes_[Grid4.tris_[z].v[1]].z + Grid4.nodes_[Grid4.tris_[z].v[2]].z) / 3.;
			X_[2] = 0;
			XY_[0] = a;
			XY_[1] = b;

			if (XY_[1] >= 19 && XY_[1] <= 36)
				Zone_ = 3;
			else if (XY_[0] <= 24)
				Zone_ = 4;
			else if (XY_[0] <= 72)
				Zone_ = 2;
			else
				Zone_ = 5;

			Rec_[1].n_Add(Tri_Index_, Weight_ * NumPar_now);
			if (this == &D)
			{
				V_[0] = sin(Vphi) * cos(Vtheta) + ua_D_1[a][b] * B[a][b][0] / vel;
				V_[1] = sin(Vphi) * sin(Vtheta) + ua_D_1[a][b] * B[a][b][1] / vel;
				V_[2] = cos(Vphi) + ua_D_1[a][b] * B[a][b][2] / vel;
				Rec_[1].Mu_Add(Tri_Index_, -mass_ * ua_D_1[a][b] * Weight_ * NumPar_now);
				Rec_[1].E_Add(Tri_Index_, -(3. / 2 * Ti[a][b] * qe + 1.0 / 2 * mass_ * ua_D_1[a][b] * ua_D_1[a][b]) * Weight_ * NumPar_now);
			}
			else if (this == &T)
			{
				V_[0] = sin(Vphi) * cos(Vtheta) + ua_T_1[a][b] * B[a][b][0] / vel;
				V_[1] = sin(Vphi) * sin(Vtheta) + ua_T_1[a][b] * B[a][b][1] / vel;
				V_[2] = cos(Vphi) + ua_T_1[a][b] * B[a][b][2] / vel;
				Rec_[1].Mu_Add(Tri_Index_, -mass_ * ua_T_1[a][b] * Weight_ * NumPar_now);
				Rec_[1].E_Add(Tri_Index_, -(3. / 2 * Ti[a][b] * qe + 1.0 / 2 * mass_ * ua_T_1[a][b] * ua_T_1[a][b]) * Weight_ * NumPar_now);
			}
			if (K_DT)
			{
				if (this == &D)
				{
					if (K_database_Pra == 1)
						Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_D_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
					if (K_database_Pra == 2)
						Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_D_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
				}
				if (this == &T)
				{
					if (K_database_Pra == 1)
						Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_T_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
					if (K_database_Pra == 2)
						Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_T_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
				}
			}
			else
			{
				if (K_database_Pra == 1)
					Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_D_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
				if (K_database_Pra == 2)
					Rec_[1].Pra_Add(Tri_Index_, ne[a][b] * n_D_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
			}
		}
		else
		{
			int a = z / N_radial;
			int b = z % N_radial;
			Tn_ = Ti[a][b];
			Weight_ = Weight_Grid_[a][b];
			NumPar_now = NumPar_Grid_[a][b];
			X_[0] = (Grid[a][b][0] + Grid[a][b][1] + Grid[a][b][2] + Grid[a][b][3]) / 4.;
			X_[1] = (Grid[a][b][4] + Grid[a][b][5] + Grid[a][b][6] + Grid[a][b][7]) / 4.;
			X_[2] = 0;
			XY_[0] = a;
			XY_[1] = b;

			if (XY_[1] >= 19 && XY_[1] <= 36)
				Zone_ = 3;
			else if (XY_[0] <= 24)
				Zone_ = 4;
			else if (XY_[0] <= 72)
				Zone_ = 2;
			else
				Zone_ = 5;

			int Index[2] = {a, b};
			Rec_[1].n_Add(Index, Weight_ * NumPar_now);
			if (this == &D)
			{
				V_[0] = sin(Vphi) * cos(Vtheta) + ua_D_1[a][b] * B[a][b][0] / vel;
				V_[1] = sin(Vphi) * sin(Vtheta) + ua_D_1[a][b] * B[a][b][1] / vel;
				V_[2] = cos(Vphi) + ua_D_1[a][b] * B[a][b][2] / vel;
				Rec_[1].Mu_Add(Index, -mass_ * ua_D_1[a][b] * Weight_ * NumPar_now);
				Rec_[1].E_Add(Index, -(3. / 2 * Ti[a][b] * qe + 1.0 / 2 * mass_ * ua_D_1[a][b] * ua_D_1[a][b]) * Weight_ * NumPar_now);
			}
			else if (this == &T)
			{
				V_[0] = sin(Vphi) * cos(Vtheta) + ua_T_1[a][b] * B[a][b][0] / vel;
				V_[1] = sin(Vphi) * sin(Vtheta) + ua_T_1[a][b] * B[a][b][1] / vel;
				V_[2] = cos(Vphi) + ua_T_1[a][b] * B[a][b][2] / vel;
				Rec_[1].Mu_Add(Index, -mass_ * ua_T_1[a][b] * Weight_ * NumPar_now);
				Rec_[1].E_Add(Index, -(3. / 2 * Ti[a][b] * qe + 1.0 / 2 * mass_ * ua_T_1[a][b] * ua_T_1[a][b]) * Weight_ * NumPar_now);
			}
			if (K_DT)
			{
				if (this == &D)
				{
					if (K_database_Pra == 1)
						Rec_[1].Pra_Add(Index, ne[a][b] * n_D_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
					if (K_database_Pra == 2)
						Rec_[1].Pra_Add(Index, ne[a][b] * n_D_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
				}
				if (this == &T)
				{
					if (K_database_Pra == 1)
						Rec_[1].Pra_Add(Index, ne[a][b] * n_T_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
					if (K_database_Pra == 2)
						Rec_[1].Pra_Add(Index, ne[a][b] * n_T_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
				}
			}
			else
			{
				if (K_database_Pra == 1)
					Rec_[1].Pra_Add(Index, ne[a][b] * n_D_1[a][b] * PRB12_H.cal(Te[a][b], ne[a][b]));
				if (K_database_Pra == 2)
					Rec_[1].Pra_Add(Index, ne[a][b] * n_D_1[a][b] * R2_1_8_H10.cal(ne[a][b], Te[a][b]) * qe);
			}
		}
	}
	else if (k == 5) // CD4 seeding
	{
		Charge_ = 0;
		ChargeTag_ = 0;
		fate_[0] = 3;
		sourcePar_[0] = 2;
		Tn_ = T_N;
		Weight_ = 1.;
		NumPar_now = Num_CD4_pump / (double)numPar_flight_CD4;
		if (this == &H2 || this == &D2 || this == &T2)
			NumPar_now *= 2;
		if (this == &H || this == &D || this == &T)
			NumPar_now *= 4;

		X_[0] = 1.75590002441406;
		X_[1] = (1.082 - 1.048) * Random() - 1.082;
		X_[2] = 0;

		XY_[0] = -1;
		XY_[1] = -1;
		Zone_ = 6;

		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);
		Vinit_cosine(V_, 0, 1);
		// std::cout << "fashe: " << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
	}
	else if (k == 6) // Calculate the velocity after the Collision
	{
		/// Conservation of energy after collision and conservation of momentum
		/*if (fate_[0] != 9 && (fate_[0] >= 7 || fate_[0] <= 14))
			if (sourcePar_[0] == 4)
			{
				double V_temp1[3], V_temp2[3] = {Random(), Random(), Random()};
				for (int i = 0; i < 3; i++)
					V_temp1[i] = intersect(V_temp2, V_, i);
				V_temp2[0] = sqrt(pow(V_temp1[0], 2) + pow(V_temp1[1], 2) + pow(V_temp1[2], 2));
				V_temp2[1] = sqrt(pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2));
				for (int i = 0; i < 3; i++)
					V_[i] = V_[i] + V_temp1[i] / V_temp2[0] * V_temp2[1];
				return;
			}*/

		/// Energy is conserved after a collision but momentum is not
		vel = sqrt(pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2));
		// std::cout << fate_[0] << '\t' << sourcePar_[0] << '\t';
		// std::cout << Tn_ << '\t' << vel << '\t';
		for (int i = 0; i < 3; i++)
			V_[i] = V_[i] / vel;
		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);
		// std::cout << vel << endl;
		if (Charge_ != 0)
			VtoVcharge();
	}
	else if (k == 7) // Cross-CX between D and T
	{
		/*m_A for Ion, m_B for Neutral, R is the same
		unit vector with a random direction			*/
		double V_ion[3], V_temp[3], m_A, m_B, Modules;
		if (this == &T || this == &T2)
		{
			m_A = Dmass;
			if (K_Vi == 1)
			{
				V_ion[0] = V_D_1_now[0];
				V_ion[1] = V_D_1_now[0];
				V_ion[2] = V_D_1_now[0];
			}
		}
		if (this == &D || this == &D2)
		{
			m_A = Tmass;
			if (K_Vi == 1)
			{
				V_ion[0] = V_T_1_now[0];
				V_ion[1] = V_T_1_now[0];
				V_ion[2] = V_T_1_now[0];
			}
		}
		m_B = mass_;
		if (Charge_ == 0)
			for (int i = 0; i < 3; i++)
				V_temp[i] = V_[i];
		else
		{
			pause();
			double Rand_temp1, Rand_temp2;
			Rand_temp1 = 2 * pi * Random();
			Rand_temp1 = pi * Random();
			double V_temp1[3], V_temp2[3];
			V_temp1[0] = sin(Rand_temp1) * sin(Rand_temp2);
			V_temp1[1] = sin(Rand_temp1) * cos(Rand_temp2);
			V_temp1[2] = cos(Rand_temp1);
			V_temp2[0] = V_temp1[1] * B[XY_[0]][XY_[1]][2] - V_temp1[2] * B[XY_[0]][XY_[1]][1];
			V_temp2[1] = V_temp1[2] * B[XY_[0]][XY_[1]][0] - V_temp1[0] * B[XY_[0]][XY_[1]][2];
			V_temp2[2] = V_temp1[0] * B[XY_[0]][XY_[1]][1] - V_temp1[1] * B[XY_[0]][XY_[1]][0];
			double Module = sqrt(pow(V_temp2[0], 2) + pow(V_temp2[1], 2) + pow(V_temp2[2], 2));
			V_temp[0] = V_Charge_[0] + V_Charge_[7] * V_temp2[0] / Module; // * erx[XY_[0]][XY_[1]];
			V_temp[1] = V_Charge_[1] + V_Charge_[7] * V_temp2[1] / Module; // * ery[XY_[0]][XY_[1]];
			V_temp[2] = V_Charge_[2] + V_Charge_[7] * V_temp2[2] / Module;
		}
		double Vtheta, Vphi;
		if (K_Vi == 2)
		{
			Vtheta = 2 * pi * Random();
			Vphi = pi * Random();
			if (K_Maxwell == 1)
				vel = Maxwell(Ti[XY_[0]][XY_[1]], m_A);
			else if (K_Maxwell == 2)
				vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / m_A);
			V_ion[0] = vel * sin(Vphi) * cos(Vtheta);
			V_ion[1] = vel * sin(Vphi) * sin(Vtheta);
			V_ion[2] = vel * cos(Vphi);
		}
		if (K_test3)
		{
			std::cout << "1V of ion: " << V_ion[0] << ", " << V_ion[1] << ", " << V_ion[2] << " T:"
					  << 1. / 2. * m_A * (pow(V_ion[0], 2) + pow(V_ion[1], 2) + pow(V_ion[2], 2)) / qe << endl;
			std::cout << "1V of neu: " << V_temp[0] << ", " << V_temp[1] << ", " << V_temp[2] << " T:"
					  << 1. / 2. * m_B * (pow(V_temp[0], 2) + pow(V_temp[1], 2) + pow(V_temp[2], 2)) / qe << endl;
		}
		Vtheta = 2 * pi * Random();
		Vphi = pi * Random();
		Modules = sqrt(pow(V_temp[0] - V_ion[0], 2) +
					   pow(V_temp[1] - V_ion[1], 2) + pow(V_temp[2] - V_ion[2], 2));
		// z = 1 for particle B, z = 2 for particle A
		V_1[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] - m_A * Modules * sin(Vphi) * cos(Vtheta));
		V_1[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] - m_A * Modules * sin(Vphi) * sin(Vtheta));
		V_1[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] - m_A * Modules * cos(Vphi));
		V_2[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] + m_B * Modules * sin(Vphi) * cos(Vtheta));
		V_2[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] + m_B * Modules * sin(Vphi) * sin(Vtheta));
		V_2[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] + m_B * Modules * cos(Vphi));
		if (K_test3)
		{
			std::cout << "2V of neu: " << V_1[0] << ", " << V_1[1] << ", " << V_1[2] << " T:"
					  << 1. / 2. * m_B * (pow(V_1[0], 2) + pow(V_1[1], 2) + pow(V_1[2], 2)) / qe << endl;
			std::cout << "2V of ion: " << V_2[0] << ", " << V_2[1] << ", " << V_2[2] << " T:"
					  << 1. / 2. * m_A * (pow(V_2[0], 2) + pow(V_2[1], 2) + pow(V_2[2], 2)) / qe << endl;
		}
		if (z == 1)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_1[i];
		}
		if (z == 2)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_2[i];
		}
		// setfate(0, "rec", "D");
		vel = 1.;
	}
	else if (k == 8)
	{
		Charge_ = 0;
		ChargeTag_ = 0;
		fate_[0] = 3;
		sourcePar_[0] = 2;
		Tn_ = T_N;
		Weight_ = 1.;
		if (this == &D2)
		{
			NumPar_now = Num_D2_pump / (double)numPar_flight_D2;
		}
		if (this == &T2)
		{
			NumPar_now = Num_T2_pump / (double)numPar_flight_T2;
		}

		// location of pumping
		double rand_temp = Random();
		X_[0] = (Grid1.Wall(116, 0) - Grid1.Wall(115, 0)) * rand_temp + Grid1.Wall(115, 0);
		X_[1] = (Grid1.Wall(116, 1) - Grid1.Wall(115, 1)) * rand_temp + Grid1.Wall(115, 1);
		X_[2] = 0;
		XY_[0] = -1;
		XY_[1] = -1;
		Zone_ = 6;

		if (K_Maxwell == 1)
			vel = Maxwell(Tn_, mass_);
		else if (K_Maxwell == 2)
			vel = sqrt((3.0 * qe * Tn_) / mass_);
		Vinit_cosine(V_, cosang[115], sinang[115]);
		// std::cout << "fashe: " << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
		// std::cout << V_[0] << V_[1] << V_[2] << endl;
	}
	else if (k == 9) // Calculation of elastic collision with neutral particle
	{
		/*m_A for Ion, m_B for Neutral, R is the same
		unit vector with a random direction			*/
		double V_ion[3], V_temp[3], m_A, m_B, Modules;
		if (z == 1)
		{
			if (this == &T || this == &T2)
			{
				m_A = Dmass;
				if (K_Vi == 1)
				{
					V_ion[0] = V_D_0_now[0];
					V_ion[1] = V_D_0_now[0];
					V_ion[2] = V_D_0_now[0];
				}
			}
			if (this == &D || this == &D2)
			{
				m_A = Tmass;
				if (K_Vi == 1)
				{
					V_ion[0] = V_T_0_now[0];
					V_ion[1] = V_T_0_now[0];
					V_ion[2] = V_T_0_now[0];
				}
			}
		}
		else if (z == 2)
		{
			if (this == &T || this == &T2)
			{
				m_A = D2mass;
				if (K_Vi == 1)
				{
					V_ion[0] = V_D2_0_now[0];
					V_ion[1] = V_D2_0_now[0];
					V_ion[2] = V_D2_0_now[0];
				}
			}
			if (this == &D || this == &D2)
			{
				m_A = T2mass;
				if (K_Vi == 1)
				{
					V_ion[0] = V_T2_0_now[0];
					V_ion[1] = V_T2_0_now[0];
					V_ion[2] = V_T2_0_now[0];
				}
			}
		}
		m_B = mass_;
		if (Charge_ == 0)
			for (int i = 0; i < 3; i++)
				V_temp[i] = V_[i];
		else
		{
			pause();
			double Rand_temp1, Rand_temp2;
			Rand_temp1 = 2 * pi * Random();
			Rand_temp1 = pi * Random();
			double V_temp1[3], V_temp2[3];
			V_temp1[0] = sin(Rand_temp1) * sin(Rand_temp2);
			V_temp1[1] = sin(Rand_temp1) * cos(Rand_temp2);
			V_temp1[2] = cos(Rand_temp1);
			V_temp2[0] = V_temp1[1] * B[XY_[0]][XY_[1]][2] - V_temp1[2] * B[XY_[0]][XY_[1]][1];
			V_temp2[1] = V_temp1[2] * B[XY_[0]][XY_[1]][0] - V_temp1[0] * B[XY_[0]][XY_[1]][2];
			V_temp2[2] = V_temp1[0] * B[XY_[0]][XY_[1]][1] - V_temp1[1] * B[XY_[0]][XY_[1]][0];
			double Module = sqrt(pow(V_temp2[0], 2) + pow(V_temp2[1], 2) + pow(V_temp2[2], 2));
			V_temp[0] = V_Charge_[0] + V_Charge_[7] * V_temp2[0] / Module; // * erx[XY_[0]][XY_[1]];
			V_temp[1] = V_Charge_[1] + V_Charge_[7] * V_temp2[1] / Module; // * ery[XY_[0]][XY_[1]];
			V_temp[2] = V_Charge_[2] + V_Charge_[7] * V_temp2[2] / Module;
		}
		double Vtheta, Vphi;
		if (K_Vi == 2)
		{
			Vtheta = 2 * pi * Random();
			Vphi = pi * Random();
			if (K_Maxwell == 1)
				vel = Maxwell(Ti[XY_[0]][XY_[1]], m_A);
			else if (K_Maxwell == 2)
				vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / m_A);
			V_ion[0] = vel * sin(Vphi) * cos(Vtheta);
			V_ion[1] = vel * sin(Vphi) * sin(Vtheta);
			V_ion[2] = vel * cos(Vphi);
		}
		if (K_test3)
		{
			std::cout << "1V of ion: " << V_ion[0] << ", " << V_ion[1] << ", " << V_ion[2] << " T:"
					  << 1. / 2. * m_A * (pow(V_ion[0], 2) + pow(V_ion[1], 2) + pow(V_ion[2], 2)) / qe << endl;
			std::cout << "1V of neu: " << V_temp[0] << ", " << V_temp[1] << ", " << V_temp[2] << " T:"
					  << 1. / 2. * m_B * (pow(V_temp[0], 2) + pow(V_temp[1], 2) + pow(V_temp[2], 2)) / qe << endl;
		}
		Vtheta = 2 * pi * Random();
		Vphi = pi * Random();
		Modules = sqrt(pow(V_temp[0] - V_ion[0], 2) + pow(V_temp[1] - V_ion[1], 2) + pow(V_temp[2] - V_ion[2], 2));
		// z = 1 for particle B, z = 2 for particle A
		V_1[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] - m_A * Modules * sin(Vphi) * cos(Vtheta));
		V_1[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] - m_A * Modules * sin(Vphi) * sin(Vtheta));
		V_1[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] - m_A * Modules * cos(Vphi));
		V_2[0] = 1.0 / (m_A + m_B) * (m_A * V_ion[0] + m_B * V_temp[0] + m_B * Modules * sin(Vphi) * cos(Vtheta));
		V_2[1] = 1.0 / (m_A + m_B) * (m_A * V_ion[1] + m_B * V_temp[1] + m_B * Modules * sin(Vphi) * sin(Vtheta));
		V_2[2] = 1.0 / (m_A + m_B) * (m_A * V_ion[2] + m_B * V_temp[2] + m_B * Modules * cos(Vphi));
		if (K_test3)
		{
			std::cout << "2V of neu: " << V_1[0] << ", " << V_1[1] << ", " << V_1[2] << " T:"
					  << 1. / 2. * m_B * (pow(V_1[0], 2) + pow(V_1[1], 2) + pow(V_1[2], 2)) / qe << endl;
			std::cout << "2V of ion: " << V_2[0] << ", " << V_2[1] << ", " << V_2[2] << " T:"
					  << 1. / 2. * m_A * (pow(V_2[0], 2) + pow(V_2[1], 2) + pow(V_2[2], 2)) / qe << endl;
		}
		if (z == 1)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_1[i];
		}
		if (z == 2)
		{
			for (int i = 0; i < 3; i++)
				V_[i] = V_2[i];
		}
		// setfate(0, "rec", "D");
		vel = 1.;
	}
	for (int i = 0; i < 3; i++)
		V_[i] *= vel;
	// std::cout << V_[0] << V_[1] << V_[2] << endl;
}

/// @brief add a collision in the Particle class
/// @param zhonglei 1 for ADAS collision database, 2 for EIRENE collision
/// database
/// @param Coll1 ADAS collision class
/// @param Coll2 EIRENE collision class
void Particle::addCollProb(int zhonglei, ADAS *Coll1, EIRENE *Coll2)
{
	double CollProb[N_poloidal][N_radial];
	if (zhonglei == 1)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				CollProb[i][j] = Coll1->cal(Te[i][j], ne[i][j], Charge_);
			}
	}
	if (zhonglei == 2)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				CollProb[i][j] = Coll2->cal(ne[i][j], Te[i][j]);
			}
	}
}

double Particle::lambda(int i, int j, int charge)
{
	return lambda_[i][j][charge];
}

double Particle::lambda_min(int charge)
{
	return lambda_min_[charge];
}

double Particle::lambda_now()
{
	return lambda_now_;
}

void Particle::CalTn()
{
	/*if (Charge_ == 0)
			Tn_ = 1.0 / 2 * mass_ * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] *
	V_[2]); else Tn_ = 1.0 / 2 * mass_ * (V_Charge_[3] * V_Charge_[3] +
	V_Charge_[7] * V_Charge_[7]);*/
	// std::cout << name_ << Charge_ << '\t' << Tn_ << '\t';
	Tn_ = 1.0 / 3.0 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) / qe;
	// std::cout << Tn_ << '\t' << fate_[0] << '\t' << sourcePar_[0] << endl;
	/*if (Tn_ < 1e-3)
	{
		std::cerr << Tn_ << " velocity is too small"
				  << "\t";
		std::cerr << name_ << '\t' << Charge_ << "+" << '\t' << fate_[0] << endl;
		// exit(3);
	}
	if (Tn_ > 1e3)
	{
		std::cerr << Tn_ << " velocity is too huge"
				  << "\t";
		std::cerr << name_ << '\t' << Charge_ << "+" << '\t' << fate_[0] << endl;
		// exit(3);
	}*/
}

void Particle::VtoVcharge()
{
	V_Charge_[3] = B[XY_[0]][XY_[1]][0] * V_[0] + B[XY_[0]][XY_[1]][1] * V_[1] + B[XY_[0]][XY_[1]][2] * V_[2];
	V_Charge_[0] = V_Charge_[3] * B[XY_[0]][XY_[1]][0];
	V_Charge_[1] = V_Charge_[3] * B[XY_[0]][XY_[1]][1];
	V_Charge_[2] = V_Charge_[3] * B[XY_[0]][XY_[1]][2];
	V_Charge_[4] = V_[0] - V_Charge_[0];
	V_Charge_[5] = V_[1] - V_Charge_[1];
	V_Charge_[6] = V_[2] - V_Charge_[2];
	V_Charge_[7] = sqrt(pow(V_Charge_[4], 2) + pow(V_Charge_[5], 2) + pow(V_Charge_[6], 2));
	ChargeTag_ = 1;
}

void Particle::VchargetoV()
{
	double Rand_temp1, Rand_temp2;
	Rand_temp1 = 2 * pi * Random();
	Rand_temp1 = pi * Random();
	double V_temp1[3], V_temp2[3];
	V_temp1[0] = sin(Rand_temp1) * sin(Rand_temp2);
	V_temp1[1] = sin(Rand_temp1) * cos(Rand_temp2);
	V_temp1[2] = cos(Rand_temp1);
	V_temp2[0] = V_temp1[1] * B[XY_[0]][XY_[1]][2] - V_temp1[2] * B[XY_[0]][XY_[1]][1];
	V_temp2[1] = V_temp1[2] * B[XY_[0]][XY_[1]][0] - V_temp1[0] * B[XY_[0]][XY_[1]][2];
	V_temp2[2] = V_temp1[0] * B[XY_[0]][XY_[1]][1] - V_temp1[1] * B[XY_[0]][XY_[1]][0];
	double Module = sqrt(pow(V_temp2[0], 2) + pow(V_temp2[1], 2) + pow(V_temp2[2], 2));
	V_[0] = V_Charge_[0] + V_Charge_[7] * V_temp2[0] / Module; // * erx[XY_[0]][XY_[1]];
	V_[1] = V_Charge_[1] + V_Charge_[7] * V_temp2[1] / Module; // * ery[XY_[0]][XY_[1]];
	V_[2] = V_Charge_[2] + V_Charge_[7] * V_temp2[2] / Module;
	ChargeTag_ = 0;
}

void Particle::Vchargefix()
{
	V_Charge_[0] = V_Charge_[3] * B[XY_[0]][XY_[1]][0];
	V_Charge_[1] = V_Charge_[3] * B[XY_[0]][XY_[1]][1];
	V_Charge_[2] = V_Charge_[3] * B[XY_[0]][XY_[1]][2];
}

void Particle::track()
{
	if (this == &D)
	{
		if (Charge_ > 0)
		{
			std::cout << fatename(fate_[0]) << ',' << sourcename(sourcePar_[0]) << endl;
			exit(5);
		}
	}
	/*if (Zone_ == 7)
	{
		std::cout << name() << ',' << Charge() << ',' << fatename(fate_[0]) << ',' << sourcename(sourcePar_[0]) << endl;
		std::cout << "X," << X(0) << ',' << X(1) << ';' << X_new(0) << ',' << X_new(1) << endl;
		std::cout << "I don't know where the Particle is." << endl;
		std::cout << "Zone has been 7." << endl;
		exit(0);
	}*/
	if (MeshMode == 1)
	{
		if (Zone_ != 1)
		{
			for (int i = 0; i < 3; i++)
				Xold_[i] = X_new_[i];
		}
		// std::cout << Zone_ << '\t' << Charge_ << "+\t" << endl;
		if (Charge_ == 0)
		{
			if (ChargeTag_ == 1)
			{
				VchargetoV();
				CalTn();
				// std::cout << name_ << '\t' << Charge_ << "+\t" << ChargeTag_ << endl;
			}
			if (K_flight == 1 || K_flight == 3)
			{
				IfColl_ = 0;
				/*if (Zone_ <= 5 && boundary_start_ && IfParticleOut(boundary_start_))
				{
					XY_[0] = XY_new_[0];
					XY_[1] = XY_new_[1];
					Zone_ = Zone_new_;
				}
				boundary_start_ = 0;
				XY_new_[0] = -2;
				XY_new_[1] = -2;*/
				d_flight_ = 1.;
				Rand_flight_ = Random();
				while (!IfColl_ && Zone_ < 6 && Zone_ > 1)
				{
					/*if (boundary_start_ && IfParticleOut(boundary_start_))
					{
						XY_[0] = XY_new_[0];
						XY_[1] = XY_new_[1];
						Zone_ = Zone_new_;
						if (Zone_ > 5)
							break;
					}
					boundary_start_ = 0;
					XY_new_[0] = -2;
					XY_new_[1] = -2;*/
					// dt_ = -log(Random()) * lambda_min_[XY_[0]][XY_[1]][Charge_];
					/*if (boundary_start_)
					{
						if (XY_new_[0] == 0)
						{
							std::cout << "why?" << endl;
						}
						XY_[0] = XY_new_[0];
						XY_[1] = XY_new_[1];
						Zone_ = Zone_new_;
						boundary_start_ = 0;
						XY_new_[0] = -2;
						XY_new_[1] = -2;
					}*/
					if (K_Tn == 1)
					{
						CalLambda();
						if (K_flight == 1)
						{
							dt_ = lambda_now_;
						}
						else if (K_flight == 3)
						{
							dt_ = -log(Random()) * lambda_min_[Charge_];
						}
						/*if (dt_ <= 0)
						{
							std::cout << name_ + " dt = 0" << endl;
							// pause();
							exit(0);
						}*/
					}
					else if (K_Tn == 2)
					{
						if (K_flight == 1)
						{
							dt_ = lambda_[XY_[0]][XY_[1]][Charge_];
						}
						else if (K_flight == 3)
						{
							dt_ = -log(Random()) * lambda_min_[Charge_];
						}
					}
					if (XY_[0] == -1 || XY_[1] == -1)
					{
						std::cout << "dt_ was calculated error!" << endl;
						std::cout << name() << ',' << Charge() << ',' << fatename(fate_[0]) << ',' << sourcename(sourcePar_[0]) << endl;
						std::cout << "X," << X(0) << ',' << X(1) << ';' << X_new(0) << ',' << X_new(1) << endl;
						// pause();
						exit(0);
					}
					for (int i = 0; i < 3; i++)
						X_new_[i] = X_[i] + V_[i] * dt_;
					// std::cout << "dt=" << dt_ << endl;
					Caltrace();
					/*if (!Grid1.IfinIndex(XY_[0], XY_[1], X_new_[0], X_new_[1]))
					{
						std::cout << name_ + ": " << X_[0] << "," << X_[1] << " " << X_new_[0] << "," << X_new_[1] << " " << XY_[0] << "," << XY_[1] << endl;
						std::cout << "the Caltrace() in Particle.cpp have some problem.\n";
					}*/
				}
				d_flight_ = 1.;
				if (K_test1)
				{
					std::cout << "see see what happened: " << endl;
					std::cout << "X," << X(0) << ',' << X(1) << ';' << X_new(0) << ',' << X_new(1) << endl;
				}
				if (Zone() > 5)
				{
					dt_ = 1000;
					double V_temp[3];
					V_temp[0] = V_[0] / sqrt(V_[0] * V_[0] + V_[1] * V_[1]);
					V_temp[1] = V_[1] / sqrt(V_[0] * V_[0] + V_[1] * V_[1]);
					for (int i = 0; i < 2; i++)
						X_new_[i] = X_[i] + dt_ * V_temp[i];

					if (K_test1)
					{
						std::cout << name_ << " X: " << X_[0] << ',' << X_[1] << '\t' << X_new_[0] << ',' << X_new_[1] << endl;
						std::cout << "V_temp: " << V_temp[0] << ',' << V_temp[1] << endl;
						std::cout << endl;
					}
				}
				/*else if (Zone() == 7)
				{
					dt_ = 1e-8;
					for (int i = 0; i < 3; i++)
						X_new_[i] = X_[i] + dt_ * V_[i];
				}*/
			}
			else if (K_flight == 2)
			{
				for (int i = 0; i < 3; i++)
					X_new_[i] = X_[i] + dt * V_[i];
				if (K_GRID)
					Grid1.Find(X_new_, &Zone_, XY_);
				else
					Find();
			}
			/*else if (K_flight == 3)
			{
				dt_ = -log(Random()) * lambda_min_[Charge_];
				for (int i = 0; i < 3; i++)
					X_new_[i] = X_[i] + dt_ * V_[i];
				if (K_GRID)
					Grid1.Find(X_new_, &Zone_, XY_);
				else
					Find();
			}*/
			if (Zone_ == 1)
			{
				int num_intersect = 0;
				if (IfFlightOut_ == 1)
				{
					num_intersect = 1;
					IfFlightOut_ = 0;
				}
				// std::cout << PP->X(0) << '\t' << PP->X(1) << '\t' << PP->X_new(0) << '\t' << PP->X_new(1) << endl;
				/*for (int i = 0; i < Grid1.Core_Boundry_num(); i++)
				{
					if (get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid1.Core_Boundry(i, 0), Grid1.Core_Boundry(i, 1), Grid1.Core_Boundry(i + 1, 0), Grid1.Core_Boundry(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
					{
						InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - X_[0]), 2) + pow((InterscePoint[num_intersect][2] - X_[1]), 2));
						InterscePoint[num_intersect][3] = i;
						InterscePoint[num_intersect++][4] = 18;
					}
				}
				if (num_intersect > 1)
				{
					for (int i = 1; i < num_intersect - 1; i++)
					{
						if (InterscePoint[i][0] < InterscePoint[0][0])
							for (int j = 0; j < 5; j++)
								InterscePoint[0][j] = InterscePoint[i][j];
					}
				}*/
				if (num_intersect < 1)
				{
					std::cout << name() << ',' << Charge() << ',' << fatename(fate_[0]) << ',' << sourcename(sourcePar_[0]) << endl;
					std::cout << "X," << X(0) << ',' << X(1) << ';' << X_new(0) << ',' << X_new(1) << endl;
					std::cout << "XY: " << XY(0) << ',' << XY(1) << '\t' << Zone_ << endl;
					std::cout << "I don't know where the Particle after fly in the Core." << endl;
					exit(0);
				}
				// Flux_CoreBoundry[(int)InterscePoint[0][3]][0] += Weight_ * NumPar_now;
				Weight_ = 0.;
				/*for (int i = 0; i < 3; i++)
				{
					X_new_[i] = Xold_[i];
				}
				Zone_ = Zone_old_;

				if (this == &H || this == &D || this == &T)
				{
					Ion_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
					Ion_[Charge_].Mu_Add(XY_, 0.);
					Ion_[Charge_].E_Add(XY_, Tn_ * qe * Weight_ * NumPar_now);
					Weight_ = 0;
					setfate(1, 5, 18);
					// std::cerr << "Particle " << name_ << " has gone in core!" << endl;
					return;
				}
				else if (this == &H2 || this == &D2 || this == &T2)
				{
					Ion_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					// Init(); // electron collision
					Charge_ = 1;
					if (!K_D2Flight)
					{
						Sum_n_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now;
						Weight_ = 0;
						return;
					}
				}
				else
				{
					Charge_ += 1;
					track();
				}*/
			}

			/// if triangle is off
			/*while (Zone_ == 6)
			{
				if (K_flight == 1)
				{
					dt_ = log(Random()) * lambda_[XY_[0]][XY_[1]][Charge_];
					for (int i = 0; i < 3; i++)
						X_new_[i] = X_[i] - V_[i] * dt_ * 100.;
				}
				else if (K_flight == 2)
				{
					for (int i = 0; i < 3; i++)
						X_new_[i] = X_[i] + dt * V_[i] * 100.;
				}
			}*/
		}
		else if (Charge_ > 0)
		{
			if (ChargeTag_ == 0)
			{
				VtoVcharge();
			}
			divimp_F();
			divimp_Ti();
			//  std::cout << name_ << '\t' << V_[0] << '\t' << V_[1] << '\t' << V_[2] <<
			//  endl; std::cout << Rand_temp << '\t' << sin(Rand_temp) << '\t' <<
			//  cos(Rand_temp) << endl;
			V_Charge_[0] = V_Charge_[3] * B[XY_[0]][XY_[1]][0];
			V_Charge_[1] = V_Charge_[3] * B[XY_[0]][XY_[1]][1];
			V_Charge_[2] = V_Charge_[3] * B[XY_[0]][XY_[1]][2];
			for (int i = 0; i < 3; i++)
				X_new_[i] = X_[i] + dt * V_Charge_[i];
			divimp_E();
			if (K_EcrossBDrift)
				divimp_drift_E();
			if (K_abnormal_transport)
				if (XY_[0] > 1 && XY_[0] < 96)
					divimp_anomalous_diffusion();
			if (K_GRID)
				Grid1.Find(X_new_, &Zone_, XY_);
			else
				Find();
			/*if (Zone_ == 1)
			{
				for (int i = 0; i < 3; i++)
				{
					X_[i] = Xold_[i];
				}
				track();
			}
			if (Zone_ == 6)
				if (Random() < 0.95 && backGridBoundry)
				{
					for (int i = 0; i < 3; i++)
					{
						X_[i] = Xold_[i];
					}
					// Charge_ -= 1;
					setfate(0, 2, 19);
					track();
				}
				else
				{
					Charge_ -= 1;
					track();
				}*/
			// std::cout << V_Charge_[0] << '\t' << V_Charge_[1] << '\t' << V_Charge_[2]
			// << '\t' << V_Charge_[3] << '\t' << V_Charge_[4] << endl;
			// std::cout << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
		}
		// FluxCal_Grid();
	}
	else if (MeshMode == 3)
	{
		if (Zone_ != 1)
		{
			for (int i = 0; i < 3; i++)
				Xold_[i] = X_new_[i];
		}
		if (Charge_ == 0)
		{
			if (ChargeTag_ == 1)
			{
				VchargetoV();
				CalTn();
				// std::cout << name_ << '\t' << Charge_ << "+\t" << ChargeTag_ << endl;
			}
			if (K_flight == 1 || K_flight == 3)
			{
				IfColl_ = 0;
				d_flight_ = 1.;
				Rand_flight_ = Random();
				while (!IfColl_ && Zone_ < 7)
				{
					if (K_Tn == 1)
					{
						CalLambda();
						if (K_flight == 1)
						{
							dt_ = lambda_now_;
						}
						else if (K_flight == 3)
						{
							dt_ = -log(Random()) * lambda_min_[Charge_];
						}
					}
					else if (K_Tn == 2)
					{
						if (K_flight == 1)
						{
							dt_ = lambda_[XY_[0]][XY_[1]][Charge_];
						}
						else if (K_flight == 3)
						{
							dt_ = -log(Random()) * lambda_min_[Charge_];
						}
					}
					for (int i = 0; i < 3; i++)
						X_new_[i] = X_[i] + V_[i] * dt_;
					Caltrace_Tri();
					if (Zone_ == 1)
					{
						Weight_ = 0.;
						return;
					}
					if (IfColl_ == 1)
					{
						return;
					}
					if (Zone_ == 7)
					{
						if (InterscePoint[0][4] == 11) // wall
						{
							return;
						}
						else if (InterscePoint[0][4] == 1) // target
						{
							return;
						}
						else
						{
							std::cout << "InterscePoint: " << InterscePoint[0][4] << endl;
							pause();
						}
					}
				}
			}
		}
	}
}

void Particle::CalLambda()
{
	double CS_Vacuum = 0.1;
	if (Zone_ < 6)
	{
		if (K_Vi == 1)
		{
			double Vtheta, Vphi, vel;
			if (Zone_ < 6)
			{
				if (K_H)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(Ti[XY_[0]][XY_[1]], Hmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / Hmass);
					V_H_1_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_H_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][0];
					V_H_1_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_H_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][1];
					V_H_1_now[2] = vel * cos(Vphi) + ua_H_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][2];
				}

				if (K_D)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(Ti[XY_[0]][XY_[1]], Dmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / Dmass);
					V_D_1_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_D_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][0];
					V_D_1_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_D_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][1];
					V_D_1_now[2] = vel * cos(Vphi) + ua_D_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][2];

					// std::cout << V_D_1_now[0] << '\t' << V_D_1_now[1] << '\t' << V_D_1_now[2] << endl;
					/// statistic the velocity of
					/*for (int i = 0; i < 3; i++)
						Sum_V_D_1_[XY_[0]][XY_[1]][i][0] += V_D_1_now[i];
					Num_V_D_1_[XY_[0]][XY_[1]][0] += 1;*/
				}

				if (K_T)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					if (K_Maxwell == 1)
						vel = Maxwell(Ti[XY_[0]][XY_[1]], Tmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * Ti[XY_[0]][XY_[1]]) / Tmass);
					V_T_1_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_T_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][0];
					V_T_1_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_T_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][1];
					V_T_1_now[2] = vel * cos(Vphi) + ua_T_1[XY_[0]][XY_[1]] * B[XY_[0]][XY_[1]][2];
				}
			}
			if (MeshMode == 3)
			{
				if (K_H)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_H_0_Tri[Tri_Index_], Hmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_H_0_Tri[Tri_Index_]) / Hmass);
					V_H_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_H_0_Tri[Tri_Index_];
					V_H_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_H_0_Tri[Tri_Index_];
					V_H_0_now[2] = vel * cos(Vphi) + ua_H_0_Tri[Tri_Index_];

					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_H2_0_Tri[Tri_Index_], H2mass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_H2_0_Tri[Tri_Index_]) / H2mass);
					V_H2_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_H2_0_Tri[Tri_Index_];
					V_H2_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_H2_0_Tri[Tri_Index_];
					V_H2_0_now[2] = vel * cos(Vphi) + ua_H2_0_Tri[Tri_Index_];
				}
				if (K_D)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_D_0_Tri[Tri_Index_], Dmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_D_0_Tri[Tri_Index_]) / Dmass);
					V_D_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_D_0_Tri[Tri_Index_];
					V_D_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_D_0_Tri[Tri_Index_];
					V_D_0_now[2] = vel * cos(Vphi) + ua_D_0_Tri[Tri_Index_];

					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_D2_0_Tri[Tri_Index_], D2mass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_D2_0_Tri[Tri_Index_]) / D2mass);
					V_D2_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_D2_0_Tri[Tri_Index_];
					V_D2_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_D2_0_Tri[Tri_Index_];
					V_D2_0_now[2] = vel * cos(Vphi) + ua_D2_0_Tri[Tri_Index_];
				}
				if (K_T)
				{
					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_T_0_Tri[Tri_Index_], Tmass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_T_0_Tri[Tri_Index_]) / Tmass);
					V_T_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_T_0_Tri[Tri_Index_];
					V_T_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_T_0_Tri[Tri_Index_];
					V_T_0_now[2] = vel * cos(Vphi) + ua_T_0_Tri[Tri_Index_];

					Vtheta = 2 * pi * Random();
					Vphi = pi * Random();
					// std::cout << Vtheta << '\t' << Vphi << '\t' << Vtheta / Vphi << endl;
					if (K_Maxwell == 1)
						vel = Maxwell(T_T2_0_Tri[Tri_Index_], T2mass);
					else if (K_Maxwell == 2)
						vel = sqrt((3.0 * qe * T_T2_0_Tri[Tri_Index_]) / T2mass);
					V_T2_0_now[0] = vel * sin(Vphi) * cos(Vtheta) + ua_T2_0_Tri[Tri_Index_];
					V_T2_0_now[1] = vel * sin(Vphi) * sin(Vtheta) + ua_T2_0_Tri[Tri_Index_];
					V_T2_0_now[2] = vel * cos(Vphi) + ua_T2_0_Tri[Tri_Index_];
				}
			}
		}
	}

	if (MeshMode == 1)
	{
		if (this == &H || this == &D || this == &T)
		{
			if (this == &H)
			{
				CX_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);

				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					CX_[0].Setcs_now(CX_[0].cs(XY_[0], XY_[1]));
					CX_DT_[0].Setcs_now(0);
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now());
				}
				else if (K_database == 2)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * 0.5 * Hmass * CX_[0].V_relative_2() / qe, Ti[XY_[0]][XY_[1]]));
					CX_DT_[0].Setcs_now(0);
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now());
				}
			}
			else if (this == &D)
			{
				CX_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);

				if (K_DT)
					CX_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);

				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					CX_[0].Setcs_now(CX_[0].cs(XY_[0], XY_[1]));
					if (K_CX_DT)
						CX_DT_[0].Setcs_now(CX_DT_[0].cs(XY_[0], XY_[1]));
					else
						CX_DT_[0].Setcs_now(0.);
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now());
				}
				else if (K_database == 2)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					if (K_Vi == 1)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * (0.5 * Dmass) * CX_[0].V_relative_2() / qe,
																			   Ti[XY_[0]][XY_[1]]));
						if (K_CX_DT)
							CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * DTmass * CX_DT_[0].V_relative_2() / qe,
																					  Ti[XY_[0]][XY_[1]]));
						else
							CX_DT_[0].Setcs_now(0.);
					}
					if (K_Vi == 2)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
						if (K_CX_DT)
							CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
						else
							CX_DT_[0].Setcs_now(0.);
					}
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now());
				}
			}
			else if (this == &T)
			{
				CX_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
					CX_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);

				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					CX_[0].Setcs_now(CX_[0].cs(XY_[0], XY_[1]));
					if (K_CX_DT)
						CX_DT_[0].Setcs_now(CX_DT_[0].cs(XY_[0], XY_[1]));
					else
						CX_DT_[0].Setcs_now(0.);
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now());
				}
				else if (K_database == 2)
				{
					if (K_Vi == 1)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * (0.5 * Tmass) * CX_[0].V_relative_2() / qe,
																			   Ti[XY_[0]][XY_[1]]));
						if (K_CX_DT)
							CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * DTmass * CX_DT_[0].V_relative_2() / qe,
																					  Ti[XY_[0]][XY_[1]]));
						else
							CX_DT_[0].Setcs_now(0.);
					}
					Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
					if (K_Vi == 2)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
						if (K_CX_DT)
							CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
						else
							CX_DT_[0].Setcs_now(0.);
					}
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now());
				}
			}
		}
		else if (this == &H2 || this == &D2 || this == &T2)
		{
			if (this == &H2)
			{
				CX_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);

				Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
				Diss1_[0].Setcs_now(Diss1_[0].cs(XY_[0], XY_[1]));
				Diss2_[0].Setcs_now(Diss2_[0].cs(XY_[0], XY_[1]));
				CX_[0].Setcs_now(n_H_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * HH2mass * (pow(V_[0] - V_H_1_now[0], 2) + pow(V_[1] - V_H_1_now[1], 2) + pow(V_[2] - V_H_1_now[2], 2)) / qe, Ti[XY_[0]][XY_[1]]));
				Ela_[0].Setcs_now(n_H_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(1. / 2. * HH2mass * (pow(V_[0] - V_H_1_now[0], 2) + pow(V_[1] - V_H_1_now[1], 2) + pow(V_[2] - V_H_1_now[2], 2)) / qe, Ti[XY_[0]][XY_[1]]));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() + CX_[0].cs_now() + Ela_[0].cs_now());
			}
			else if (this == &D2)
			{
				CX_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
				{
					CX_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
					Ela_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				}

				Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
				Diss1_[0].Setcs_now(Diss1_[0].cs(XY_[0], XY_[1]));
				Diss2_[0].Setcs_now(Diss2_[0].cs(XY_[0], XY_[1]));
				if (K_Vi == 1)
				{
					CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * DD2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																		   Ti[XY_[0]][XY_[1]]));
					Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * DD2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																		   Ti[XY_[0]][XY_[1]]));
					if (K_DT)
					{
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * D2Tmass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]]));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * D2Tmass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]]));
					}
				}
				if (K_Vi == 2)
				{
					CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
					// Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
					// Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
					Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
					if (K_DT)
					{
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
					}
				}
				DS_[1][0].Setcs_now(DS_[1][0].cs(XY_[0], XY_[1]));
				DS_[1][1].Setcs_now(DS_[1][1].cs(XY_[0], XY_[1]));
				DS_[1][2].Setcs_now(DS_[1][2].cs(XY_[0], XY_[1]));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() + CX_[0].cs_now() + Ela_[0].cs_now() + CX_DT_[0].cs_now() + Ela_DT_[0].cs_now());
			}
			else if (this == &T2)
			{
				CX_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
				{
					CX_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
					Ela_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				}

				Ion_[0].Setcs_now(Ion_[0].cs(XY_[0], XY_[1]));
				Diss1_[0].Setcs_now(Diss1_[0].cs(XY_[0], XY_[1]));
				Diss2_[0].Setcs_now(Diss2_[0].cs(XY_[0], XY_[1]));
				if (K_Vi == 1)
				{
					CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * TT2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																		   Ti[XY_[0]][XY_[1]]));
					Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * TT2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																		   Ti[XY_[0]][XY_[1]]));
					if (K_DT)
					{
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * DT2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]]));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * DT2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]]));
					}
				}
				if (K_Vi == 2)
				{
					CX_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
					Ela_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
					if (K_DT)
					{
						CX_DT_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
						Ela_DT_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
					}
				}
				DS_[1][0].Setcs_now(DS_[1][0].cs(XY_[0], XY_[1]));
				DS_[1][1].Setcs_now(DS_[1][1].cs(XY_[0], XY_[1]));
				DS_[1][2].Setcs_now(DS_[1][2].cs(XY_[0], XY_[1]));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() + CX_[0].cs_now() + Ela_[0].cs_now() + CX_DT_[0].cs_now() + Ela_DT_[0].cs_now());
			}
		}
	}
	else if (MeshMode == 3)
	{
		if (this == &H || this == &D || this == &T)
		{
			if (this == &H)
			{
				CX_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_H_0_now[0], V_H_0_now[1], V_H_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_H2_0_now[0], V_H2_0_now[1], V_H2_0_now[2], V_[0], V_[1], V_[2]);

				R_with_H_[0].Setcs_now(R_with_H_[0].cs(Tri_Index_));
				R_with_H2_[0].Setcs_now(R_with_H2_[0].cs(Tri_Index_));
				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
					CX_[0].Setcs_now(CX_[0].cs(Tri_Index_));
					CX_DT_[0].Setcs_now(CS_Vacuum);
				}
				else if (K_database == 2)
				{
					if (Zone_ < 6)
					{
						Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * 0.5 * Hmass * (pow(V_[0] - V_H_1_now[0], 2) + pow(V_[1] - V_H_1_now[1], 2) + pow(V_[2] - V_H_1_now[2], 2)) / qe, Ti[XY_[0]][XY_[1]]));
						CX_DT_[0].Setcs_now(CS_Vacuum);
					}
					else
					{
						Ion_[0].Setcs_now(CS_Vacuum);
						CX_[0].Setcs_now(CS_Vacuum);
						CX_DT_[0].Setcs_now(CS_Vacuum);
					}
				}
				lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H2_[0].cs_now());
			}
			else if (this == &D)
			{
				CX_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
					CX_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_D_0_now[0], V_D_0_now[1], V_D_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_D2_0_now[0], V_D2_0_now[1], V_D2_0_now[2], V_[0], V_[1], V_[2]);

				R_with_H_[0].Setcs_now(R_with_H_[0].cs(Tri_Index_));
				R_with_H2_[0].Setcs_now(R_with_H2_[0].cs(Tri_Index_));

				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
					CX_[0].Setcs_now(CX_[0].cs(Tri_Index_));
					if (K_CX_DT)
						CX_DT_[0].Setcs_now(CX_DT_[0].cs(Tri_Index_));
					else
						CX_DT_[0].Setcs_now(CS_Vacuum);
				}
				else if (K_database == 2)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
					if (K_Vi == 1)
					{
						if (Zone_ < 6)
						{
							CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * 0.5 * Dmass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																				   Ti[XY_[0]][XY_[1]]));
							if (K_CX_DT)
								CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * DTmass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																						  Ti[XY_[0]][XY_[1]]));
							else
								CX_DT_[0].Setcs_now(CS_Vacuum);
						}
						else
						{
							CX_[0].Setcs_now(CS_Vacuum);
							CX_DT_[0].Setcs_now(CS_Vacuum);
						}
					}
					if (K_Vi == 2)
					{
						if (Zone_ < 6)
						{
							CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
							if (K_CX_DT)
								CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
							else
								CX_DT_[0].Setcs_now(CS_Vacuum);
						}
						else
						{
							CX_[0].Setcs_now(CS_Vacuum);
							CX_DT_[0].Setcs_now(CS_Vacuum);
						}
					}
				}
				lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H2_[0].cs_now());
			}
			else if (this == &T)
			{
				CX_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
					CX_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_T_0_now[0], V_T_0_now[1], V_T_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_T2_0_now[0], V_T2_0_now[1], V_T2_0_now[2], V_[0], V_[1], V_[2]);

				R_with_H_[0].Setcs_now(R_with_H_[0].cs(Tri_Index_));
				R_with_H2_[0].Setcs_now(R_with_H2_[0].cs(Tri_Index_));
				if (K_database == 1)
				{
					Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
					CX_[0].Setcs_now(CX_[0].cs(Tri_Index_));
					if (K_CX_DT)
						CX_DT_[0].Setcs_now(CX_DT_[0].cs(Tri_Index_));
					else
						CX_DT_[0].Setcs_now(CS_Vacuum);
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H_[0].cs_now());
				}
				else if (K_database == 2)
				{
					if (Zone_ < 6)
					{
						if (K_Vi == 1)
						{
							CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * Tmass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																				   Ti[XY_[0]][XY_[1]] / coefficient_T));
							if (K_CX_DT)
								CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(1. / 2. * Tmass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																						  Ti[XY_[0]][XY_[1]] / coefficient_D));
							else
								CX_DT_[0].Setcs_now(CS_Vacuum);
						}
						else if (K_Vi == 2)
						{
							CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
							if (K_CX_DT)
								CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_1_8_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
							else
								CX_DT_[0].Setcs_now(CS_Vacuum);
						}
					}
					else
					{
						CX_[0].Setcs_now(CS_Vacuum);
						CX_DT_[0].Setcs_now(CS_Vacuum);
					}
					lambda_now_ = 1. / (Ion_[0].cs_now() + CX_[0].cs_now() + CX_DT_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H2_[0].cs_now());
				}
			}
		}
		else if (this == &H2 || this == &D2 || this == &T2)
		{
			if (this == &H2)
			{
				CX_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_H_1_now[0], V_H_1_now[1], V_H_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_H_0_now[0], V_H_0_now[1], V_H_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_H2_0_now[0], V_H2_0_now[1], V_H2_0_now[2], V_[0], V_[1], V_[2]);

				Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
				Diss1_[0].Setcs_now(Diss1_[0].cs(Tri_Index_));
				Diss2_[0].Setcs_now(Diss2_[0].cs(Tri_Index_));
				CX_[0].Setcs_now(n_H_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_, Ti[XY_[0]][XY_[1]]));
				Ela_[0].Setcs_now(n_H_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(3. / 2. * Tn_, Ti[XY_[0]][XY_[1]]));
				R_with_H_[0].Setcs_now(R_with_H_[0].cs(Tri_Index_));
				R_with_H2_[0].Setcs_now(R_with_H2_[0].cs(Tri_Index_));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() + CX_[0].cs_now() + Ela_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H_[0].cs_now());
			}
			else if (this == &D2)
			{
				CX_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_D_0_now[0], V_D_0_now[1], V_D_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_D2_0_now[0], V_D2_0_now[1], V_D2_0_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
				{
					CX_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
					Ela_DT_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				}

				Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
				Diss1_[0].Setcs_now(Diss1_[0].cs(Tri_Index_));
				Diss2_[0].Setcs_now(Diss2_[0].cs(Tri_Index_));
				if (Zone_ < 6)
				{
					if (K_Vi == 1)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																			   Ti[XY_[0]][XY_[1]] / coefficient_D));
						Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																			   Ti[XY_[0]][XY_[1]] / coefficient_D));
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]] / coefficient_T));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]] / coefficient_T));
					}
					if (K_Vi == 2)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
						Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_D));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_D, Ti[XY_[0]][XY_[1]] / coefficient_T));
					}
				}
				else
				{
					CX_[0].Setcs_now(CS_Vacuum);
					CX_DT_[0].Setcs_now(CS_Vacuum);
					Ela_[0].Setcs_now(CS_Vacuum);
					Ela_DT_[0].Setcs_now(CS_Vacuum);
				}
				DS_[1][0].Setcs_now(DS_[1][0].cs(Tri_Index_));
				DS_[1][1].Setcs_now(DS_[1][1].cs(Tri_Index_));
				DS_[1][2].Setcs_now(DS_[1][2].cs(Tri_Index_));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() +
									CX_[0].cs_now() + Ela_[0].cs_now() + CX_DT_[0].cs_now() + Ela_DT_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H2_[0].cs_now());
			}
			else if (this == &T2)
			{
				CX_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				Ela_[0].Set_V_relative(V_T_1_now[0], V_T_1_now[1], V_T_1_now[2], V_[0], V_[1], V_[2]);
				R_with_H_[0].Set_V_relative(V_T_0_now[0], V_T_0_now[1], V_T_0_now[2], V_[0], V_[1], V_[2]);
				R_with_H2_[0].Set_V_relative(V_T2_0_now[0], V_T2_0_now[1], V_T2_0_now[2], V_[0], V_[1], V_[2]);
				if (K_DT)
				{
					CX_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
					Ela_DT_[0].Set_V_relative(V_D_1_now[0], V_D_1_now[1], V_D_1_now[2], V_[0], V_[1], V_[2]);
				}

				Ion_[0].Setcs_now(Ion_[0].cs(Tri_Index_));
				Diss1_[0].Setcs_now(Diss1_[0].cs(Tri_Index_));
				Diss2_[0].Setcs_now(Diss2_[0].cs(Tri_Index_));
				if (Zone_ < 6)
				{
					if (K_Vi == 1)
					{
						CX_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																			   Ti[XY_[0]][XY_[1]] / coefficient_T));
						Ela_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_T_1_now[0], 2) + pow(V_[1] - V_T_1_now[1], 2) + pow(V_[2] - V_T_1_now[2], 2)) / qe,
																			   Ti[XY_[0]][XY_[1]] / coefficient_T));
						CX_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]] / coefficient_D));
						Ela_DT_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(1. / 2. * D2mass * (pow(V_[0] - V_D_1_now[0], 2) + pow(V_[1] - V_D_1_now[1], 2) + pow(V_[2] - V_D_1_now[2], 2)) / qe,
																				  Ti[XY_[0]][XY_[1]] / coefficient_D));
					}
					if (K_Vi == 2)
					{
						CX_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
						CX_DT_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R3_2_3_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
						// Ela_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
						// Ela_DT_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3T_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
						Ela_[0].Setcs_now(n_T_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_T));
						Ela_DT_[0].Setcs_now(n_D_1[XY_[0]][XY_[1]] * R0_3D_H3.cal(3. / 2. * Tn_ / coefficient_T, Ti[XY_[0]][XY_[1]] / coefficient_D));
					}
				}
				else
				{
					CX_[0].Setcs_now(CS_Vacuum);
					CX_DT_[0].Setcs_now(CS_Vacuum);
					Ela_[0].Setcs_now(CS_Vacuum);
					Ela_DT_[0].Setcs_now(CS_Vacuum);
				}
				DS_[1][0].Setcs_now(DS_[1][0].cs(Tri_Index_));
				DS_[1][1].Setcs_now(DS_[1][1].cs(Tri_Index_));
				DS_[1][2].Setcs_now(DS_[1][2].cs(Tri_Index_));
				lambda_now_ = 1. / (Ion_[0].cs_now() + Diss1_[0].cs_now() + Diss2_[0].cs_now() +
									CX_[0].cs_now() + Ela_[0].cs_now() + CX_DT_[0].cs_now() + Ela_DT_[0].cs_now() + R_with_H_[0].cs_now() + R_with_H_[0].cs_now());
			}
		}
	}
}

/// @brief calculate the trace of the Particle
/// @param tag 0 for collision in this grid, 1 for no collision in the grid
void Particle::Caltrace()
{
	/*if (!Grid1.IfinIndex1(XY_[0], XY_[1], X_[0], X_[1]))
	{
		std::cout << name_ + ": " << X_[0] << "," << X_[1] << " " << X_new_[0] << "," << X_new_[1] << " " << XY_[0] << "," << XY_[1] << endl;
		std::cout << "the X is not in the XY in Caltrace() of Particle.cpp have some problem.";
		std::cout << endl;
	}*/
	int FindIt = 0;
	boundary_start_ = 0;
	IfColl_ = 1;
	double secx, secy, distance;
	// judge if the particle cross through the first boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 1), &secx, &secy))
	{
		if (!PointandLine(Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 1), X_new_[0], X_new_[1]))
		{
			FindIt = 1;
			boundary_start_ = 1;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_)
				{
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 1 in Caltrace()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}
				}
			}

			NumParStat();
			n_Flux_Grid_[XY_[0]][XY_[1]][0] += Weight_ * NumPar_now;
			T_Flux_Grid_[XY_[0]][XY_[1]][0] += Tn_ * Weight_ * NumPar_now;
			if (XY_[1] == 1)
			{
				if (XY_[0] < 25 || XY_[0] > 72)
				{
					if (XY_[0] < 12)
						FT_.accumulate(0, 8, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					else if (XY_[0] < 25)
						FT_.accumulate(1, 8, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] > 84)
						FT_.accumulate(7, 8, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					else if (XY_[0] > 75)
						FT_.accumulate(6, 8, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = -1;
					XY_[1] = -1;
					Zone_ = 6;
					IfColl_ = 0;
					for (int i = 0; i < Grid1.PLasma_Grid_Boundry_num(); i++)
					{
						if (X_[0] >= min(Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i + 1, 0)) && X_[0] <= max(Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i + 1, 0)) && X_[1] >= min(Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 1)) && X_[1] <= max(Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 1)))
						{
							AddPlasmaBoundaryEro(i);
							break;
						}
					}
				}
				else
				{
					IfFlightOut_ = 1;
					InterscePoint[0][0] = 0;
					InterscePoint[0][1] = X_[0];
					InterscePoint[0][2] = X_[1];
					for (int i = 0; i < Grid1.Core_Boundry_num(); i++)
					{
						if (X_[0] >= min(Grid1.Core_Boundry(i, 0), Grid1.Core_Boundry(i + 1, 0)) && X_[0] <= max(Grid1.Core_Boundry(i, 0), Grid1.Core_Boundry(i + 1, 0)) && X_[1] >= min(Grid1.Core_Boundry(i, 1), Grid1.Core_Boundry(i + 1, 1)) && X_[1] <= max(Grid1.Core_Boundry(i, 1), Grid1.Core_Boundry(i + 1, 1)))
						{
							InterscePoint[0][3] = i;
							break;
						}
					}
					InterscePoint[0][4] = 18;

					boundary_start_ = 0;
					XY_[0] = -1;
					XY_[1] = -1;
					Zone_ = 1;
					IfColl_ = 0;
				}
			}
			else
			{
				if (XY_[1] == 19)
				{
					if (XY_[0] < 25)
					{
						Zone_ = 4;
					}
					else if (XY_[0] > 72)
					{
						Zone_ = 5;
					}
					else
					{
						Zone_ = 2;
					}
				}
				XY_[0] = XY_[0];
				XY_[1] = XY_[1] - 1;
			}
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (K_test1)
			{
				std::cout << FindIt << "(after): " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
				std::cout << "boundary_start: " << boundary_start_ << endl;
			}
			return;
		}
	}

	// judge if the particle cross through the second boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 1), &secx, &secy))
	{
		if (!PointandLine(Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(1, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 1), X_new_[0], X_new_[1]))
		{
			FindIt = 2;
			if (K_test1)
			{
				std::cout << FindIt << ": " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
			}
			boundary_start_ = 2;
			IfColl_ = 0;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_)
				{
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 2 in Caltrace()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}
				}
			}
			NumParStat();
			n_Flux_Grid_[XY_[0]][XY_[1]][1] += Weight_ * NumPar_now;
			T_Flux_Grid_[XY_[0]][XY_[1]][1] += Tn_ * Weight_ * NumPar_now;

			if (XY_[0] == 96)
			{
				IfFlightOut_ = 7;
				boundary_start_ = 0;
				IfColl_ = 0;

				InterscePoint[0][0] = 0;
				InterscePoint[0][1] = secx;
				InterscePoint[0][2] = secy;
				InterscePoint[0][3] = XY_[1] + N_radial;
				InterscePoint[0][4] = 1;

				XY_[0] = -1;
				XY_[1] = -1;
				Zone_ = 7;
			}
			else
			{
				if (XY_[1] < 19 && XY_[0] == 72)
				{
					FT_.accumulate(5, 2, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = 25;
				}
				else if (XY_[1] < 19 && XY_[0] == 24)
				{
					FT_.accumulate(1, 6, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = 73;
					Zone_ = 5;
				}
				else
				{
					if (XY_[0] == 12)
						FT_.accumulate(0, 1, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[1] >= 19 && XY_[0] == 24)
						FT_.accumulate(1, 2, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 32)
						FT_.accumulate(2, 3, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 48)
						FT_.accumulate(3, 4, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 64)
						FT_.accumulate(4, 5, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[1] >= 19 && XY_[0] == 72)
						FT_.accumulate(5, 6, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 84)
						FT_.accumulate(6, 7, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);

					XY_[0] = XY_[0] + 1;
				}
				// XY_[1] = XY_[1];
			}
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (K_test1)
			{
				std::cout << FindIt << "(after): " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
				std::cout << "boundary_start: " << boundary_start_ << endl;
			}
			return;
		}
	}

	// judge if the particle cross through the third boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 1), &secx, &secy))
	{
		if (!PointandLine(Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(2, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 1), X_new_[0], X_new_[1]))
		{
			FindIt = 3;
			if (K_test1)
			{
				std::cout << FindIt << ": " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
			}
			IfColl_ = 0;
			boundary_start_ = 3;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_)
				{
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 3 in Caltrace()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}
				}
			}
			NumParStat();
			n_Flux_Grid_[XY_[0]][XY_[1]][2] += Weight_ * NumPar_now;
			T_Flux_Grid_[XY_[0]][XY_[1]][2] += Tn_ * Weight_ * NumPar_now;

			if (XY_[1] == 36)
			{
				if (XY_[0] < 12)
					FT_.accumulate(0, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 24)
					FT_.accumulate(1, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 32)
					FT_.accumulate(2, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 48)
					FT_.accumulate(3, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 64)
					FT_.accumulate(4, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 72)
					FT_.accumulate(5, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else if (XY_[0] < 84)
					FT_.accumulate(6, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				else
					FT_.accumulate(7, 9, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
				XY_[0] = -1;
				XY_[1] = -1;
				Zone_ = 6;
				IfColl_ = 0;
				for (int i = 0; i < Grid1.PLasma_Grid_Boundry_num(); i++)
				{
					if (X_[0] >= min(Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i + 1, 0)) && X_[0] <= max(Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i + 1, 0)) && X_[1] >= min(Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 1)) && X_[1] <= max(Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 1)))
					{
						AddPlasmaBoundaryEro(i);
						break;
					}
				}
			}
			else
			{
				if (XY_[1] == 18)
				{
					Zone_ = 3;
				}
				XY_[1] = XY_[1] + 1;
				XY_[0] = XY_[0];
			}
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (K_test1)
			{
				std::cout << FindIt << "(after): " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
				std::cout << "boundary_start: " << boundary_start_ << endl;
			}
			return;
		}
	}

	// judge if the particle cross through the third boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 1), &secx, &secy))
	{
		if (!PointandLine(Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(3, 1), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 0), Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 1), X_new_[0], X_new_[1]))
		{
			FindIt = 4;
			if (K_test1)
			{
				std::cout << FindIt << ": " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
			}
			IfColl_ = 0;
			boundary_start_ = 4;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{ // Rand_flight_ > d_flight_
				if (Rand_flight_ > d_flight_)
				{
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 4 in Caltrace()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}
				}
			}
			NumParStat();
			n_Flux_Grid_[XY_[0]][XY_[1]][3] += Weight_ * NumPar_now;
			T_Flux_Grid_[XY_[0]][XY_[1]][3] += Tn_ * Weight_ * NumPar_now;

			if (XY_[0] == 1)
			{
				IfColl_ = 0;
				IfFlightOut_ = 7;
				boundary_start_ = 0;
				InterscePoint[0][0] = 0;
				InterscePoint[0][1] = X_new_[0];
				InterscePoint[0][2] = X_new_[1];
				InterscePoint[0][3] = XY_[1];
				InterscePoint[0][4] = 1;

				XY_[0] = -1;
				XY_[1] = -1;
				Zone_ = 7;
			}
			else
			{
				if (XY_[1] < 19 && XY_[0] == 25)
				{
					FT_.accumulate(2, 5, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = 72;
				}
				else if (XY_[1] < 19 && XY_[0] == 73)
				{
					FT_.accumulate(6, 1, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = 24;
					Zone_ = 4;
				}
				else
				{
					if (XY_[0] == 13)
						FT_.accumulate(1, 0, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[1] >= 19 && XY_[0] == 25)
						FT_.accumulate(2, 1, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 33)
						FT_.accumulate(3, 2, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 49)
						FT_.accumulate(4, 3, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 65)
						FT_.accumulate(5, 4, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[1] >= 19 && XY_[0] == 73)
						FT_.accumulate(6, 5, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					if (XY_[0] == 85)
						FT_.accumulate(7, 6, Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now);
					XY_[0] = XY_[0] - 1;
				}
				XY_[1] = XY_[1];
			}
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (K_test1)
			{
				std::cout << FindIt << "(after): " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
				std::cout << "boundary_start: " << boundary_start_ << endl;
			}
			return;
		}
	}

	if (!FindIt)
	{
		if (K_test1)
		{
			std::cout << FindIt << "(no)): " << X_[0] << ',' << X_[1] << " " << X_new_[0] << ',' << X_new_[1] << " " << secx << ',' << secy << endl;
			std::cout << "boundary_start: " << boundary_start_ << endl;
		}
		IfColl_ = 1;
		boundary_start_ = 0;
		dt_trace_ = dt_;
		d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
		if (K_flight == 1)
		{ // Rand_flight_ > d_flight_
			if (Rand_flight_ > d_flight_)
			{
				IfColl_ = 1;
				boundary_start_ = 0;
				/*if (this == &D2)
					std::cout << "dt_trace5: " << dt_trace_ << " ";*/
				if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
				{
					std::cout << "error 5 in Caltrace()" << endl;
					pause();
				}
				dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
				/*if (!Grid1.IfinIndex1(XY_[0], XY_[1], X_new_[0], X_new_[1]))
				{
					std::cout << name_ + ": " << X_[0] << "," << X_[1] << " " << X_new_[0] << "," << X_new_[1] << endl;
					std::cout << "the Caltrace() in Particle.cpp have some problem.";
					std::cout << endl;
				}*/
				/*if (this == &D2)
					std::cout << dt_trace_ << endl;*/
				for (int i = 0; i < 3; i++)
				{
					X_new_[i] = X_[i] + V_[i] * dt_trace_;
				}
				NumParStat();
				return;
			}
			else
			{
				// std::cout << "no collision: " << Rand_flight_ << " " << d_flight_ << endl;
				IfColl_ = 0;
			}
		}
		//  std::cout << "IfColl: " << IfColl_ << endl;
		NumParStat();
	}
}

/// @brief calculate the trace of the Particle in triangle mesh
/// @param tag 0 for collision in this grid, 1 for no collision in the grid
void Particle::Caltrace_Tri()
{
	/*if (!Grid4.Ifingrid(Tri_Index_, X_[0], X_[1]))
	{
		std::cout << "ni zai na?" << endl;
		pause();
	}*/
	bool Track_test = 0;
	if (Track_test)
	{
		std::cout << "tiaoshi1: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
		std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
		std::cout << Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r << ", " << Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z << ", " << Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r << ", " << Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z << endl;
	}
	int FindIt = 0;
	boundary_start_ = 0;
	IfColl_ = 1;
	double secx, secy, distance;
	// judge if the particle cross through the FIRST boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z, &secx, &secy))
	{
		if (!PointandLine(Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z, X_new_[0], X_new_[1]))
		{
			FindIt = 1;
			boundary_start_ = 1;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_) // collision occur
				{
					if (Track_test)
					{
						std::cout << "tiaoshi2: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 1 in Caltrace_Tri()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					if (Track_test)
					{
						std::cout << "tiaoshi6: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}

					NumParStat();
					if (Grid4.lines_info(Tri_Index_, 0) == 0) // next
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[0];
					}
					else if (Grid4.lines_info(Tri_Index_, 0) == -1 || Grid4.lines_info(Tri_Index_, 0) == -2) // Wall
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = 0;
						InterscePoint[0][4] = 11;
						Zone_ = 7;
					}
					else if (Grid4.lines_info(Tri_Index_, 0) == 1) // Core Boundary
					{
						IfFlightOut_ = 1;
						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = XY_[0] - 24;
						InterscePoint[0][4] = 18;

						boundary_start_ = 0;
						Tri_Index_ = -1;
						Zone_ = 1;
						return;
					}
					else if (Grid4.lines_info(Tri_Index_, 0) == 2) // Plasma Boundary
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[0];
						if (Zone_ < 6)
						{
							Zone_ = 6;
						}
						else if (Zone_ >= 6)
						{
							Zone_ = 5;
						}
					}
					else if (Grid4.lines_info(Tri_Index_, 0) == 3 || Grid4.lines_info(Tri_Index_, 0) == 4) // Inner Target or Outer Target
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;
						IfColl_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = 0;
						/*if (Grid4.lines_info(Tri_Index_, 0) == 3)
							InterscePoint[0][3] = XY_[1];
						else if (Grid4.lines_info(Tri_Index_, 0) == 4)
							InterscePoint[0][3] = XY_[1] + N_radial;*/
						InterscePoint[0][4] = 1;
						Zone_ = 7;
						return;
					}

					if (Track_test)
					{
						std::cout << "New index: " << Tri_Index_ << endl;
					}
					if (Zone_ < 6)
					{
						XY_[0] = Tri_B2_[Tri_Index_][0];
						XY_[1] = Tri_B2_[Tri_Index_][1];
					}
					else
					{
						XY_[0] = -1;
						XY_[1] = -1;
					}
					return;
				}
			}
		}
		else
		{
			if (Track_test)
				std::cout << "1weicengchuqu?" << endl;
		}
	}

	// judge if the particle cross through the SECOND boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].z, &secx, &secy))
	{
		if (!PointandLine(Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[1]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].z, X_new_[0], X_new_[1]))
		{
			FindIt = 1;
			boundary_start_ = 2;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_)
				{
					if (Track_test)
					{
						std::cout << "tiaoshi3: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 2 in Caltrace_Tri()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					if (Track_test)
					{
						std::cout << "tiaoshi7: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}

					NumParStat();
					if (Grid4.lines_info(Tri_Index_, 1) == 0) // next
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[3];
					}
					else if (Grid4.lines_info(Tri_Index_, 1) == -1 || Grid4.lines_info(Tri_Index_, 1) == -2) // Wall
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = 1;
						InterscePoint[0][4] = 11;
						Zone_ = 7;
					}
					else if (Grid4.lines_info(Tri_Index_, 1) == 1) // Core Boundary
					{
						IfFlightOut_ = 1;
						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = XY_[0] - 24;
						InterscePoint[0][4] = 18;

						boundary_start_ = 0;
						XY_[0] = -1;
						XY_[1] = -1;
						Tri_Index_ = -1;
						Zone_ = 1;
						return;
					}
					else if (Grid4.lines_info(Tri_Index_, 1) == 2) // Plasma Boundary
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[3];
						if (Zone_ < 6)
						{
							Zone_ = 6;
							XY_[0] = -1;
							XY_[1] = -1;
						}
						else if (Zone_ >= 6)
						{
							Zone_ = 5;
						}
					}
					else if (Grid4.lines_info(Tri_Index_, 1) == 3 || Grid4.lines_info(Tri_Index_, 1) == 4) // Inner Target or Outer Target
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;
						IfColl_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						/*if (Grid4.lines_info(Tri_Index_, 1) == 3)
							InterscePoint[0][3] = XY_[1];
						else if (Grid4.lines_info(Tri_Index_, 1) == 4)
							InterscePoint[0][3] = XY_[1] + N_radial;*/
						InterscePoint[0][3] = 1;
						InterscePoint[0][4] = 1;
						Zone_ = 7;
						return;
					}

					if (Track_test)
					{
						std::cout << "New index: " << Tri_Index_ << endl;
					}
					if (Zone_ < 6)
					{
						XY_[0] = Tri_B2_[Tri_Index_][0];
						XY_[1] = Tri_B2_[Tri_Index_][1];
					}
					else
					{
						XY_[0] = -1;
						XY_[1] = -1;
					}
					return;
				}
			}
		}
		else
		{
			if (Track_test)
				std::cout << "2weicengchuqu?" << endl;
		}
	}

	// judge if the particle cross through the THIRD boundry of grid
	if (!FindIt && get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1], Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z, &secx, &secy))
	{
		if (!PointandLine(Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[2]].z, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].r, Grid4.nodes_[Grid4.tris_[Tri_Index_].v[0]].z, X_new_[0], X_new_[1]))
		{
			FindIt = 1;
			boundary_start_ = 3;
			X_new_[0] = secx;
			X_new_[1] = secy;
			dt_trace_ = (X_new_[0] - X_[0]) / V_[0];
			d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
			if (K_flight == 1)
			{
				if (Rand_flight_ > d_flight_)
				{
					if (Track_test)
					{
						std::cout << "tiaoshi4: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 1;
					boundary_start_ = 0;
					if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
					{
						std::cout << "error 3 in Caltrace_Tri()" << endl;
						pause();
					}
					dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
					/*if (this == &D2)
						std::cout << "dt_trace4: " << dt_trace_ << endl;*/
					for (int i = 0; i < 3; i++)
					{
						X_new_[i] = X_[i] + V_[i] * dt_trace_;
					}
					NumParStat();
					return;
				}
				else
				{
					if (Track_test)
					{
						std::cout << "tiaoshi8: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
						std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
					}
					IfColl_ = 0;
					if (dt_ < 0)
					{
						std::cout << "dt < 0" << endl;
						exit(0);
					}

					NumParStat();
					if (Grid4.lines_info(Tri_Index_, 2) == 0) // next
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[6];
					}
					else if (Grid4.lines_info(Tri_Index_, 2) == -1 || Grid4.lines_info(Tri_Index_, 2) == -2) // Wall
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = 2;
						InterscePoint[0][4] = 11;
						Zone_ = 7;
					}
					else if (Grid4.lines_info(Tri_Index_, 2) == 1) // Core Boundary
					{
						IfFlightOut_ = 1;
						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						InterscePoint[0][3] = XY_[0] - 24;
						InterscePoint[0][4] = 18;

						boundary_start_ = 0;
						XY_[0] = -1;
						XY_[1] = -1;
						Tri_Index_ = -1;
						Zone_ = 1;
						return;
					}
					else if (Grid4.lines_info(Tri_Index_, 2) == 2) // Plasma Boundary
					{
						Tri_Index_ = Grid4.tris_[Tri_Index_].neigh[6];
						if (Zone_ < 6)
						{
							Zone_ = 6;
							XY_[0] = -1;
							XY_[1] = -1;
						}
						else if (Zone_ >= 6)
						{
							Zone_ = 5;
						}
					}
					else if (Grid4.lines_info(Tri_Index_, 2) == 3 || Grid4.lines_info(Tri_Index_, 2) == 4) // Inner Target or Outer Target
					{
						IfFlightOut_ = 7;
						boundary_start_ = 0;
						IfColl_ = 0;

						InterscePoint[0][0] = 0;
						InterscePoint[0][1] = secx;
						InterscePoint[0][2] = secy;
						/*if (Grid4.lines_info(Tri_Index_, 2) == 3)
							InterscePoint[0][3] = XY_[1];
						else if (Grid4.lines_info(Tri_Index_, 2) == 4)
							InterscePoint[0][3] = XY_[1] + N_radial;*/
						InterscePoint[0][3] = 2;
						InterscePoint[0][4] = 1;
						Zone_ = 7;
						return;
					}

					if (Track_test)
					{
						std::cout << "New index: " << Tri_Index_ << endl;
					}
					if (Zone_ < 6)
					{
						XY_[0] = Tri_B2_[Tri_Index_][0];
						XY_[1] = Tri_B2_[Tri_Index_][1];
					}
					else
					{
						XY_[0] = -1;
						XY_[1] = -1;
					}
					return;
				}
			}
		}
		else
		{
			if (Track_test)
				std::cout << "3weicengchuqu?" << endl;
		}
	}

	if (!FindIt)
	{
		if (!Grid4.Ifingrid(Tri_Index_, X_new_[0], X_new_[1]))
		{
			std::cout << "nimadewsm?" << endl;
			pause();
		}
		IfColl_ = 1;
		boundary_start_ = 0;
		dt_trace_ = dt_;
		d_flight_ = d_flight_ * exp(-dt_trace_ / dt_);
		if (K_flight == 1)
		{ // Rand_flight_ > d_flight_
			if (Rand_flight_ > d_flight_)
			{
				if (Track_test)
				{
					std::cout << "tiaoshi5: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
					std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
				}
				IfColl_ = 1;
				boundary_start_ = 0;
				/*if (this == &D2)
					std::cout << "dt_trace5: " << dt_trace_ << " ";*/
				if (log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_))) < -1)
				{
					std::cout << "error 4 in Caltrace_Tri()" << endl;
					pause();
				}
				dt_trace_ = -dt_trace_ * log(Rand_flight_ / (d_flight_ / exp(-dt_trace_ / dt_)));
				for (int i = 0; i < 3; i++)
				{
					X_new_[i] = X_[i] + V_[i] * dt_trace_;
				}
				NumParStat();
				return;
			}
			else
			{
				if (Track_test)
				{
					std::cout << "tiaoshi9: " << name_ << "\t" << Tri_Index_ << "\t" << X_[0] << ", " << X_[1] << "\t" << X_[2] << "\t";
					std::cout << X_new_[0] << ", " << X_new_[1] << "\t" << X_new_[2] << endl;
				}
				// std::cout << "no collision: " << Rand_flight_ << " " << d_flight_ << endl;
				IfColl_ = 0;
			}
		}
	}
}

void Particle::NumParStat()
{
	X_[0] = X_new_[0];
	X_[1] = X_new_[1];
	X_[2] = X_new_[2];
	// std::cout << "dt_trace: " << dt_trace_ << endl;
	if (Zone_ < 6 && Zone_ > 1)
	{
		if (Charge_ == 0)
		{
			if (K_flight == 1 || K_flight == 3)
			{
				if (K_T_array == 1)
				{
					FT_.accumulateGrid(XY_[0], XY_[1], 1.5 * Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now * dt_trace_);
				}
				Sum_n_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now * dt_trace_;
				Sum_E_[XY_[0]][XY_[1]][Charge_] += 1.5 * Weight_ * NumPar_now * Tn_ * dt_trace_;
				for (int i = 0; i < 3; i++)
					Sum_V_[XY_[0]][XY_[1]][i][Charge_] += V_[i] * Weight_ * NumPar_now * dt_trace_;
			}
			else if (K_flight == 2)
			{
				Sum_n_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now * dt;
				T_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now * Tn_ * dt;
				for (int i = 0; i < 3; i++)
					Sum_V_[XY_[0]][XY_[1]][i][Charge_] += V_[i] * Weight_ * NumPar_now * dt;
			}
		}

		if (Charge_ > 0)
		{
			Sum_V_[XY_[0]][XY_[1]][0][Charge_] += V_Charge_[3] * Weight_ * NumPar_now * dt;
			Sum_V_[XY_[0]][XY_[1]][1][Charge_] += V_Charge_[7] * Weight_ * NumPar_now * dt;
		}
	}
	if (MeshMode == 3)
	{
		// std::cout << "dt_trace: " << dt_trace_ << endl
		if (Charge_ == 0)
		{
			if (K_flight == 1 || K_flight == 3)
			{
				/*if (K_T_array == 1)
				{
					FT_.accumulateGrid(XY_[0], XY_[1], 1.5 * Tn_, V_[0], V_[1], V_[2], Weight_ * NumPar_now * dt_trace_);
				}*/
				Tri_Sum_n_[Tri_Index_][Charge_] += Weight_ * NumPar_now * dt_trace_;
				Tri_Sum_E_[Tri_Index_][Charge_] += 1.5 * Weight_ * NumPar_now * Tn_ * dt_trace_;
				for (int i = 0; i < 3; i++)
					Tri_Sum_V_[Tri_Index_][i][Charge_] += V_[i] * Weight_ * NumPar_now * dt_trace_;
			}
			else if (K_flight == 2)
			{
				Tri_Sum_n_[Tri_Index_][Charge_] += Weight_ * NumPar_now * dt;
				Tri_Sum_E_[Tri_Index_][Charge_] += Weight_ * NumPar_now * Tn_ * dt;
				for (int i = 0; i < 3; i++)
					Tri_Sum_V_[Tri_Index_][i][Charge_] += V_[i] * Weight_ * NumPar_now * dt;
			}
		}

		if (Charge_ > 0)
		{
			Tri_Sum_V_[Tri_Index_][0][Charge_] += V_Charge_[3] * Weight_ * NumPar_now * dt;
			Tri_Sum_V_[Tri_Index_][1][Charge_] += V_Charge_[7] * Weight_ * NumPar_now * dt;
		}
	}
}

void Particle::CalProb()
{
	if (this == &H || this == &D || this == &T)
	{
		if (MeshMode == 3)
		{
			Ion_[0].SetTrics_cs(Tri_B2_);
			CX_[0].SetTrics_cs(Tri_B2_);
			if (K_CX_DT)
			{
				CX_DT_[0].SetTrics_cs(Tri_B2_);
			}
			for (int i = 0; i < num_trimesh_; i++)
			{
				Tri_Totalcs_[i][0] = Ion_[0].cs(i) + CX_[0].cs(i) + CX_DT_[0].cs(i) + R_with_H_[0].cs(i) + R_with_H2_[0].cs(i);
				Tri_CollProb_[i][0] = 1 - exp(-dt * Tri_Totalcs_[i][0]);
				Tri_lambda_[i][0] = 1 / Tri_Totalcs_[i][0];
			}
		}
		else
		{
			for (int i = 1; i < N_poloidal - 1; i++)
			{
				for (int j = 1; j < N_radial - 1; j++)
				{
					if (K_Prob == 1)
					{
						Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) + CX_DT_[0].cs(i, j);
						// std::cout << name_ + " cross: " << i << ", " << j << endl;
						// std::cout << Ion_[0].cs(i, j) << ", " << CX_[0].cs(i, j) << ", " << CX_DT_[0].cs(i, j) << endl;
						CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
						lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
						if (K_flight == 2)
						{
							if (i == 1 && j == 1)
							{
								lambda_min_[0] = lambda_[i][j][0];
							}
							else
							{
								lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
							}
						}
					}
					else if (K_Prob == 2)
					{
						CollProb_[i][j][0] = 1 - exp(-dt * Ion_[0].cs(i, j)) + 1 - exp(-dt * CX_[0].cs(i, j));
					}
					// CollProb_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j);
					//  std::cout << CollProb_[i][j][0] << '\t';
				}
				// std::cout << endl;
			}
			for (int i = 0; i < N_poloidal; i++)
			{
				for (int j = 0; j < N_radial; j++)
				{
					if (i == 0)
					{
						if (j == 0)
						{
							lambda_[i][j][0] = lambda_[i + 1][j + 1][0];
						}
						else if (j == N_radial - 1)
						{
							lambda_[i][j][0] = lambda_[i + 1][j - 1][0];
						}
						else
							lambda_[i][j][0] = lambda_[i + 1][j][0];
					}
					else
					{
						if (j == 0)
						{
							lambda_[i][j][0] = lambda_[i][j + 1][0];
						}
						if (j == N_radial - 1)
						{
							lambda_[i][j][0] = lambda_[i][j - 1][0];
						}
						if (i == N_poloidal - 1)
						{
							lambda_[i][j][0] = lambda_[i - 1][j][0];
						}
					}
				}
			}
		}
		CoreTotalcs_[0] = Ion_[0].Cor_cs() + CX_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
	}
	else if (this == &H2 || this == &D2 || this == &T2)
	{
		if (MeshMode == 3)
		{
			if (K_Prob == 1)
			{
				Ion_[0].SetTrics_cs(Tri_B2_);
				CX_[0].SetTrics_cs(Tri_B2_);
				Diss1_[0].SetTrics_cs(Tri_B2_);
				Diss2_[0].SetTrics_cs(Tri_B2_);
				Ela_[0].SetTrics_cs(Tri_B2_);

				DS_[1][0].SetTrics_cs(Tri_B2_);
				DS_[1][1].SetTrics_cs(Tri_B2_);
				DS_[1][2].SetTrics_cs(Tri_B2_);

				for (int i = 0; i < num_trimesh_; i++)
				{
					Tri_Totalcs_[i][0] = Ion_[0].cs(i) + CX_[0].cs(i) + Diss1_[0].cs(i) + Diss2_[0].cs(i) + Ela_[0].cs(i) + R_with_H_[0].cs(i) + R_with_H2_[0].cs(i);
					Tri_CollProb_[i][0] = 1 - exp(-dt * Tri_Totalcs_[i][0]);
					Tri_lambda_[i][0] = 1 / Tri_Totalcs_[i][0];

					Tri_Totalcs_[i][1] = DS_[1][0].cs(i) + DS_[1][1].cs(i) + DS_[1][2].cs(i);
					Tri_CollProb_[i][1] = 1 - exp(-dt * Tri_Totalcs_[i][1]);
					Tri_lambda_[i][1] = 1 / Tri_Totalcs_[i][1];
				}
			}
		}
		else if (MeshMode == 1)
		{
			for (int i = 1; i < N_poloidal - 1; i++)
			{
				for (int j = 1; j < N_radial - 1; j++)
				{
					if (K_Prob == 1)
					{
						Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) + Diss1_[0].cs(i, j) + Diss2_[0].cs(i, j) + Ela_[0].cs(i, j);
						CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
						lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
						Totalcs_[i][j][1] = DS_[1][0].cs(i, j) + DS_[1][1].cs(i, j) + DS_[1][2].cs(i, j);
						CollProb_[i][j][1] = 1 - exp(-dt * Totalcs_[i][j][1]);
						lambda_[i][j][1] = 1 / Totalcs_[i][j][1];
						if (i == 1 && j == 1)
						{
							lambda_min_[0] = lambda_[i][j][0];
						}
						else
						{
							lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
						}
					}
					else if (K_Prob == 2)
					{
						CollProb_[i][j][0] =
							1 - exp(-dt * Ion_[0].cs(i, j)) + 1 -
							exp(-dt * CX_[0].cs(i, j)) + 1 - exp(-dt * Diss1_[0].cs(i, j)) + 1 -
							exp(-dt * Diss2_[0].cs(i, j)) + 1 - exp(-dt * Ela_[0].cs(i, j));
						CollProb_[i][j][1] = 1 - exp(-dt * DS_[1][0].cs(i, j)) + 1 -
											 exp(-dt * DS_[1][1].cs(i, j)) + 1 -
											 exp(-dt * DS_[1][2].cs(i, j));
					}
				}
			}
			for (int i = 0; i < N_poloidal; i++)
			{
				for (int j = 0; j < N_radial; j++)
				{
					if (i == 0)
					{
						lambda_[i][j][0] = lambda_[i + 1][j][0];
						lambda_[i][j][1] = lambda_[i + 1][j][1];
						if (j == 0)
						{
							lambda_[i][j][0] = lambda_[i + 1][j + 1][0];
							lambda_[i][j][1] = lambda_[i + 1][j + 1][1];
						}
						if (j == N_radial - 1)
						{
							lambda_[i][j][0] = lambda_[i + 1][j - 1][0];
							lambda_[i][j][1] = lambda_[i + 1][j - 1][1];
						}
					}
					else
					{
						if (j == 0)
						{
							lambda_[i][j][0] = lambda_[i][j + 1][0];
							lambda_[i][j][1] = lambda_[i][j + 1][1];
						}
						if (j == N_radial - 1)
						{
							lambda_[i][j][0] = lambda_[i][j - 1][0];
							lambda_[i][j][1] = lambda_[i][j - 1][1];
						}
						if (i == N_poloidal - 1)
						{
							lambda_[i][j][0] = lambda_[i - 1][j][0];
							lambda_[i][j][1] = lambda_[i - 1][j][1];
						}
					}
				}
			}
		}
		CoreTotalcs_[0] = Ion_[0].Cor_cs() + CX_[0].Cor_cs() + Diss1_[0].Cor_cs() +
						  Diss2_[0].Cor_cs() + Ela_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
		CoreTotalcs_[1] = DS_[1][0].Cor_cs() + DS_[1][1].Cor_cs() + DS_[1][2].Cor_cs();
		CoreCollProb_[1] = 1 - exp(-dt * CoreTotalcs_[1]);
	}
	else if (this == &CD4)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) +
									Diss1_[0].cs(i, j) + Diss2_[0].cs(i, j);
				CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
				lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
				Totalcs_[i][j][1] = DS_[1][0].cs(i, j) + DS_[1][1].cs(i, j) +
									DS_[1][2].cs(i, j) + DS_[1][3].cs(i, j);
				CollProb_[i][j][1] = 1 - exp(-dt * Totalcs_[i][j][1]);
				lambda_[i][j][1] = 1 / Totalcs_[i][j][1];
				if (i == 1 && j == 1)
				{
					lambda_min_[0] = lambda_[i][j][0];
				}
				else
				{
					lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
				}
			}
		CoreTotalcs_[0] = Ion_[0].Cor_cs() + CX_[0].Cor_cs() + Diss1_[0].Cor_cs() + Diss2_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
		CoreTotalcs_[1] = DS_[1][0].Cor_cs() + DS_[1][1].Cor_cs() + DS_[1][2].Cor_cs() + DS_[1][3].Cor_cs();
		CoreCollProb_[1] = 1 - exp(-dt * CoreTotalcs_[1]);
	}
	else if (this == &CD3)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) +
									Diss1_[0].cs(i, j) + Diss2_[0].cs(i, j);
				CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
				lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
				Totalcs_[i][j][1] = DS_[1][0].cs(i, j) + DS_[1][1].cs(i, j) + DS_[1][2].cs(i, j);
				CollProb_[i][j][1] = 1 - exp(-dt * Totalcs_[i][j][1]);
				lambda_[i][j][1] = 1 / Totalcs_[i][j][1];
				if (i == 1 && j == 1)
				{
					lambda_min_[0] = lambda_[i][j][0];
				}
				else
				{
					lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
				}
			}
		CoreTotalcs_[0] =
			Ion_[0].Cor_cs() + CX_[0].Cor_cs() + Diss1_[0].Cor_cs() + Diss2_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
		CoreTotalcs_[1] = DS_[1][0].Cor_cs() + DS_[1][1].Cor_cs() + DS_[1][2].Cor_cs();
		CoreCollProb_[1] = 1 - exp(-dt * CoreTotalcs_[1]);
	}
	else if (this == &CD2)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) +
									Diss1_[0].cs(i, j) + Diss2_[0].cs(i, j);
				CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
				lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
				Totalcs_[i][j][1] = DS_[1][0].cs(i, j) + DS_[1][1].cs(i, j) + DS_[1][2].cs(i, j);
				CollProb_[i][j][1] = 1 - exp(-dt * Totalcs_[i][j][1]);
				lambda_[i][j][1] = 1 / Totalcs_[i][j][1];
				if (i == 1 && j == 1)
				{
					lambda_min_[0] = lambda_[i][j][0];
				}
				else
				{
					lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
				}
			}
		CoreTotalcs_[0] = Ion_[0].Cor_cs() + CX_[0].Cor_cs() + Diss1_[0].Cor_cs() + Diss2_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
		CoreTotalcs_[1] = DS_[1][0].Cor_cs() + DS_[1][1].Cor_cs() + DS_[1][2].Cor_cs();
		CoreCollProb_[1] = 1 - exp(-dt * CoreTotalcs_[1]);
	}
	else if (this == &CD1)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				Totalcs_[i][j][0] = Ion_[0].cs(i, j) + CX_[0].cs(i, j) +
									Diss1_[0].cs(i, j) + Diss2_[0].cs(i, j) + Diss3_[0].cs(i, j);
				CollProb_[i][j][0] = 1 - exp(-dt * Totalcs_[i][j][0]);
				lambda_[i][j][0] = 1 / Totalcs_[i][j][0];
				Totalcs_[i][j][1] = DS_[1][0].cs(i, j) + DS_[1][1].cs(i, j) + DS_[1][2].cs(i, j);
				CollProb_[i][j][1] = 1 - exp(-dt * Totalcs_[i][j][1]);
				lambda_[i][j][1] = 1 / Totalcs_[i][j][1];
				if (i == 1 && j == 1)
				{
					lambda_min_[0] = lambda_[i][j][0];
				}
				else
				{
					lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
				}
			}
		CoreTotalcs_[0] = Ion_[0].Cor_cs() + CX_[0].Cor_cs() + Diss1_[0].Cor_cs() + Diss2_[0].Cor_cs() + Diss3_[0].Cor_cs();
		CoreCollProb_[0] = 1 - exp(-dt * CoreTotalcs_[0]);
		CoreTotalcs_[1] = DS_[1][0].Cor_cs() + DS_[1][1].Cor_cs() + DS_[1][2].Cor_cs();
		CoreCollProb_[1] = 1 - exp(-dt * CoreTotalcs_[1]);
	}
	else if (this == &C)
	{
		for (int k = 0; k < MaxCharge_ + 1; k++)
		{
			for (int i = 0; i < N_poloidal; i++)
				for (int j = 0; j < N_radial; j++)
				{
					Totalcs_[i][j][k] = Ion_[k].cs(i, j) + CX_[k].cs(i, j) + Rec_[k].cs(i, j);
					CollProb_[i][j][k] = 1 - exp(-dt * Totalcs_[i][j][k]);
					lambda_[i][j][k] = 1 / Totalcs_[i][j][k];
					if (i == 1 && j == 1)
					{
						lambda_min_[0] = lambda_[i][j][0];
					}
					else
					{
						lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
					}
				}
			CoreTotalcs_[k] = Ion_[k].Cor_cs() + CX_[k].Cor_cs() + Rec_[k].Cor_cs();
			CoreCollProb_[k] = 1 - exp(-dt * CoreTotalcs_[k]);
		}
	}
	else if (this == &Ar)
	{
		for (int k = 0; k < MaxCharge_ + 1; k++)
		{
			for (int i = 0; i < N_poloidal; i++)
				for (int j = 0; j < N_radial; j++)
				{
					Totalcs_[i][j][k] = Ion_[k].cs(i, j) + CX_[k].cs(i, j) + Rec_[k].cs(i, j);
					CollProb_[i][j][k] = 1 - exp(-dt * Totalcs_[i][j][k]);
					lambda_[i][j][k] = 1 / Totalcs_[i][j][k];
					if (i == 1 && j == 1)
					{
						lambda_min_[0] = lambda_[i][j][0];
					}
					else
					{
						lambda_min_[0] = min(lambda_min_[0], lambda_[i][j][0]);
					}
				}
			CoreTotalcs_[k] = Ion_[k].Cor_cs() + CX_[k].Cor_cs() + Rec_[k].Cor_cs();
			CoreCollProb_[k] = 1 - exp(-dt * CoreTotalcs_[k]);
		}
	}
}

double Particle::CollProb()
{
	if (Zone_ == 1)
		return CoreCollProb_[Charge_];
	else if (Zone_ < 7)
		return CollProb_[XY_[0]][XY_[1]][Charge_];
	else
		throw std::logic_error("Invalid Zone_ in CollProb()");
}

void Particle::Coll()
{
	/*if (Zone_ > 5)
	{
		std::cout << "Particle shouldn't collde here!" << endl;
		exit(6);
	}*/
	/*if (X_[0] != X_new_[0])
	{
		std::cout << "collision error;" << endl;
		exit(0);
		pause();
	}*/

	double V_temp1; // parallel velocity of D ion particle before collision
	double V_temp2; // parallel velocity of neutral particle before collision
	double V_temp3; // parallel velocity of T ion particle before collision
	if (Zone_ < 6)
	{
		if (K_D)
		{
			V_temp1 = B[XY_[0]][XY_[1]][0] * V_D_1_now[0] + B[XY_[0]][XY_[1]][1] * V_D_1_now[1] + B[XY_[0]][XY_[1]][2] * V_D_1_now[2];
		}
		if (K_T)
		{
			V_temp3 = B[XY_[0]][XY_[1]][0] * V_T_1_now[0] + B[XY_[0]][XY_[1]][1] * V_T_1_now[1] + B[XY_[0]][XY_[1]][2] * V_T_1_now[2];
		}
	}
	if (this == &H || this == &D || this == &T)
	{
		if (Charge_ == MaxCharge_)
		{
			Ion_factor = 0;
			CX_factor = 0;
			R_with_H_factor = 0;
			R_with_H2_factor = 0;
		}
		else
		{
			if (MeshMode == 1)
			{
				Ion_factor = Ion_[Charge_].cs(XY_, Zone_);
				if (K_Tn == 1)
				{
					CX_factor = CX_[Charge_].cs_now();
					if (K_CX_DT)
						CX_DT_factor = CX_DT_[Charge_].cs_now();
					else
						CX_DT_factor = 0.;
				}
				else if (K_Tn == 2)
				{
					CX_factor = CX_[Charge_].cs(XY_, Zone_);
					if (K_CX_DT)
						CX_DT_factor = CX_DT_[Charge_].cs(XY_, Zone_);
					else
						CX_DT_factor = 0.;
				}
			}
			else if (MeshMode == 3)
			{
				if (Zone_ < 6)
				{
					Ion_factor = Ion_[Charge_].cs(Tri_Index_);
					if (K_Tn == 1)
					{
						CX_factor = CX_[Charge_].cs_now();
						if (K_CX_DT)
							CX_DT_factor = CX_DT_[Charge_].cs_now();
						else
							CX_DT_factor = 0.;
					}
					else if (K_Tn == 2)
					{
						CX_factor = CX_[Charge_].cs(Tri_Index_);
						if (K_CX_DT)
							CX_DT_factor = CX_DT_[Charge_].cs(Tri_Index_);
						else
							CX_DT_factor = 0.;
					}
				}
				else
				{
					Ion_factor = 0;
					CX_factor = 0;
					CX_DT_factor = 0;
				}
				R_with_H_factor = R_with_H_[Charge_].cs(Tri_Index_);
				R_with_H2_factor = R_with_H2_[Charge_].cs(Tri_Index_);
			}
			// std::cout << Ion_factor << '\t' << CX_factor << '\t' << CX_DT_factor << endl;
		}
		double rand_temp, Var_temp;
		if (K_flight == 1)
		{
			rand_temp = 0;
			Var_temp = 1;
		}
		else if (K_flight == 2)
		{
			rand_temp = Random();
			Var_temp = CollProb();
		}
		else if (K_flight == 3)
		{
			rand_temp = Random();
			if (K_Tn == 1)
			{
				Var_temp = lambda_min_[0] / lambda_now_;
			}
			else if (K_Tn == 2)
			{
				Var_temp = lambda_min_[0] / lambda_[XY_[0]][XY_[1]][0];
			}
			// std::cout << "CollProb: " << rand_temp << " " << Var_temp << endl;
		}
		double sigma_all;
		sigma_all = Ion_factor + CX_factor + CX_DT_factor + R_with_H_factor + R_with_H2_factor;
		// std::cout << rand_temp << '\t' << CollProb() << endl;
		if (rand_temp < Var_temp)
		{
			// std::cout << rand_temp << '\t' << CollProb() << endl;
			// std::cout << this->name() << '\t' << Tn_ << '\t';
			double rand_one = Random();
			if (rand_one < Ion_factor / sigma_all)
			{

				// std::cout << name_ << ": " << XY_[0] << ", " << XY_[1] << ", " << XY_[2] << endl;
				// std::cout << "X: " << ": " << X_[0] << ", " << X_[1] << ", " << X_[2] << endl;
				// std::cout << "1" << '\t';
				if (Zone_ > 5)
				{
					std::cout << "Collsion error 1 in Ion" << endl;
				}
				V_temp2 = B[XY_[0]][XY_[1]][0] * V_[0] + B[XY_[0]][XY_[1]][1] * V_[1] + B[XY_[0]][XY_[1]][2] * V_[2];

				if (MeshMode == 1)
				{
					Ion_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
					if (K_mu == 2 || K_mu == 3)
					{
						Ion_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
						Ion_[Charge_].E_Add(XY_, 1.5 * Tn_ * qe * Weight_ * NumPar_now);
					}
					else if (K_mu == 1)
					{
						Ion_[Charge_].Mu_Add(XY_, mass_ * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
						Ion_[Charge_].E_Add(XY_, (1.5 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
					}
				}
				else if (MeshMode == 3)
				{
					Ion_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
					if (K_mu == 2 || K_mu == 3)
					{
						Ion_[Charge_].Mu_Add(Tri_Index_, mass_ * V_temp2 * Weight_ * NumPar_now);
						Ion_[Charge_].E_Add(Tri_Index_, 1.5 * Tn_ * qe * Weight_ * NumPar_now);
					}
					else if (K_mu == 1)
					{
						Ion_[Charge_].Mu_Add(Tri_Index_, mass_ * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
						Ion_[Charge_].E_Add(Tri_Index_, (1.5 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
					}
				}

				setfate(1, 5, Index_);
				Weight_ = 0;
				// std::cout << Weight_ << endl;
				return;
			}
			else if (rand_one < (Ion_factor + CX_factor) / sigma_all)
			{
				if (Zone_ > 5)
				{
					std::cout << "Collsion error 1 in Ion" << endl;
				}
				if (MeshMode == 1)
				{
					CX_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
					for (int i = 0; i < 3; i++)
					{
						Sum_V_CX_Ion_Be_[XY_[0]][XY_[1]][i][0] += V_D_1_now[i] * Weight_ * NumPar_now * dt_trace_;
						Sum_V_CX_Neu_Be_[XY_[0]][XY_[1]][i][0] += V_[i] * Weight_ * NumPar_now * dt_trace_;
					}
					if (this == &D)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_[Charge_].Mu_Add(XY_, (-Dmass * ua_D_1[XY_[0]][XY_[1]]) * Weight_ * NumPar_now);
							CX_[Charge_].E_Add(XY_, -(1.5 * Ti[XY_[0]][XY_[1]] * qe + 0.5 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							CX_[Charge_].Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							CX_[Charge_].E_Add(XY_, (-1.5 * Ti[XY_[0]][XY_[1]] * qe - 0.5 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
					}
					if (this == &T)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_[Charge_].Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_[Charge_].E_Add(XY_, -(1.5 * Ti[XY_[0]][XY_[1]] * qe + 0.5 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							CX_[Charge_].Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							CX_[Charge_].E_Add(XY_, (-1.5 * Ti[XY_[0]][XY_[1]] * qe + -0.5 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
					}
					// std::cout << Index_ << endl;
					setfate(0, 1, Index_); // CX with D
					// std::cout << name_ + " before CX V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					Init(2, 2);
					if (K_mu == 2)
					{
						V_temp2 = B[XY_[0]][XY_[1]][0] * V_1[0] + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2];
						CX_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
						CX_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
					}
					if (K_mu == 3)
					{
						if (this == &D)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * (V_1[0] - V_temp1) + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2] - V_temp1;
						}
						if (this == &T)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * (V_1[0] - V_temp1) + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2] - V_temp3;
						}
						CX_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
						if (V_temp2 < 0)
						{
							CX_[Charge_].E_Add(XY_, -(1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (V_temp2 >= 0)
						{
							CX_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
					}
				}
				else if (MeshMode == 3)
				{
					CX_[Charge_].n_Add(Tri_Index_, Weight_ * NumPar_now);
					if (K_mu == 3)
					{
						if (this == &D)
						{
							CX_[Charge_].Mu_Add(Tri_Index_, -Dmass * Weight_ * NumPar_now * (CX_[Charge_].V_relative(0) * B[XY_[0]][XY_[1]][0] + CX_[Charge_].V_relative(1) * B[XY_[0]][XY_[1]][1] + CX_[Charge_].V_relative(2) * B[XY_[0]][XY_[1]][2]));
							CX_[Charge_].E_Add(Tri_Index_, -0.5 * Dmass * CX_[0].V_2_relative() * Weight_ * NumPar_now);
						}
						if (this == &T)
						{
							CX_[Charge_].Mu_Add(Tri_Index_, -Tmass * Weight_ * NumPar_now * (CX_[Charge_].V_relative(0) * B[XY_[0]][XY_[1]][0] + CX_[Charge_].V_relative(1) * B[XY_[0]][XY_[1]][1] + CX_[Charge_].V_relative(2) * B[XY_[0]][XY_[1]][2]));
							CX_[Charge_].E_Add(Tri_Index_, -0.5 * Tmass * CX_[0].V_2_relative() * Weight_ * NumPar_now);
						}
					}
					for (int i = 0; i < 3; i++)
						SetV(i, V_D_1_now[i]);
					setfate(0, 1, Index_); // CX with D
				}
				CalTn();
				// std::cout << name_ + " after CX V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
				//  std::cout << "CX" << endl;
				//   std::cout << P->Tn() << '\t' << P->name() << endl;
				for (int i = 0; i < 3; i++)
				{
					Sum_V_CX_Ion_Af_[XY_[0]][XY_[1]][i][0] += V_1[i] * Weight_ * NumPar_now * dt_trace_;
					Sum_V_CX_Neu_Af_[XY_[0]][XY_[1]][i][0] += V_[i] * Weight_ * NumPar_now * dt_trace_;
				}
				return;
			}
			else if (rand_one < (Ion_factor + CX_factor + CX_DT_factor) / sigma_all)
			{
				// std::cout << "???" << endl;
				if (MeshMode == 1)
				{
					if (this == &D)
					{
						CX_DT_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_DT_[Charge_].Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_DT_[Charge_].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							CX_DT_[Charge_].Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							CX_DT_[Charge_].E_Add(XY_, (-3. / 2 * Ti[XY_[0]][XY_[1]] * qe - 1. / 2 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
						// std::cout << Index_ << endl;
						setfate(0, 1, 15); // CX with T
						// std::cout << name_ + " before CXDT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(7, 2);
						if (K_mu == 2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_1[0] + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2];
							CX_DT_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								CX_DT_[Charge_].E_Add(XY_, -(1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								CX_DT_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						if (K_mu == 3)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_1[0] + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2] - ua_T_1[XY_[0]][XY_[1]];
							CX_DT_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								CX_DT_[Charge_].E_Add(XY_, -(1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								CX_DT_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						P = &T;
						P->Particlefrom(this, 1, 0);
						P->CalTn();
						// std::cout << name_ + " after CXDT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						//  std::cout << P->Tn() << '\t' << P->name() << endl;
						return;
					}
					else if (this == &T)
					{
						CX_DT_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_DT_[Charge_].Mu_Add(XY_, -Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_DT_[Charge_].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							CX_DT_[Charge_].Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							CX_DT_[Charge_].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
						// std::cout << Index_ << endl;
						setfate(0, 1, 3); // CX with D
						// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(7, 2);
						if (K_mu == 2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_1[0] + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2];
							CX_DT_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								CX_DT_[Charge_].E_Add(XY_, -(1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								CX_DT_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (K_mu == 3)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_1[0] + B[XY_[0]][XY_[1]][1] * V_1[1] + B[XY_[0]][XY_[1]][2] * V_1[2] - V_temp1;
							CX_DT_[Charge_].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								CX_DT_[Charge_].E_Add(XY_, -(1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								CX_DT_[Charge_].E_Add(XY_, (1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						P = &D;
						P->Particlefrom(this, 1, 0);
						P->CalTn();
						// std::cout << name_ + " after CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						//  std::cout << P->Tn() << '\t' << P->name() << endl;
						return;
					}
				}
				else if (MeshMode == 3)
				{
					CX_[Charge_].n_Add(Tri_Index_, Weight_ * NumPar_now);
					if (K_mu == 3)
					{
						if (this == &D)
						{
							CX_[Charge_].Mu_Add(Tri_Index_, Weight_ * NumPar_now * (Dmass * V_temp1 - Tmass * (V_T_1_now[0] * B[XY_[0]][XY_[1]][0] + V_T_1_now[1] * B[XY_[0]][XY_[1]][1] + V_T_1_now[2] * B[XY_[0]][XY_[1]][2])));
							CX_[Charge_].E_Add(Tri_Index_, 0.5 * Weight_ * NumPar_now * (Dmass * (V_[1] * V_[1] + V_[1] * V_[1] + V_[2] * V_[2]) - Tmass * (V_T_1_now[1] * V_T_1_now[1] + V_T_1_now[1] * V_T_1_now[1] + V_T_1_now[2] * V_T_1_now[2])));
						}
						if (this == &T)
						{
							CX_[Charge_].Mu_Add(Tri_Index_, Weight_ * NumPar_now * (Tmass * V_temp1 - Dmass * (V_D_1_now[0] * B[XY_[0]][XY_[1]][0] + V_D_1_now[1] * B[XY_[0]][XY_[1]][1] + V_D_1_now[2] * B[XY_[0]][XY_[1]][2])));
							CX_[Charge_].E_Add(Tri_Index_, 0.5 * Weight_ * NumPar_now * (Tmass * (V_[1] * V_[1] + V_[1] * V_[1] + V_[2] * V_[2]) - Dmass * (V_D_1_now[1] * V_D_1_now[1] + V_D_1_now[1] * V_D_1_now[1] + V_D_1_now[2] * V_D_1_now[2])));
						}
					}
					for (int i = 0; i < 3; i++)
						SetV(i, V_T_1_now[i]);
					setfate(0, 1, Index_); // CX with D
				}
			}
			else if (rand_one < (Ion_factor + CX_factor + CX_DT_factor + R_with_H_factor) / sigma_all)
			{
				if (MeshMode == 3)
				{
					if (this == &D)
					{
						R_with_H_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
						R_with_H_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * Dmass * V_[0], -Weight_ * NumPar_now * Dmass * V_[0], -Weight_ * NumPar_now * Dmass * V_[0]);
						R_with_H_[0].E_Add(Tri_Index_, -(0.5 * Dmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						setfate(0, 9, 3); // CX with D
						// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(9, 1);
						R_with_H_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * Dmass * V_[0], +Weight_ * NumPar_now * Dmass * V_[0], +Weight_ * NumPar_now * Dmass * V_[0]);
						R_with_H_[0].E_Add(Tri_Index_, +(0.5 * Dmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
					}
					if (this == &T)
					{
						R_with_H_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
						R_with_H_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * Tmass * V_[0], -Weight_ * NumPar_now * Tmass * V_[0], -Weight_ * NumPar_now * Tmass * V_[0]);
						R_with_H_[0].E_Add(Tri_Index_, -(0.5 * Tmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						setfate(0, 9, 15); // CX with D
						// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(9, 1);
						R_with_H_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * Tmass * V_[0], +Weight_ * NumPar_now * Tmass * V_[0], +Weight_ * NumPar_now * Tmass * V_[0]);
						R_with_H_[0].E_Add(Tri_Index_, +(0.5 * Tmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
					}
				}
			}
			else if (rand_one < (Ion_factor + CX_factor + CX_DT_factor + R_with_H_factor + R_with_H2_factor) / sigma_all)
			{
				if (MeshMode == 3)
				{
					if (this == &D)
					{
						R_with_H2_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
						R_with_H2_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * D2mass * V_[0], -Weight_ * NumPar_now * D2mass * V_[0], -Weight_ * NumPar_now * D2mass * V_[0]);
						R_with_H2_[0].E_Add(Tri_Index_, -(0.5 * D2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						setfate(0, 9, 4); // CX with D
						// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(9, 2);
						R_with_H2_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * D2mass * V_[0], +Weight_ * NumPar_now * D2mass * V_[0], +Weight_ * NumPar_now * D2mass * V_[0]);
						R_with_H2_[0].E_Add(Tri_Index_, +(0.5 * D2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
					}
					if (this == &T)
					{
						R_with_H2_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
						R_with_H2_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * T2mass * V_[0], -Weight_ * NumPar_now * T2mass * V_[0], -Weight_ * NumPar_now * T2mass * V_[0]);
						R_with_H2_[0].E_Add(Tri_Index_, -(0.5 * T2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						setfate(0, 9, 16); // CX with D
						// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						Init(9, 2);
						R_with_H2_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * T2mass * V_[0], +Weight_ * NumPar_now * T2mass * V_[0], +Weight_ * NumPar_now * T2mass * V_[0]);
						R_with_H2_[0].E_Add(Tri_Index_, +(0.5 * T2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
					}
				}
			}
		}
	}
	else if (this == &H2 || this == &D2 || this == &T2)
	{
		if (MeshMode == 1)
		{
			Ion_factor = Ion_[0].cs(XY_, Zone_);
			if (K_Tn == 1)
			{
				CX_factor = CX_[0].cs_now();
				Ela_factor = Ela_[0].cs_now();
				CX_DT_factor = CX_DT_[0].cs_now();
				Ela_DT_factor = Ela_DT_[0].cs_now();
			}
			if (K_Tn == 2)
			{
				CX_factor = CX_[0].cs(XY_, Zone_);
				Ela_factor = Ela_[0].cs(XY_, Zone_);
				CX_DT_factor = CX_DT_[0].cs(XY_, Zone_);
				Ela_DT_factor = Ela_DT_[0].cs(XY_, Zone_);
			}
			MAR_factor = MAR_[0].cs(XY_[0], XY_[1]);
			Diss1_factor = Diss1_[0].cs(XY_[0], XY_[1]);
			Diss2_factor = Diss2_[0].cs(XY_[0], XY_[1]);
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].cs(XY_[0], XY_[1]);
		}
		else if (MeshMode == 3)
		{
			if (Zone_ < 6)
			{
				Ion_factor = Ion_[0].cs(Tri_Index_);
				if (K_Tn == 1)
				{
					CX_factor = CX_[0].cs_now();
					Ela_factor = Ela_[0].cs_now();
					if (K_DT)
					{
						CX_DT_factor = CX_DT_[0].cs_now();
						Ela_DT_factor = Ela_DT_[0].cs_now();
					}
					else
					{
						CX_DT_factor = 0.;
						Ela_DT_factor = 0.;
					}
				}
				if (K_Tn == 2)
				{
					CX_factor = CX_[0].cs(Tri_Index_);
					Ela_factor = Ela_[0].cs(Tri_Index_);
					if (K_DT)
					{
						CX_DT_factor = CX_DT_[0].cs(Tri_Index_);
						Ela_DT_factor = Ela_DT_[0].cs(Tri_Index_);
					}
					else
					{
						CX_DT_factor = 0.;
						Ela_DT_factor = 0.;
					}
				}
				MAR_factor = MAR_[0].cs(Tri_Index_);
				Diss1_factor = Diss1_[0].cs(Tri_Index_);
				Diss2_factor = Diss2_[0].cs(Tri_Index_);
				for (int i = 0; i < 3; i++)
					DS_factor[i] = DS_[1][i].cs(Tri_Index_);
			}
			else
			{
				Ion_factor = 0;
				CX_factor = 0;
				Ela_factor = 0;
				CX_DT_factor = 0;
				Ela_DT_factor = 0;
				MAR_factor = 0;
				Diss1_factor = 0;
				Diss2_factor = 0;
			}
			R_with_H_factor = R_with_H_[0].cs(Tri_Index_);
			R_with_H2_factor = R_with_H2_[0].cs(Tri_Index_);
		}
		double T_D_0_temp;
		double Sum_factor[10] = {0.}, rand_one = Random();
		int Num_Collision_D2 = 10;
		Sum_factor[0] = Ion_factor;
		Sum_factor[1] = Sum_factor[0] + 0.; // MAR_factor;
		Sum_factor[2] = Sum_factor[1] + Diss1_factor;
		Sum_factor[3] = Sum_factor[2] + Diss2_factor;
		Sum_factor[4] = Sum_factor[3] + CX_factor;
		Sum_factor[5] = Sum_factor[4] + Ela_factor;
		Sum_factor[6] = Sum_factor[5] + CX_DT_factor;
		Sum_factor[7] = Sum_factor[6] + Ela_DT_factor;
		Sum_factor[8] = Sum_factor[7] + R_with_H_factor;
		Sum_factor[9] = Sum_factor[8] + R_with_H2_factor;

		if (Charge_ == 0)
		{
			double rand_temp, Var_temp;
			if (K_flight == 1)
			{
				rand_temp = 0;
				Var_temp = 1;
			}
			else if (K_flight == 2)
			{
				rand_temp = Random();
				Var_temp = CollProb();
			}
			else if (K_flight == 3)
			{
				rand_temp = Random();
				if (K_Tn == 1)
				{
					Var_temp = lambda_min_[0] / lambda_now_;
				}
				else if (K_Tn == 2)
				{
					Var_temp = lambda_min_[0] / lambda_[XY_[0]][XY_[1]][0];
				}
				// std::cout << "CollProb: " << rand_temp << " " << Var_temp << endl;
			}
			if (rand_temp < Var_temp)
			{
				int i = 0;
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ion_[0].n_Add(XY_, Weight_ * NumPar_now);
					// Init(); // electron collision
					setfate(0, 5, Index_);
					T_D_0_temp = 0;
					if (!K_D2Flight)
					{
						double beilv = 0 + H2PCollCal();
						if (beilv == 0)
						{
							Weight_ = 0;
							return;
						}

						if (this == &H2)
							P = &H;
						else if (this == &D2)
							P = &D;
						else if (this == &T2)
							P = &T;

						Vector3 v_H2(V_[0], V_[0], V_[0]);
						double E_FrankCondon = 3.0;											// 3 eV
						Vector3 v_H = calculate_dissociation_velocity(v_H2, E_FrankCondon); // 模拟一次分裂
						V_[0] = v_H.x;
						V_[1] = v_H.y;
						V_[2] = v_H.z;
						P->Particlefrom(this, beilv, 0);
						return;
					}
					else
						Charge_ = 1;

					if (K_test1 || K_test2)
					{
						std::cout << "coll_1 Ion" << endl;
					}
				}
				else if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					MAR_[0].n_Add(XY_, Weight_ * NumPar_now);
					if (this == &H2)
						P = &H;
					else if (this == &D2)
						P = &D;
					else if (this == &T2)
						P = &T;
					P->Particlefrom(this, 3, 0);
					P->setfate(0, 6, this->Index_);
					// P->Init(6);

					if (K_test1)
					{
						std::cout << "coll_2 Mar" << endl;
					}
					// PartoPar_Temp.store(this, 3, &D, 0);
					return;
				}
				else if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss1_[0].n_Add(XY_, Weight_ * NumPar_now);
					if (this == &H2)
						P = &H;
					else if (this == &D2)
						P = &D;
					else if (this == &T2)
						P = &T;

					Vector3 v_H2(V_[0], V_[0], V_[0]);									// 假设分子沿 X 轴运动
					double E_FrankCondon = 3.0;											// 3 eV
					Vector3 v_H = calculate_dissociation_velocity(v_H2, E_FrankCondon); // 模拟一次分裂
					V_[0] = v_H.x;
					V_[1] = v_H.y;
					V_[2] = v_H.z;
					P->Particlefrom(this, 2, 0);
					P->setfate(0, 7, this->Index_);
					// P->Init(6);

					if (K_test1 || K_test2)
					{
						std::cout << "coll_3 Diss" << endl;
					}
					// PartoPar_Temp.store(this, 2, &D, 0);
					return;
				}
				else if (rand_one < Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					V_temp2 = B[XY_[0]][XY_[1]][0] * V_[0] + B[XY_[0]][XY_[1]][1] * V_[1] + B[XY_[0]][XY_[1]][2] * V_[2];
					Diss2_[0].n_Add(XY_, Weight_ * NumPar_now);
					if (K_mu == 1)
					{
						Diss2_[0].Mu_Add(XY_, 0.5 * Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
						Diss2_[0].E_Add(XY_, 0.5 * (1.5 * Ti[XY_[0]][XY_[1]] * qe + 0.5 * Dmass * ua_D_1[XY_[0]][XY_[1]]) * Weight_ * NumPar_now);
					}
					else if (K_mu == 2 || K_mu == 3)
						if (this == &H2)
						{
							Diss2_[0].Mu_Add(XY_, 0.5 * Hmass * V_temp2 * Weight_ * NumPar_now);
							Diss2_[0].E_Add(XY_, 0.5 * (0.5 * Hmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2])) * Weight_ * NumPar_now);
							P = &H;
						}
						else if (this == &D2)
						{
							Diss2_[0].Mu_Add(XY_, 0.5 * Dmass * V_temp2 * Weight_ * NumPar_now);
							Diss2_[0].E_Add(XY_, 0.5 * (0.5 * Dmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2])) * Weight_ * NumPar_now);
							P = &D;
						}
						else if (this == &T2)
						{
							Diss2_[0].Mu_Add(XY_, 0.5 * Tmass * V_temp2 * Weight_ * NumPar_now);
							Diss2_[0].E_Add(XY_, 0.5 * (0.5 * Tmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2])) * Weight_ * NumPar_now);
							P = &T;
						}

					Vector3 v_H2(V_[0], V_[0], V_[0]);									// 假设分子沿 X 轴运动
					double E_FrankCondon = 3.0;											// 3 eV
					Vector3 v_H = calculate_dissociation_velocity(v_H2, E_FrankCondon); // 模拟一次分裂
					V_[0] = v_H.x;
					V_[1] = v_H.y;
					V_[2] = v_H.z;
					P->Particlefrom(this, 1, 0);
					P->CalTn();
					P->setfate(0, 8, this->Index_);
					// P->Init(6);
					if (K_test1)
					{
						std::cout << "coll_4 Diss2" << endl;
					}
					// PartoPar_Temp.store(this, 1, &D, 0);
					// Weight_ = 0;
					return;
				}
				else if (rand_one < Sum_factor[4] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_[0].n_Add(XY_, Weight_ * NumPar_now);
					if (this == &D2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_[0].Mu_Add(XY_, -Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2 || K_mu == 3)
						{
							CX_[0].Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							CX_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
					}
					else if (this == &T2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_[0].Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2 || K_mu == 3)
						{
							CX_[0].Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							CX_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
					}
					setfate(0, 1, Index_);
					// std::cout << name_ + " before CX V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					// Init(2, 1);
					CalTn();

					if (!K_D2Flight)
					{
						double beilv = 1 + H2PCollCal();

						if (this == &H2)
							P = &H;
						else if (this == &D2)
							P = &D;
						else if (this == &T2)
							P = &T;

						if (Random() * beilv < 1)
						{
							if (this == &H2)
							{
								V_[0] = V_H_1_now[0];
								V_[1] = V_H_1_now[1];
								V_[2] = V_H_1_now[2];
								Tn_ = Ti[XY_[0]][XY_[1]];
							}
							if (this == &D2)
							{
								V_[0] = V_D_1_now[0];
								V_[1] = V_D_1_now[1];
								V_[2] = V_D_1_now[2];
								Tn_ = Ti[XY_[0]][XY_[1]];
							}
							if (this == &T2)
							{
								V_[0] = V_T_1_now[0];
								V_[1] = V_T_1_now[1];
								V_[2] = V_T_1_now[2];
								Tn_ = Ti[XY_[0]][XY_[1]];
							}
						}
						else
						{
							Vector3 v_H2(V_[0], V_[0], V_[0]);									// 假设分子沿 X 轴运动
							double E_FrankCondon = 3.0;											// 3 eV
							Vector3 v_H = calculate_dissociation_velocity(v_H2, E_FrankCondon); // 模拟一次分裂
							V_[0] = v_H.x;
							V_[1] = v_H.y;
							V_[2] = v_H.z;
						}
						P->Particlefrom(this, beilv, 0);
						P->CalTn();
						return;
					}
					else
					{
						Charge_ = 1;
						// std::cout << name_ + " after CX V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
						if (K_test1 || K_test2)
						{
							std::cout << "CX: " << X_[0] << ',' << X_[1] << '\t' << X_new_[0] << ',' << X_new_[1] << endl;
						}
						// std::cout << Index_ << endl;
						if (this == &H2)
							PartoPar_Temp.store(this, 1, &H, V_2, 1, Index_, 0);
						else if (this == &D2)
							PartoPar_Temp.store(this, 1, &D, V_2, 1, Index_, 0);
						else if (this == &T2)
							PartoPar_Temp.store(this, 1, &T, V_2, 1, Index_, 0);
						if (K_test1)
						{
							std::cout << "coll_5 CX" << endl;
						}
					}
				}
				else if (rand_one < Sum_factor[5] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ela_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					if (this == &D2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							Ela_.begin()->Mu_Add(XY_, -Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							Ela_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							Ela_.begin()->Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							Ela_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
					}
					else if (this == &T2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							Ela_.begin()->Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							Ela_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							Ela_.begin()->Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							Ela_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
					}
					setfate(0, 9, Index_);
					// std::cout << name_ + " before Ela V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					Init(2, 1);
					if (K_mu == 2)
					{
						V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2];
						if (this == &D2)
						{
							Ela_[Charge_].Mu_Add(XY_, Dmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								Ela_[Charge_].E_Add(XY_, -(1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_[Charge_].E_Add(XY_, (1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (this == &T2)
						{
							Ela_[Charge_].Mu_Add(XY_, Tmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 < 0)
								Ela_[Charge_].E_Add(XY_, -(1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_[Charge_].E_Add(XY_, (1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
					}
					else if (K_mu == 3)
					{
						if (this == &D2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2] - V_temp1;
							Ela_[Charge_].Mu_Add(XY_, Dmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 * ua_D_1[XY_[0]][XY_[1]] < 0)
								Ela_[Charge_].E_Add(XY_, -(1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_[Charge_].E_Add(XY_, (1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (this == &T2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2] - V_temp3;
							Ela_[Charge_].Mu_Add(XY_, Tmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 * ua_T_1[XY_[0]][XY_[1]] < 0)
								Ela_[Charge_].E_Add(XY_, -(1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_[Charge_].E_Add(XY_, (1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
					}
					CalTn();
					// std::cout << name_ + " after Ela V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					if (K_test1 || K_test2)
					{
						std::cout << "coll_6 ela" << endl;
					}
					return;
				}
				else if (rand_one < Sum_factor[6] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_DT_[0].n_Add(XY_, Weight_ * NumPar_now);
					if (this == &D2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_DT_[0].Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_DT_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2 || K_mu == 3)
						{
							CX_DT_[0].Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							CX_DT_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
					}
					else if (this == &T2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							CX_DT_[0].Mu_Add(XY_, -Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							CX_DT_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2 || K_mu == 3)
						{
							CX_DT_[0].Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							CX_DT_[0].E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
					}
					// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					setfate(0, 1, Index_);
					Init(7, 1);
					CalTn();
					Charge_ = 1;
					// std::cout << name_ + " after CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					if (K_test1 || K_test2)
					{
						std::cout << "CX: " << X_[0] << ',' << X_[1] << '\t' << X_new_[0] << ',' << X_new_[1] << endl;
					}
					// std::cout << Index_ << endl;
					if (this == &H2)
						PartoPar_Temp.store(this, 1, &H, V_2, 1, Index_, 0);
					else if (this == &D2)
						PartoPar_Temp.store(this, 1, &T, V_2, 1, Index_, 0);
					else if (this == &T2)
						PartoPar_Temp.store(this, 1, &D, V_2, 1, Index_, 0);
					if (K_test1)
					{
						std::cout << "coll_7 CX_DT" << endl;
					}
				}
				else if (rand_one < Sum_factor[7] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ela_DT_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					if (this == &D2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							Ela_DT_.begin()->Mu_Add(XY_, -Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							Ela_DT_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(ua_T_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							Ela_DT_.begin()->Mu_Add(XY_, -Tmass * V_temp3 * Weight_ * NumPar_now);
							Ela_DT_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Tmass * pow(V_temp3, 2)) * Weight_ * NumPar_now);
						}
					}
					else if (this == &T2)
					{
						if (K_mu == 1 || K_Vi == 2)
						{
							Ela_DT_.begin()->Mu_Add(XY_, -Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							Ela_DT_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						else if (K_mu == 2)
						{
							Ela_DT_.begin()->Mu_Add(XY_, -Dmass * V_temp1 * Weight_ * NumPar_now);
							Ela_DT_.begin()->E_Add(XY_, -(3. / 2 * Ti[XY_[0]][XY_[1]] * qe + 1. / 2 * Dmass * pow(V_temp1, 2)) * Weight_ * NumPar_now);
						}
					}
					setfate(0, 9, Index_);
					// std::cout << name_ + " before Ela_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					Init(7, 1);
					if (K_mu == 2)
					{
						V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2];
						if (this == &D2)
						{
							Ela_DT_[Charge_].Mu_Add(XY_, Tmass * V_temp2 * Weight_ * NumPar_now);
							Ela_DT_[Charge_].E_Add(XY_, (1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (this == &T2)
						{
							Ela_DT_[Charge_].Mu_Add(XY_, Dmass * V_temp2 * Weight_ * NumPar_now);
							Ela_DT_[Charge_].E_Add(XY_, (1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
					}
					if (K_mu == 3)
					{
						if (this == &D2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2] - V_temp3;
							Ela_DT_[Charge_].Mu_Add(XY_, Tmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 * ua_T_1[XY_[0]][XY_[1]] < 0)
								Ela_DT_[Charge_].E_Add(XY_, -(1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_DT_[Charge_].E_Add(XY_, (1. / 2 * Tmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
						else if (this == &T2)
						{
							V_temp2 = B[XY_[0]][XY_[1]][0] * V_2[0] + B[XY_[0]][XY_[1]][1] * V_2[1] + B[XY_[0]][XY_[1]][2] * V_2[2] - V_temp1;
							Ela_DT_[Charge_].Mu_Add(XY_, Dmass * V_temp2 * Weight_ * NumPar_now);
							if (V_temp2 * ua_D_1[XY_[0]][XY_[1]] < 0)
								Ela_DT_[Charge_].E_Add(XY_, -(1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							else
								Ela_DT_[Charge_].E_Add(XY_, (1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
						}
					}
					// std::cout << name_ + " after Ela_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
					CalTn();
					if (K_test1 || K_test2)
					{
						std::cout << "coll_8 ela_DT" << endl;
					}
					return;
				}
				else if (rand_one < Sum_factor[8] / Sum_factor[Num_Collision_D2 - 1])
				{
					if (MeshMode == 3)
					{
						if (this == &D)
						{
							R_with_H_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
							R_with_H_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * Dmass * V_[0], -Weight_ * NumPar_now * Dmass * V_[0], -Weight_ * NumPar_now * Dmass * V_[0]);
							R_with_H_[0].E_Add(Tri_Index_, -(0.5 * Dmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
							setfate(0, 9, 3); // CX with D
							// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
							Init(9, 1);
							R_with_H_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * Dmass * V_[0], +Weight_ * NumPar_now * Dmass * V_[0], +Weight_ * NumPar_now * Dmass * V_[0]);
							R_with_H_[0].E_Add(Tri_Index_, +(0.5 * Dmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						}
						if (this == &T)
						{
							R_with_H_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
							R_with_H_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * Tmass * V_[0], -Weight_ * NumPar_now * Tmass * V_[0], -Weight_ * NumPar_now * Tmass * V_[0]);
							R_with_H_[0].E_Add(Tri_Index_, -(0.5 * Tmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
							setfate(0, 9, 15); // CX with D
							// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
							Init(9, 1);
							R_with_H_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * Tmass * V_[0], +Weight_ * NumPar_now * Tmass * V_[0], +Weight_ * NumPar_now * Tmass * V_[0]);
							R_with_H_[0].E_Add(Tri_Index_, +(0.5 * Tmass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						}
					}
				}
				else if (rand_one < Sum_factor[9] / Sum_factor[Num_Collision_D2 - 1])
				{
					if (MeshMode == 3)
					{
						if (this == &D)
						{
							R_with_H2_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
							R_with_H2_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * D2mass * V_[0], -Weight_ * NumPar_now * D2mass * V_[0], -Weight_ * NumPar_now * D2mass * V_[0]);
							R_with_H2_[0].E_Add(Tri_Index_, -(0.5 * D2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
							setfate(0, 9, 4); // CX with D
							// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
							Init(9, 2);
							R_with_H2_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * D2mass * V_[0], +Weight_ * NumPar_now * D2mass * V_[0], +Weight_ * NumPar_now * D2mass * V_[0]);
							R_with_H2_[0].E_Add(Tri_Index_, +(0.5 * D2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						}
						if (this == &T)
						{
							R_with_H2_[0].n_Add(Tri_Index_, Weight_ * NumPar_now);
							R_with_H2_[0].Mu_Add(Tri_Index_, -Weight_ * NumPar_now * T2mass * V_[0], -Weight_ * NumPar_now * T2mass * V_[0], -Weight_ * NumPar_now * T2mass * V_[0]);
							R_with_H2_[0].E_Add(Tri_Index_, -(0.5 * T2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
							setfate(0, 9, 16); // CX with D
							// std::cout << name_ + " before CX_DT V: " << V_[0] << ", " << V_[1] << ", " << V_[2] << " T: " << Tn_ << endl;
							Init(9, 2);
							R_with_H2_[0].Mu_Add(Tri_Index_, +Weight_ * NumPar_now * T2mass * V_[0], +Weight_ * NumPar_now * T2mass * V_[0], +Weight_ * NumPar_now * T2mass * V_[0]);
							R_with_H2_[0].E_Add(Tri_Index_, +(0.5 * T2mass * (V_[0] * V_[0] + V_[1] * V_[1] + V_[2] * V_[2]) * Weight_ * NumPar_now));
						}
					}
				}
			}
		}
		if (Charge_ == 1)
		{
			// std::cout << X_[0] << '\t' << X_[1] << '\t' << name_ << '\t';
			// std::cout << XY_[0] << '\t' << XY_[1] << '\t' << Zone_ << endl;
			Sum_n_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now * dt;
			double Sum_factor[10] = {0.}, rand_one = Random(), rand_temp;
			if (K_flight == 1 || !K_D2Flight || K_flight == 3)
				rand_temp = -1;
			else
				rand_temp = Random();
			int Num_Collision_D2 = 3;
			Sum_factor[0] = DS_factor[0];
			Sum_factor[1] = Sum_factor[0] + DS_factor[1];
			Sum_factor[2] = Sum_factor[1] + DS_factor[2];
			/*for (int i = 0; i < 3; i++)
			{
				std::cout << Sum_factor[i] << " ";
			}*/
			// std::cout << endl;
			// std::cout << CollProb() << endl;
			if (rand_temp < CollProb())
			{
				V_temp2 = B[XY_[0]][XY_[1]][0] * V_[0] + B[XY_[0]][XY_[1]][1] * V_[1] + B[XY_[0]][XY_[1]][2] * V_[2];
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
					if (K_mu == 1)
					{
						if (this == &D2)
						{
							DS_[1][0].Mu_Add(XY_, 2. * Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							DS_[1][0].E_Add(XY_, 2. * (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
						if (this == &T2)
						{
							DS_[1][0].Mu_Add(XY_, 2. * Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							DS_[1][0].E_Add(XY_, 2. * (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
						}
					}
					else if (K_mu == 2 || K_mu == 3)
					{
						DS_[1][0].Mu_Add(XY_, 2. * Dmass * V_temp2 * Weight_ * NumPar_now);
						DS_[1][0].E_Add(XY_, 2 * (3. / 2 * Tn_ * qe + 1. / 2 * mass_ * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
					}
					Charge_ = 0;
					Weight_ = 0;
					if (K_test1 || K_test2)
					{
						std::cout << "coll_7 DS1" << endl;
					}
					return;
				}
				else if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
					if (K_mu == 1)
					{
						if (this == &D2)
						{
							DS_[1][1].Mu_Add(XY_, Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_D_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
						}
						if (this == &T2)
						{
							DS_[1][1].Mu_Add(XY_, Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
							DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_T_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
						}
					}
					else if (K_mu == 2 || K_mu == 3)
						if (this == &H2)
						{
							DS_[1][1].Mu_Add(XY_, Hmass * V_temp2 * Weight_ * NumPar_now);
							DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							P = &H;
						}
						else if (this == &D2)
						{
							DS_[1][1].Mu_Add(XY_, Dmass * V_temp2 * Weight_ * NumPar_now);
							DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							P = &D;
						}
						else if (this == &T2)
						{
							DS_[1][1].Mu_Add(XY_, Tmass * V_temp2 * Weight_ * NumPar_now);
							DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(V_temp2, 2))) * Weight_ * NumPar_now);
							P = &T;
						}
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 11, Index_);
					P->Init(6);
					Charge_ = 0;
					Weight_ = 0;
					if (K_test1 || K_test2)
					{
						std::cout << "coll_8 DS2" << endl;
					}
					return;
				}
				else if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
					if (this == &H2)
						P = &H;
					if (this == &D2)
						P = &D;
					if (this == &T2)
						P = &T;
					P->Particlefrom(this, 2, 0);
					P->setfate(0, 12, Index_);
					P->Init(6);
					if (K_test1 || K_test2)
					{
						std::cout << "coll_9 DS3" << endl;
					}
					return;
				}
			}
			if (!K_D2Flight)
			{
				Charge_ = 0;
				return;
			}
			/*else
			{
					if (PartoPar_Temp.ifstore())
							PartoPar_Temp.store(this, 1, this, Charge_);
			}*/
			// if (PartoPar_Temp.ifstore() == 0)
		}
	}
	else if (this == &CD4)
	{
		Ion_factor = Ion_.begin()->cs(XY_, Zone_);
		CX_factor = CX_.begin()->cs(XY_, Zone_);
		if (Zone_ == 1)
		{
			Diss1_factor = Diss1_[0].Cor_cs();
			Diss2_factor = Diss2_[0].Cor_cs();
			for (int i = 0; i < 4; i++)
				DS_factor[i] = DS_[1][i].Cor_cs();
		}
		else if (Zone_ == 6)
		{
			Diss1_factor = 0;
			Diss2_factor = 0;
			for (int i = 0; i < 4; i++)
				DS_factor[i] = 0;
		}
		else
		{
			Diss1_factor = Diss1_[0].cs(XY_[0], XY_[1]);
			Diss2_factor = Diss2_[0].cs(XY_[0], XY_[1]);
			for (int i = 0; i < 4; i++)
				DS_factor[i] = DS_[1][i].cs(XY_[0], XY_[1]);
		}
		if (Charge_ == 0)
		{
			double Sum_factor[4] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 4;
			Sum_factor[0] = Ion_factor;
			Sum_factor[1] = Sum_factor[0] + CX_factor;
			Sum_factor[2] = Sum_factor[1] + Diss1_factor;
			Sum_factor[3] = Sum_factor[2] + Diss2_factor;
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ion_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Charge_ = 1;
					setfate(0, 5, 5);
					return;
				}
				else if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Init(2, 1);
					Charge_ = 1;
					setfate(0, 1, 5);
					CalTn();
					PartoPar_Temp.store(this, 1, &D, V_2, 1, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 1, &H, V_2, 1, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 1, &D, V_2, 1, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 1, &T, V_2, 1, 5, 0);*/
					return;
				}
				if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss1_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 7, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 1, &H, V_, 7, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 1, &D, V_, 7, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 1, &T, V_, 7, 5, 0);*/
					P = &CD3;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 7, 5);
					P->Init(6);
					return;
				}
				if (rand_one < Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss2_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 14, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 1, &H, V_, 14, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 1, &D, V_, 14, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 1, &T, V_, 14, 5, 0);*/
					P = &CD3;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 8, 5);
					P->Init(6);
					return;
				}
			}
		}
		else if (Charge_ == 1)
		{
			double Sum_factor[4] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 4;
			Sum_factor[0] = DS_factor[0];
			Sum_factor[1] = DS_factor[1] + Sum_factor[0];
			Sum_factor[2] = DS_factor[2] + Sum_factor[1];
			Sum_factor[3] = DS_factor[3] + Sum_factor[2];
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
					P = &CD3;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 10, 5);
					P->Init(6);
					return;
				}
				if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 11, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 1, &H, V_, 11, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 1, &D, V_, 11, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 1, &T, V_, 11, 5, 0);*/
					P = &CD3;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 11, 5);
					P->Init(6);
					return;
				}
				if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 12, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 1, &H, V_, 12, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 1, &D, V_, 12, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 1, &T, V_, 12, 5, 0);*/
					P = &CD3;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 12, 5);
					P->Init(6);
					return;
				}
				if (rand_one < Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][3].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 2, &D, V_, 13, 5, 0);
					/*if (K_H)
						PartoPar_Temp.store(this, 2, &H, V_, 13, 5, 0);
					if (K_D)
						PartoPar_Temp.store(this, 2, &D, V_, 13, 5, 0);
					if (K_T)
						PartoPar_Temp.store(this, 2, &T, V_, 13, 5, 0);*/
					P = &CD2;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 13, 5);
					P->Init(6);
					return;
				}
			}
		}
	}
	else if (this == &CD3)
	{
		Ion_factor = Ion_.begin()->cs(XY_, Zone_);
		CX_factor = CX_.begin()->cs(XY_, Zone_);
		if (Zone_ == 1)
		{
			Diss1_factor = Diss1_[0].Cor_cs();
			Diss2_factor = Diss2_[0].Cor_cs();
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].Cor_cs();
		}
		else if (Zone_ == 6)
		{
			Diss1_factor = 0;
			Diss2_factor = 0;
			for (int i = 0; i < 3; i++)
				DS_factor[i] = 0;
		}
		else
		{
			Diss1_factor = Diss1_[0].cs(XY_[0], XY_[1]);
			Diss2_factor = Diss2_[0].cs(XY_[0], XY_[1]);
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].cs(XY_[0], XY_[1]);
		}
		if (Charge_ == 0)
		{
			double Sum_factor[4] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 4;
			Sum_factor[0] = Ion_factor;
			Sum_factor[1] = Sum_factor[0] + CX_factor;
			Sum_factor[2] = Sum_factor[1] + Diss1_factor;
			Sum_factor[3] = Sum_factor[2] + Diss2_factor;
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ion_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Charge_ = 1;
					setfate(0, 5, 6);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Init(2, 1);
					Charge_ = 1;
					VtoVcharge();
					setfate(0, 1, 6);
					CalTn();
					PartoPar_Temp.store(this, 1, &D, V_2, 1, 6, 0);
					/*
					if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_2, 1, 6, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_2, 1, 6, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_2, 1, 6, 0);*/
					return;
				}
				else if (rand_one <
						 Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss1_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 7, 6, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 7, 6, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 7, 6, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 7, 6, 0);*/
					P = &CD2;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 7, 6);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss2_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 8, 6, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 8, 6, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 8, 6, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 8, 6, 0);*/
					P = &CD2;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 8, 6);
					P->Init(6);
					return;
				}
			}
		}
		else if (Charge_ == 1)
		{
			double Sum_factor[3] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 3;
			Sum_factor[0] = DS_factor[0];
			Sum_factor[1] = DS_factor[1] + Sum_factor[0];
			Sum_factor[2] = DS_factor[2] + Sum_factor[1];
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
					P = &CD2;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 10, 6);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 11, 6, 0);
					/*
					if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 11, 6, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 11, 6, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 11, 6, 0);*/
					P = &CD2;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 11, 6);
					P->Init(6);
					return;
				}
				else if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 12, 6, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 12, 6, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 12, 6, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 12, 6, 0);*/
					P = &CD2;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 12, 6);
					P->Init(6);
					return;
				}
			}
		}
	}
	else if (this == &CD2)
	{
		Ion_factor = Ion_.begin()->cs(XY_, Zone_);
		CX_factor = CX_.begin()->cs(XY_, Zone_);
		if (Zone_ == 1)
		{
			Diss1_factor = Diss1_[0].Cor_cs();
			Diss2_factor = Diss2_[0].Cor_cs();
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].Cor_cs();
		}
		else if (Zone_ == 6)
		{
			Diss1_factor = 0;
			Diss2_factor = 0;
			for (int i = 0; i < 3; i++)
				DS_factor[i] = 0;
		}
		else
		{
			Diss1_factor = Diss1_[0].cs(XY_[0], XY_[1]);
			Diss2_factor = Diss2_[0].cs(XY_[0], XY_[1]);
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].cs(XY_[0], XY_[1]);
		}
		if (Charge_ == 0)
		{
			double Sum_factor[4] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 4;
			Sum_factor[0] = Ion_factor;
			Sum_factor[1] = Sum_factor[0] + CX_factor;
			Sum_factor[2] = Sum_factor[1] + Diss1_factor;
			Sum_factor[3] = Sum_factor[2] + Diss2_factor;
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ion_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Charge_ = 1;
					setfate(0, 5, 7);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Init(2, 1);
					Charge_ = 1;
					VtoVcharge();
					setfate(0, 1, 7);
					CalTn();
					PartoPar_Temp.store(this, 1, &D, V_2, 1, 7, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_2, 1, 7, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_2, 1, 7, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_2, 1, 7, 0);*/
					return;
				}
				else if (rand_one <
						 Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss1_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 7, 7, 0);
					/*
					if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 7, 7, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 7, 7, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 7, 7, 0);*/
					P = &CD1;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 7, 7);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss2_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 8, 7, 0);
					/*
					if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 8, 7, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 8, 7, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 8, 7, 0);*/
					P = &CD1;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 8, 7);
					P->Init(6);
					return;
				}
			}
		}
		else if (Charge_ == 1)
		{
			double Sum_factor[3] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 3;
			Sum_factor[0] = DS_factor[0];
			Sum_factor[1] = DS_factor[1] + Sum_factor[0];
			Sum_factor[2] = DS_factor[2] + Sum_factor[1];
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
					P = &CD1;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 10, 7);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 11, 7, 0);
					/*
					if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 11, 7, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 11, 7, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 11, 7, 0);*/
					P = &CD1;
					P->Particlefrom(this, 1, 1);
					P->setfate(0, 11, 7);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 12, 7, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 12, 7, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 12, 7, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 12, 7, 0);*/
					P = &CD1;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 12, 7);
					P->Init(6);
					return;
				}
			}
		}
	}
	else if (this == &CD1)
	{
		Ion_factor = Ion_.begin()->cs(XY_, Zone_);
		CX_factor = CX_.begin()->cs(XY_, Zone_);
		if (Zone_ == 1)
		{
			Diss1_factor = Diss1_[0].Cor_cs();
			Diss2_factor = Diss2_[0].Cor_cs();
			Diss3_factor = Diss3_[0].Cor_cs();
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].Cor_cs();
		}
		else if (Zone_ == 6)
		{
			Diss1_factor = 0;
			Diss2_factor = 0;
			Diss3_factor = 0;
			for (int i = 0; i < 3; i++)
				DS_factor[i] = 0;
		}
		else
		{
			Diss1_factor = Diss1_[0].cs(XY_[0], XY_[1]);
			Diss2_factor = Diss2_[0].cs(XY_[0], XY_[1]);
			Diss3_factor = Diss3_[0].cs(XY_[0], XY_[1]);
			for (int i = 0; i < 3; i++)
				DS_factor[i] = DS_[1][i].cs(XY_[0], XY_[1]);
		}
		if (Charge_ == 0)
		{
			double Sum_factor[5] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 5;
			Sum_factor[0] = Ion_factor;
			Sum_factor[1] = Sum_factor[0] + CX_factor;
			Sum_factor[2] = Sum_factor[1] + Diss1_factor;
			Sum_factor[3] = Sum_factor[2] + Diss2_factor;
			Sum_factor[4] = Sum_factor[3] + Diss3_factor;
			if (Random() < CollProb())
			{
				if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
				{
					Ion_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Charge_ = 1;
					setfate(0, 5, 8);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					CX_.begin()->n_Add(XY_, Weight_ * NumPar_now);
					Init(2, 1);
					Charge_ = 1;
					VtoVcharge();
					setfate(0, 1, 8);
					CalTn();
					PartoPar_Temp.store(this, 1, &D, V_2, 1, 8, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_2, 1, 8, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_2, 1, 8, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_2, 1, 8, 0);*/
					return;
				}
				else if (rand_one <
						 Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss1_[0].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 7, 8, 0);
					/*if (K_back == 1)
						PartoPar_Temp.store(this, 1, &H, V_, 7, 8, 0);
					if (K_back == 2)
						PartoPar_Temp.store(this, 1, &D, V_, 7, 8, 0);
					if (K_back == 3)
						PartoPar_Temp.store(this, 1, &T, V_, 7, 8, 0);*/
					P = &C;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 7, 8);
					P->Init(6);
					return;
				}
				else if (rand_one < Sum_factor[3] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss2_[0].n_Add(XY_, Weight_ * NumPar_now);
					P = &D;
					/*if (K_back == 1)
						P = &H;
					if (K_back == 2)
						P = &D;
					if (K_back == 3)
						P = &T;*/
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 14, 8);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[4] / Sum_factor[Num_Collision_D2 - 1])
				{
					Diss3_[0].n_Add(XY_, Weight_ * NumPar_now);
					P = &C;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 8, 8);
					P->Init(6);
					return;
				}
			}
		}
		else if (Charge_ == 1)
		{
			double Sum_factor[3] = {0.}, rand_one = Random();
			int Num_Collision_D2 = 3;
			Sum_factor[0] = DS_factor[0];
			Sum_factor[1] = DS_factor[1] + Sum_factor[0];
			Sum_factor[2] = DS_factor[2] + Sum_factor[1];
			if (Random() < CollProb())
			{
				if (rand_one < CollProb())
				{
					DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
					P = &C;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 10, 8);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
					P = &D;
					/*if (K_back == 1)
						P = &H;
					if (K_back == 2)
						P = &D;
					if (K_back == 3)
						P = &T;*/
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 11, 8);
					P->Init(6);
					return;
				}
				else if (rand_one <
						 Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
				{
					DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
					PartoPar_Temp.store(this, 1, &D, V_, 12, 8, 0);
					P = &C;
					P->Particlefrom(this, 1, 0);
					P->setfate(0, 12, 8);
					P->Init(6);
					return;
				}
			}
		}
	}
	else
	{
		if (Charge_ == MaxCharge_)
		{
			Ion_factor = 0;
			CX_factor = 0;
			Rec_factor = Rec_[Charge_ - 1].cs(XY_, Zone_);
		}
		else if (Charge_ == 0)
		{
			Ion_factor = Ion_[Charge_].cs(XY_, Zone_);
			CX_factor = CX_[Charge_].cs(XY_, Zone_);
			Rec_factor = 0;
		}
		else
		{
			Ion_factor = Ion_[Charge_].cs(XY_, Zone_);
			CX_factor = CX_[Charge_].cs(XY_, Zone_);
			Rec_factor = Rec_[Charge_ - 1].cs(XY_, Zone_);
		}
		double Sum_factor[3] = {0.}, rand_one = Random();
		int Num_Collision_D2 = 3;
		Sum_factor[0] = Ion_factor;
		Sum_factor[1] = Sum_factor[0] + CX_factor;
		Sum_factor[2] = Sum_factor[1] + Rec_factor;
		if (Random() < CollProb())
		{
			if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
			{
				Ion_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
				Charge_ += 1;
				setfate(0, 7, Index_);
				return;
			}
			else if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
			{
				CX_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
				Init(2, 1);
				Charge_ += 1;
				setfate(0, 1, Index_);
				CalTn();
				PartoPar_Temp.store(this, 1, &D, V_2, 1, Index_, 0);
				/*if (K_back == 1)
					PartoPar_Temp.store(this, 1, &H, V_2, 1, Index_, 0);
				if (K_back == 2)
					PartoPar_Temp.store(this, 1, &D, V_2, 1, Index_, 0);
				if (K_back == 3)
					PartoPar_Temp.store(this, 1, &T, V_2, 1, Index_, 0);*/
				return;
			}
			else if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
			{
				Rec_[Charge_].n_Add(XY_, Weight_ * NumPar_now);
				Charge_ -= 1;
				setfate(0, 4, Index_);
				return;
			}
		}
	}
}

int Particle::H2PCollCal()
{
	// std::cout << X_[0] << '\t' << X_[1] << '\t' << name_ << '\t';
	// std::cout << XY_[0] << '\t' << XY_[1] << '\t' << Zone_ << endl;
	Sum_n_[XY_[0]][XY_[1]][Charge_] += Weight_ * NumPar_now * dt_;
	double Sum_factor[10] = {0.}, rand_one = Random(), rand_temp;
	double E_new;
	int Num_Collision_D2 = 3;
	Sum_factor[0] = DS_factor[0];
	Sum_factor[1] = Sum_factor[0] + DS_factor[1];
	Sum_factor[2] = Sum_factor[1] + DS_factor[2];
	/*for (int i = 0; i < 3; i++)
	{
		std::cout << Sum_factor[i] << " ";
	}*/
	// std::cout << endl;
	// std::cout << CollProb() << endl;
	// std::cout << "Zone: " << Zone_ << endl;
	double V_temp2;
	if (Zone_ < 6)
		V_temp2 = B[XY_[0]][XY_[1]][0] * V_[0] + B[XY_[0]][XY_[1]][1] * V_[1] + B[XY_[0]][XY_[1]][2] * V_[2];
	else
		V_temp2 = 0;

	if (rand_one < Sum_factor[0] / Sum_factor[Num_Collision_D2 - 1])
	{
		DS_[1][0].n_Add(XY_, Weight_ * NumPar_now);
		if (K_mu == 1)
		{
			if (this == &D2)
			{
				DS_[1][0].Mu_Add(XY_, 2. * Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 2. * (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
			}
			if (this == &T2)
			{
				DS_[1][0].Mu_Add(XY_, 2. * Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 2. * (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * pow(ua_D_1[XY_[0]][XY_[1]], 2)) * Weight_ * NumPar_now);
			}
		}
		else if (K_mu == 2 || K_mu == 3)
		{
			if (this == &H2)
			{
				DS_[1][0].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &D2)
			{
				DS_[1][0].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &T2)
			{
				DS_[1][0].Mu_Add(XY_, mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
		}
		if (K_test1 || K_test2)
		{
			std::cout << "coll_7 DS1" << endl;
		}

		return 0;
	}
	else if (rand_one < Sum_factor[1] / Sum_factor[Num_Collision_D2 - 1])
	{
		DS_[1][1].n_Add(XY_, Weight_ * NumPar_now);
		if (K_mu == 1)
		{
			if (this == &D2)
			{
				DS_[1][1].Mu_Add(XY_, Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_D_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
			}
			if (this == &T2)
			{
				DS_[1][1].Mu_Add(XY_, Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_T_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
			}
		}
		else if (K_mu == 2 || K_mu == 3)
		{
			if (this == &H2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &D2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &T2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
		}
		return 1;
	}
	else if (rand_one < Sum_factor[2] / Sum_factor[Num_Collision_D2 - 1])
	{
		DS_[1][2].n_Add(XY_, Weight_ * NumPar_now);
		if (K_mu == 1)
		{
			if (this == &D2)
			{
				DS_[1][1].Mu_Add(XY_, Dmass * ua_D_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_D_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
			}
			if (this == &T2)
			{
				DS_[1][1].Mu_Add(XY_, Tmass * ua_T_1[XY_[0]][XY_[1]] * Weight_ * NumPar_now);
				DS_[1][1].E_Add(XY_, (3. / 2 * Tn_ * qe + 1. / 2 * Dmass * (pow(ua_T_1[XY_[0]][XY_[1]], 2))) * Weight_ * NumPar_now);
			}
		}
		else if (K_mu == 2 || K_mu == 3)
		{
			// there is no source to Ion
			/*if (this == &H2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &D2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}
			else if (this == &T2)
			{
				DS_[1][0].Mu_Add(XY_, 0.5 * mass_ * V_temp2 * Weight_ * NumPar_now);
				DS_[1][0].E_Add(XY_, 0.5 * 1. / 2 * mass_ * (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * Weight_ * NumPar_now);
			}*/
		}
		return 2;
	}
	else
	{
		return -1;
	}
}

void Particle::Pump_add(int i, double Sum)
{
	Pump_[i] += Sum;
}

double Particle::Pump(int i)
{
	return Pump_[i];
}

void Particle::Stat(int n)
{
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
		{
			if (this == &H || this == &D || this == &T)
			{
				if (Sum_n_[i][j][0] != 0)
				{
					for (int m = 0; m < 3; m++)
					{
						V_Grid_[i][j][m][0] = Sum_V_[i][j][m][0] / Sum_n_[i][j][0];
						V_Grid_CX_Ion_Af_[i][j][m][0] = Sum_V_CX_Ion_Af_[i][j][m][0] / Sum_n_[i][j][0];
						V_Grid_CX_Ion_Be_[i][j][m][0] = Sum_V_CX_Ion_Be_[i][j][m][0] / Sum_n_[i][j][0];
						V_Grid_CX_Neu_Af_[i][j][m][0] = Sum_V_CX_Neu_Af_[i][j][m][0] / Sum_n_[i][j][0];
						V_Grid_CX_Neu_Be_[i][j][m][0] = Sum_V_CX_Neu_Be_[i][j][m][0] / Sum_n_[i][j][0];
						V_D_1_[i][j][m][0] = Sum_V_D_1_[i][j][m][0] / Num_V_D_1_[i][j][0];
					}
					E_[i][j][0] = Sum_E_[i][j][0] / Sum_n_[i][j][0];
					T_[i][j][0] = 2. / 3. * (E_[i][j][0] - 0.5 * mass_ * (pow(V_Grid_[i][j][0][0], 2) + pow(V_Grid_[i][j][1][0], 2) + pow(V_Grid_[i][j][2][0], 2)) / 1.6e-19);
				}
				n_[i][j][0] = Sum_n_[i][j][0] / Volume[i][j];
				for (int m = 0; m < 4; m++)
					if (n_Flux_Grid_[i][j][m] != 0)
						T_Flux_Grid_[i][j][m] /= n_Flux_Grid_[i][j][m];
				// FT_.normalizeByVolume();
			}
			else
				for (int k = 0; k <= MaxCharge_; k++)
				{
					if (Sum_n_[i][j][k] != 0)
					{
						for (int m = 0; m < 3; m++)
						{
							V_Grid_[i][j][m][k] = Sum_V_[i][j][m][k] / Sum_n_[i][j][k];
						}
						E_[i][j][k] = Sum_E_[i][j][k] / Sum_n_[i][j][k];
						T_[i][j][k] = 2. / 3. * (E_[i][j][k] - 0.5 * mass_ * (pow(V_Grid_[i][j][0][k], 2) + pow(V_Grid_[i][j][1][k], 2) + pow(V_Grid_[i][j][2][k], 2)) / 1.6e-19);
					}
					n_[i][j][k] = Sum_n_[i][j][k] / Volume[i][j];
					for (int m = 0; m < Num_array_ + 1; m++)
					{
						Sum_n_array_[i][j][k][m] /= Volume[i][j];
					}
				}
		}
	for (int k = 0; k <= MaxCharge_; k++)
	{
		if (this == &H || this == &D || this == &T)
		{
			if (k != MaxCharge_)
			{
				if (this == &D)
				{
					CX_[k].Stat(k, n_, &PRC96_D, n_D_1);
					CX_[k].Stat1(k, n_, &PRC96_D, n_D_1);
					if (K_DT)
						CX_DT_[k].Stat(k, n_, &PRC96_D, n_T_1);
				}
				else if (this == &T)
				{
					CX_[k].Stat(k, n_, &PRC96_D, n_T_1);
					CX_[k].Stat1(k, n_, &PRC96_D, n_T_1);
					if (K_DT)
						CX_DT_[k].Stat(k, n_, &PRC96_D, n_D_1);
				}
				if (K_database_Pra == 1)
				{
					Ion_[k].Stat(k, n_, &PLT12_H);
				}
				else
				{
					Ion_[k].StatEIRENE(k, n_, &R2_1_5_H10, 13.6);
				}
			}
			if (k != 0)
			{
				for (int m = 0; m < N_poloidal; m++)
					for (int n = 0; n < N_radial; n++)
					{
						if (this == &D)
							n_[m][n][k] = n_D_1[m][n];
						else if (this == &T)
							n_[m][n][k] = n_T_1[m][n];
					}
				Rec_[k].RecStat(k, n_, &PRB12_H);
			}
		}
		if (this == &H2 || this == &D2 || this == &T2)
		{
			Ion_[k].Stat(0, n_, NULL);
			CX_[k].Stat(0, n_, NULL);
			Ela_[k].Stat(0, n_, NULL);
			MAR_[k].Stat(0, n_, NULL);
			Diss1_[k].Stat(0, n_, NULL);
			Diss2_[k].Stat(0, n_, NULL);
			for (int i = 0; i < 3; i++)
			{
				DS_[1][i].Stat(0, n_, NULL);
			}
			for (int m = 1; m < N_poloidal - 1; m++)
			{
				for (int n = 1; n < N_radial - 1; n++)
				{
					double Var_temp = R2_2_14_H2.cal(ne[m][n], Te[m][n]);
					DS_[1][2].SetPra(m, n, ne[m][n] * n_[m][n][1] * R2_2_14_H8.cal(ne[m][n], Te[m][n], &Var_temp));
				}
			}
		}
		if (this == &C)
		{
			if (k != MaxCharge_)
			{
				Ion_[k].Stat(k, n_, &PLT96_C);
				CX_[k].Stat(k, n_, &PRC96_C);
			}
			if (k != 0)
				Rec_[k].RecStat(k, n_, &PRB96_C);
		}
		if (this == &Ar)
		{
			if (k != MaxCharge_)
			{
				Ion_[k].Stat(k, n_, &PLT89_Ar);
				CX_[k].Stat(k, n_, &PRC89_Ar);
			}
			if (k != 0)
				Rec_[k].RecStat(k, n_, &PRB89_Ar);
		}
		if (this == &CD4)
		{
			if (k == 0)
			{
				Ion_[k].Stat(k, n_, NULL);
				CX_[k].Stat(k, n_, NULL);
				Diss1_[k].Stat(0, n_, NULL);
				Diss2_[k].Stat(0, n_, NULL);
			}
			else if (Charge_ == 1)
			{
				for (int m = 0; m < 4; m++)
					DS_[k][m].Stat(0, n_, NULL);
			}
		}
		if (this == &CD3)
		{
			if (k == 0)
			{
				Ion_.begin()->Stat(k, n_, NULL);
				CX_.begin()->Stat(k, n_, NULL);
				Diss1_[k].Stat(0, n_, NULL);
				Diss2_[k].Stat(0, n_, NULL);
			}
			else if (Charge_ == 1)
			{
				for (int m = 0; m < 3; m++)
					DS_[k][m].Stat(0, n_, NULL);
			}
		}
		if (this == &CD2)
		{
			if (k == 0)
			{
				Ion_.begin()->Stat(k, n_, NULL);
				CX_.begin()->Stat(k, n_, NULL);
				Diss1_[k].Stat(0, n_, NULL);
				Diss2_[k].Stat(0, n_, NULL);
			}
			else if (Charge_ == 1)
			{
				for (int m = 0; m < 3; m++)
					DS_[k][m].Stat(0, n_, NULL);
			}
		}
		if (this == &CD1)
		{
			if (k == 0)
			{
				Ion_.begin()->Stat(k, n_, NULL);
				CX_.begin()->Stat(k, n_, NULL);
				Diss1_[k].Stat(0, n_, NULL);
				Diss2_[k].Stat(0, n_, NULL);
				Diss3_[k].Stat(0, n_, NULL);
			}
			else if (Charge_ == 1)
			{
				for (int m = 0; m < 3; m++)
					DS_[k][m].Stat(0, n_, NULL);
			}
		}
	}
	for (int i = 0; i < CXtoWall_.size(); i++)
	{
		CXtoWall_.at(i).stat();
	}
	for (int i = 0; i < RectoWall_.size(); i++)
	{
		RectoWall_.at(i).stat();
	}
	for (int i = 0; i < AlltoWall_.size(); i++)
	{
		AlltoWall_.at(i).stat();
	}

	for (int i = 0; i < CXtoTarget_.size(); i++)
	{
		CXtoTarget_.at(i).stat();
	}
	for (int i = 0; i < RectoTarget_.size(); i++)
	{
		RectoTarget_.at(i).stat();
	}
	for (int i = 0; i < AlltoTarget_.size(); i++)
	{
		AlltoTarget_.at(i).stat();
	}

	for (int i = 0; i < CXtoPlasmaBoundary_.size(); i++)
	{
		CXtoPlasmaBoundary_.at(i).stat();
	}
	for (int i = 0; i < RectoPlasmaBoundary_.size(); i++)
	{
		RectoPlasmaBoundary_.at(i).stat();
	}
	for (int i = 0; i < AlltoPlasmaBoundary_.size(); i++)
	{
		AlltoPlasmaBoundary_.at(i).stat();
	}
	/*for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < num_GridBoundry; i++)
			EneFlux_GridBoundry[i][j] /= Flux_GridBoundry[i][j];

		for (int i = 0; i < num_PFRBoundry; i++)
			EneFlux_PFRBoundry[i][j] /= Flux_PFRBoundry[i][j];

		for (int i = 0; i < num_CoreBoundry; i++)
			EneFlux_CoreBoundry[i][j] /= Flux_CoreBoundry[i][j];

		for (int i = 0; i < 36; i++)
		{
			EneFlux_Target_l[i][j] /= Flux_Target_l[i][j];
			EneFlux_Target_r[i][j] /= Flux_Target_r[i][j];
		}
	}*/
}

void Particle::Stat_Tri(int n)
{
	double a, b;
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		a = Tri_B2_[i][0];
		b = Tri_B2_[i][1];
		if (this == &H || this == &D || this == &T)
		{
			if (Tri_Sum_n_[i][0] != 0)
			{
				for (int m = 0; m < 3; m++)
				{
					Tri_V_Grid_[i][m][0] = Tri_Sum_V_[i][m][0] / Tri_Sum_n_[i][0];
				}
				Tri_E_[i][0] = Tri_Sum_E_[i][0] / Tri_Sum_n_[i][0];
				Tri_T_[i][0] = 2. / 3. * (Tri_E_[i][0] - 0.5 * mass_ * (pow(Tri_V_Grid_[i][0][0], 2) + pow(Tri_V_Grid_[i][1][0], 2) + pow(Tri_V_Grid_[i][2][0], 2)) / 1.6e-19);
			}
			Tri_n_[i][0] = Tri_Sum_n_[i][0] / Grid4.triVolume(i);
		}
		else
		{
			for (int k = 0; k < MaxCharge_; k++)
			{
				if (Tri_Sum_n_[i][k] != 0)
				{
					for (int m = 0; m < 3; m++)
					{
						Tri_V_Grid_[i][m][k] = Tri_Sum_V_[i][m][k] / Tri_Sum_n_[i][k];
					}
					Tri_E_[i][k] = Tri_Sum_E_[i][k] / Tri_Sum_n_[i][k];
					Tri_T_[i][k] = 2. / 3. * (Tri_E_[i][k] - 0.5 * mass_ * (pow(Tri_V_Grid_[i][0][k], 2) + pow(Tri_V_Grid_[i][1][k], 2) + pow(Tri_V_Grid_[i][2][k], 2)) / 1.6e-19);
				}
			}
		}
	}
	for (int k = 0; k < MaxCharge_; k++)
	{
		if (this == &H || this == &D || this == &T)
		{
			if (k != MaxCharge_)
			{
				if (this == &D)
				{
					CX_[k].Stat_Tri(k, Tri_n_, &PRC96_D, n_D_1);
					if (K_DT)
					{
						CX_DT_[0].Stat_Tri(k, Tri_n_, &PRC96_D, n_T_1);
					}
				}
				else if (this == &T)
				{
					CX_[k].Stat_Tri(k, Tri_n_, &PRC96_D, n_T_1);
					if (K_DT)
					{
						CX_DT_[0].Stat_Tri(k, Tri_n_, &PRC96_D, n_D_1);
					}
				}
				if (K_database_Pra == 1)
				{
					Ion_[k].Stat_Tri(k, Tri_n_, &PLT12_H);
				}
				else
				{
					Ion_[k].StatEIRENE_Tri(k, Tri_n_, &R2_1_5_H10, 13.6);
				}
			}
			if (k != 0)
			{
				for (int i = 0; i < Grid4.num_tris(); i++)
				{
					if (this == &D)
					{
						Tri_n_[i][k] = n_D_1[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)];
					}
					else if (this == &D)
					{
						Tri_n_[i][k] = n_T_1[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)];
					}
				}
				Rec_[k].RecStat_Tri(k, Tri_n_, &PRB12_H);
			}
		}
		if (this == &H2 || this == &D2 || this == &T2)
		{
			if (k == 0)
			{
				Ion_[k].Stat_Tri(0, Tri_n_, NULL);
				CX_[k].Stat_Tri(0, Tri_n_, NULL);
				Ela_[k].Stat_Tri(0, Tri_n_, NULL);
				MAR_[k].Stat_Tri(0, Tri_n_, NULL);
				Diss1_[k].Stat_Tri(0, Tri_n_, NULL);
				Diss2_[k].Stat_Tri(0, Tri_n_, NULL);
			}
			else if (k == 1)
			{
				for (int i = 0; i < 3; i++)
				{
					DS_[k][i].Stat_Tri(0, Tri_n_, NULL);
				}
				for (int i = 0; i < Grid4.num_tris(); i++)
				{
					double Var_temp = R2_2_14_H2.cal(ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Te[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)]);
					DS_[1][2].SetPra_Tri(i, ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * Tri_n_[i][1] * R2_2_14_H8.cal(ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Te[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], &Var_temp));
				}
			}
		}
	}
}
int Particle::NowZone()
{
	double X_Temp[3] = {X_new_[0], X_new_[1], X_new_[2]};
	int Z1 = Zone_, Z2;
	for (int i = 0; i < 3; i++)
		X_new_[i] = X_[i];
	if (K_GRID)
		Grid1.Find(X_new_, &Zone_, XY_);
	else
		Find();
	Z2 = Zone_;
	Zone_ = Z1;
	for (int i = 0; i < 3; i++)
		X_new_[i] = X_Temp[i];
	return Z2;
}

void Particle::Clear(int n)
{
	if (n == 0)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
				for (int k = 0; k < MaxCharge_ + 1; k++)
				{
					n_[i][j][k] = 0.;
					T_[i][j][k] = 0.;
					E_[i][j][k] = 0.;
					Sum_n_[i][j][k] = 0.;
					Sum_E_[i][j][k] = 0.;
					for (int k = 0; k < 7; k++)
					{
						Sum_n_array_[i][j][0][k] = 0;
					}
					for (int m = 0; m < 3; m++)
						V_Grid_[i][j][m][0] = 0.;
					for (int m = 0; m < 3; m++)
						Sum_V_[i][j][m][0] = 0.;
				}
		Ion_[0].Clear();
		// std::cout << "Ion!" << endl;
		Rec_[1].Clear();
		// std::cout << "Rec!" << endl;
		CX_[0].Clear();
		// std::cout << "CX!" << endl;
		if (this == &H || this == &D || this == &T)
		{
			CX_DT_[0].Clear();
			R_with_H_[0].Clear();
			R_with_H2_[0].Clear();
		}
		if (this == &H2 || this == &D2 || this == &T2)
		{
			Ela_DT_[0].Clear();
			Ela_[0].Clear();
			// std::cout << "Ela!" << endl;
			Diss1_[0].Clear();
			// std::cout << "Diss1!" << endl;
			Diss2_[0].Clear();
			// std::cout << "Diss2!" << endl;
			Diss3_[0].Clear();
			// std::cout << "Diss3!" << endl;
			DS_[1][0].Clear();
			// std::cout << "DS1!" << endl;
			DS_[1][1].Clear();
			// std::cout << "DS2!" << endl;
			DS_[1][2].Clear();
			// std::cout << "DS3!" << endl;
		}
		R_with_H_[0].Clear();
		R_with_H2_[0].Clear();
	}
	if (K_Methane)
	{
		DS_[1][0].Clear();
		// std::cout << "DS1!" << endl;
		DS_[1][1].Clear();
		// std::cout << "DS2!" << endl;
		DS_[1][2].Clear();
		DS_[1][3].Clear();
	}
}

void Particle::Setdt(double dt1)
{
	dt_ = dt1;
}

double Particle::FlightTime()
{
	return dt_;
}

int Particle::Index(int i)
{
	return XY_[i];
}

double Particle::X(int i)
{
	return X_[i];
}

double Particle::X_new(int i)
{
	return X_new_[i];
}

double Particle::V(int i)
{
	return V_[i];
}

double Particle::Vreal(int i)
{
	if (Charge_ == 0)
		return V_[i];
	else
		return V_Charge_[i];
}

int Particle::XY(int i)
{
	return XY_[i];
}

int Particle::Tri_Index()
{
	return Tri_Index_;
}

/*int Particle::XY_new(int i)
{
	return XY_new_[i];
}*/

int Particle::boundary_start()
{
	return boundary_start_;
}

void Particle::SetXY(int i, int j)
{
	XY_[i] = j;
}

string Particle::name() { return name_; }

int Particle::fate(int i) { return fate_[i]; }

int Particle::sourcePar(int i) { return sourcePar_[i]; }

double Particle::CollProb(int i, int j, int k) { return CollProb_[i][j][k]; }

double Particle::CoreTotalcs(int i) { return CoreTotalcs_[i]; }

void Particle::calIonRateADAS(ADAS *ParColl)
{
	for (int i = 0; i < MaxCharge_; i++)
	{
		if (this == &H)
			Ion_[i].ADASInput(ParColl, i);
		if (this == &D)
			Ion_[i].ADASInput(ParColl, i, 1);
		if (this == &T)
			Ion_[i].ADASInput(ParColl, i, 1);
	}
}

void Particle::calRecRateADAS(ADAS *ParColl)
{
	for (int i = 1; i < MaxCharge_ + 1; i++)
	{
		if (this == &H)
			Rec_[i].ADASInput(ParColl, i - 1);
		if (this == &D)
			Rec_[i].ADASInput(ParColl, i - 1, 1);
		if (this == &T)
			Rec_[i].ADASInput(ParColl, i - 1, 1);
	}
}
void Particle::calCXRateADAS(ADAS *ParColl, int CX_cross)
{
	if (CX_cross == 0)
	{
		for (int i = 0; i < MaxCharge_; i++)
			CX_[i].ADASInput(ParColl, i);
	}
	else if (CX_cross == 1) // DT cross CX
		for (int i = 0; i < MaxCharge_; i++)
		{
			if (this == &D)
				CX_DT_[i].ADASInput(ParColl, i, 3, CX_cross);
			if (this == &T)
				CX_DT_[i].ADASInput(ParColl, i, 2, CX_cross);
		}
}

void Particle::AddWallEro(int num_Ero_wall)
{
	// std::cout << "all" << endl;
	AlltoWall_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	if (fate_[0] == 1) //"CX"
	{
		if (CXtoWall_.begin()->CollSource() == 0)
		{
			CXtoWall_.begin()->setCollSource(sourcePar_[0]);
			// std::cout << "CXtowall1" << ", " << sourcePar_[0] << endl;
			CXtoWall_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		else
		{
			for (int i = 0; i < CXtoWall_.size(); i++)
			{
				if (sourcePar_[0] == CXtoWall_[i].CollSource())
				{
					// std::cout << "CXtowall2" << ", " << sourcePar_[0] << endl;
					CXtoWall_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
					return;
				}
			}
			CXtoWall_.push_back(WallEro());
			CXtoWall_.back().WallEroInit(CXtoWall_.begin()->num_wall(), CXtoWall_.begin()->N_p(), CXtoWall_.begin()->N_r());
			CXtoWall_.back().setCollSource(sourcePar_[0]);
			// std::cout << "CXtowall3" << ", " << sourcePar_[0] << endl;
			CXtoWall_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
		}
	}
	if (fate_[0] == 4) //"rec"
	{
		if (RectoWall_.begin()->CollSource() == 0)
		{
			RectoWall_.begin()->setCollSource(sourcePar_[0]);
			// std::cout << "Rectowall1" << endl;
			RectoWall_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		else
		{
			for (int i = 0; i < RectoWall_.size(); i++)
			{
				if (sourcePar_[0] == RectoWall_[i].CollSource())
				{
					// std::cout << "Rectowall2" << endl;
					RectoWall_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
					return;
				}
			}
			RectoWall_.push_back(WallEro());
			RectoWall_.back().WallEroInit(RectoWall_.begin()->num_wall(), RectoWall_.begin()->N_p(), RectoWall_.begin()->N_r());
			RectoWall_.back().setCollSource(sourcePar_[0]);
			// std::cout << "Rectowall3" << endl;
			RectoWall_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
		}
	}
}

void Particle::OutWallEro(int fate)
{
	string filename = Outputpath + fatename(fate) + "WallEro.h5";
	std::cout << filename << endl;
	H5File *file = new H5File(filename, H5F_ACC_TRUNC);
	hsize_t dimsf1[3], *dimsf2; /* dataset dimensions */
	dimsf2 = (hsize_t *)malloc(sizeof(hsize_t));
	herr_t status;
	dimsf1[2] = N_radial;
	dimsf1[1] = N_poloidal;
	if (fate == 1)
	{
		dimsf1[0] = CXtoWall_.begin()->num_wall();
		dimsf2[0] = CXtoWall_.begin()->num_wall();
		std::vector<double> data1(dimsf1[0] * dimsf1[1] * dimsf1[2]);
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		std::cout << CXtoWall_.size() << endl;
		for (int i = 0; i < CXtoWall_.size(); i++)
		{
			string groupname = "/" + sourcename(CXtoWall_[i].CollSource());
			std::cout << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = CXtoWall_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1.data(), PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = CXtoWall_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1.data(), PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoWall_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoWall_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == 4)
	{
		dimsf1[0] = RectoWall_.begin()->num_wall();
		dimsf2[0] = RectoWall_.begin()->num_wall();
		std::vector<double> data1(dimsf1[0] * dimsf1[1] * dimsf1[2]);
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int i = 0; i < RectoWall_.size(); i++)
		{
			string groupname = "/" + sourcename(RectoWall_[i].CollSource());
			std::cout << RectoWall_[i].CollSource() << '\t' << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = RectoWall_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(
				file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1.data(), PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = RectoWall_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1.data(), PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoWall_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoWall_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == -1)
	{
		dimsf1[0] = AlltoWall_.begin()->num_wall();
		dimsf2[0] = AlltoWall_.begin()->num_wall();
		std::vector<double> data1(dimsf1[0] * dimsf1[1] * dimsf1[2]);
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = AlltoWall_.begin()->PartoWall[m][j][k];
		DataSet *dataset1 = new DataSet(
			file->createDataSet("n", PredType::NATIVE_DOUBLE, dataspace1));
		dataset1->write(data1.data(), PredType::NATIVE_DOUBLE);
		delete dataset1;

		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m * dimsf1[1] * dimsf1[2] + j * dimsf1[2] + k] = AlltoWall_.begin()->E_PartoWall[m][j][k];
		DataSet *dataset2 = new DataSet(
			file->createDataSet("T", PredType::NATIVE_DOUBLE, dataspace1));
		dataset2->write(data1.data(), PredType::NATIVE_DOUBLE);
		delete dataset2;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoWall_.begin()->All_PartoWall[m];
		DataSet *dataset3 = new DataSet(
			file->createDataSet("n_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset3->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset3;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoWall_.begin()->All_E_PartoWall[m];
		DataSet *dataset4 = new DataSet(
			file->createDataSet("T_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset4->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset4;
	}
	delete file;
}

void Particle::AddTargetEro(int num_Ero_wall)
{
	// std::cout << "AlltoTarget1" << ", " << sourcePar_[0] << endl;
	AlltoTarget_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	if (fate_[0] == 1) //"CX"
	{
		if (CXtoTarget_.begin()->CollSource() == 0)
		{
			CXtoTarget_.begin()->setCollSource(sourcePar_[0]);

			// std::cout << "CXtoTarget1" << ", " << sourcePar_[0] << endl;
			CXtoTarget_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		for (int i = 0; i < CXtoTarget_.size(); i++)
		{
			if (sourcePar_[0] == CXtoTarget_[i].CollSource())
			{
				// std::cout << "CXtoTarget2" << ", " << sourcePar_[0] << endl;
				CXtoTarget_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
				return;
			}
		}
		CXtoTarget_.push_back(WallEro());
		CXtoTarget_.back().WallEroInit(CXtoTarget_.begin()->num_wall(), CXtoTarget_.begin()->N_p(), CXtoTarget_.begin()->N_r());
		CXtoTarget_.back().setCollSource(sourcePar_[0]);
		// std::cout << "CXtoTarget3" << ", " << sourcePar_[0] << endl;
		CXtoTarget_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	}
	if (fate_[0] == 4) //"rec"
	{
		if (RectoTarget_.begin()->CollSource() == 0)
		{
			RectoTarget_.begin()->setCollSource(sourcePar_[0]);
			// std::cout << "RectoTarget1" << ", " << sourcePar_[0] << endl;
			RectoTarget_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		for (int i = 0; i < RectoTarget_.size(); i++)
		{
			if (sourcePar_[0] == RectoTarget_[i].CollSource())
			{
				// std::cout << "RectoTarget2" << ", " << sourcePar_[0] << endl;
				RectoTarget_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
				return;
			}
		}
		RectoTarget_.push_back(WallEro());
		RectoTarget_.back().WallEroInit(RectoTarget_.begin()->num_wall(), RectoTarget_.begin()->N_p(), RectoTarget_.begin()->N_r());
		RectoTarget_.back().setCollSource(sourcePar_[0]);
		// std::cout << "RectoTarget3" << ", " << sourcePar_[0] << endl;
		RectoTarget_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	}
}

void Particle::OutTargetEro(int fate)
{
	string filename = Outputpath + fatename(fate) + "TargetEro.h5";
	std::cout << filename << endl;
	H5File *file = new H5File(filename, H5F_ACC_TRUNC);
	hsize_t dimsf1[3], *dimsf2; /* dataset dimensions */
	dimsf2 = (hsize_t *)malloc(sizeof(hsize_t));
	herr_t status;
	dimsf1[2] = N_radial;
	dimsf1[1] = N_poloidal;
	if (fate == 1)
	{
		dimsf1[0] = CXtoTarget_.begin()->num_wall();
		dimsf2[0] = CXtoTarget_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		std::cout << CXtoTarget_.size() << endl;
		for (int i = 0; i < CXtoTarget_.size(); i++)
		{
			string groupname = "/" + sourcename(CXtoTarget_[i].CollSource());
			std::cout << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = CXtoTarget_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(
				file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = CXtoTarget_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoTarget_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoTarget_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == 4)
	{
		dimsf1[0] = RectoTarget_.begin()->num_wall();
		dimsf2[0] = RectoTarget_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int i = 0; i < RectoTarget_.size(); i++)
		{
			string groupname = "/" + sourcename(RectoTarget_[i].CollSource());
			std::cout << RectoTarget_[i].CollSource() << '\t' << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = RectoTarget_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(
				file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = RectoTarget_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoTarget_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoTarget_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == -1)
	{
		dimsf1[0] = AlltoTarget_.begin()->num_wall();
		dimsf2[0] = AlltoTarget_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m][j][k] = AlltoTarget_.begin()->PartoWall[m][j][k];
		DataSet *dataset1 = new DataSet(
			file->createDataSet("n", PredType::NATIVE_DOUBLE, dataspace1));
		dataset1->write(data1, PredType::NATIVE_DOUBLE);
		delete dataset1;

		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m][j][k] = AlltoTarget_.begin()->E_PartoWall[m][j][k];
		DataSet *dataset2 = new DataSet(
			file->createDataSet("T", PredType::NATIVE_DOUBLE, dataspace1));
		dataset2->write(data1, PredType::NATIVE_DOUBLE);
		delete dataset2;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoTarget_.begin()->All_PartoWall[m];
		DataSet *dataset3 = new DataSet(
			file->createDataSet("n_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset3->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset3;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoTarget_.begin()->All_E_PartoWall[m];
		DataSet *dataset4 = new DataSet(
			file->createDataSet("T_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset4->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset4;
	}
	delete file;
}

void Particle::AddPlasmaBoundaryEro(int num_Ero_wall)
{
	AlltoPlasmaBoundary_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	if (fate_[0] == 1) //"CX"
	{
		if (CXtoPlasmaBoundary_.begin()->CollSource() == 0)
		{
			CXtoPlasmaBoundary_.begin()->setCollSource(sourcePar_[0]);
			CXtoPlasmaBoundary_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		for (int i = 0; i < CXtoPlasmaBoundary_.size(); i++)
		{
			if (sourcePar_[0] == CXtoPlasmaBoundary_[i].CollSource())
			{
				CXtoPlasmaBoundary_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
				return;
			}
		}
		CXtoPlasmaBoundary_.push_back(WallEro());
		CXtoPlasmaBoundary_.back().WallEroInit(CXtoPlasmaBoundary_.begin()->num_wall(), CXtoPlasmaBoundary_.begin()->N_p(), CXtoPlasmaBoundary_.begin()->N_r());
		CXtoPlasmaBoundary_.back().setCollSource(sourcePar_[0]);
		CXtoPlasmaBoundary_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	}
	if (fate_[0] == 4) //"rec"
	{
		if (RectoPlasmaBoundary_.begin()->CollSource() == 0)
		{
			RectoPlasmaBoundary_.begin()->setCollSource(sourcePar_[0]);
			RectoPlasmaBoundary_.begin()->AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
			return;
		}
		for (int i = 0; i < RectoPlasmaBoundary_.size(); i++)
		{
			if (sourcePar_[0] == RectoPlasmaBoundary_[i].CollSource())
			{
				RectoPlasmaBoundary_[i].AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
				return;
			}
		}
		RectoPlasmaBoundary_.push_back(WallEro());
		RectoPlasmaBoundary_.back().WallEroInit(RectoPlasmaBoundary_.begin()->num_wall(), RectoPlasmaBoundary_.begin()->N_p(), RectoPlasmaBoundary_.begin()->N_r());
		RectoPlasmaBoundary_.back().setCollSource(sourcePar_[0]);
		RectoPlasmaBoundary_.back().AddPar(num_Ero_wall, Weight_ * NumPar_now, Tn_, fate_[0], sourceGrid_);
	}
}

void Particle::OutPlasmaBoundaryEro(int fate)
{
	string filename = Outputpath + fatename(fate) + "PlasmaBoundaryEro.h5";
	std::cout << filename << endl;
	H5File *file = new H5File(filename, H5F_ACC_TRUNC);
	hsize_t dimsf1[3], *dimsf2; /* dataset dimensions */
	dimsf2 = (hsize_t *)malloc(sizeof(hsize_t));
	herr_t status;
	dimsf1[2] = N_radial;
	dimsf1[1] = N_poloidal;
	if (fate == 1)
	{
		dimsf1[0] = CXtoPlasmaBoundary_.begin()->num_wall();
		dimsf2[0] = CXtoPlasmaBoundary_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		std::cout << CXtoPlasmaBoundary_.size() << endl;
		for (int i = 0; i < CXtoPlasmaBoundary_.size(); i++)
		{
			string groupname = "/" + sourcename(CXtoPlasmaBoundary_[i].CollSource());
			std::cout << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = CXtoPlasmaBoundary_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(
				file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = CXtoPlasmaBoundary_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoPlasmaBoundary_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = CXtoPlasmaBoundary_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == 4)
	{
		dimsf1[0] = RectoPlasmaBoundary_.begin()->num_wall();
		dimsf2[0] = RectoPlasmaBoundary_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int i = 0; i < RectoPlasmaBoundary_.size(); i++)
		{
			string groupname = "/" + sourcename(RectoPlasmaBoundary_[i].CollSource());
			std::cout << RectoPlasmaBoundary_[i].CollSource() << '\t' << groupname << endl;
			Group *group = new Group(file->createGroup(groupname));

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = RectoPlasmaBoundary_[i].PartoWall[m][j][k];
			DataSet *dataset1 = new DataSet(
				file->createDataSet(groupname + "/n", PredType::NATIVE_DOUBLE, dataspace1));
			dataset1->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset1;

			for (int m = 0; m < dimsf1[0]; m++)
				for (int j = 0; j < dimsf1[1]; j++)
					for (int k = 0; k < dimsf1[2]; k++)
						data1[m][j][k] = RectoPlasmaBoundary_[i].E_PartoWall[m][j][k];
			DataSet *dataset2 = new DataSet(
				file->createDataSet(groupname + "/T", PredType::NATIVE_DOUBLE, dataspace1));
			dataset2->write(data1, PredType::NATIVE_DOUBLE);
			delete dataset2;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoPlasmaBoundary_[i].All_PartoWall[m];
			DataSet *dataset3 = new DataSet(
				file->createDataSet(groupname + "/n_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset3->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset3;

			for (int m = 0; m < dimsf1[0]; m++)
				data2[m] = RectoPlasmaBoundary_[i].All_E_PartoWall[m];
			DataSet *dataset4 = new DataSet(
				file->createDataSet(groupname + "/T_all", PredType::NATIVE_DOUBLE, dataspace2));
			dataset4->write(data2, PredType::NATIVE_DOUBLE);
			delete dataset4;
			delete group;
		}
	}
	else if (fate == -1)
	{
		dimsf1[0] = AlltoPlasmaBoundary_.begin()->num_wall();
		dimsf2[0] = AlltoPlasmaBoundary_.begin()->num_wall();
		double data1[dimsf1[0]][dimsf1[1]][dimsf1[2]];
		double data2[dimsf2[0]];
		DataSpace dataspace1(3, dimsf1);
		DataSpace dataspace2(1, dimsf2);
		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m][j][k] = AlltoPlasmaBoundary_.begin()->PartoWall[m][j][k];
		DataSet *dataset1 = new DataSet(
			file->createDataSet("n", PredType::NATIVE_DOUBLE, dataspace1));
		dataset1->write(data1, PredType::NATIVE_DOUBLE);
		delete dataset1;

		for (int m = 0; m < dimsf1[0]; m++)
			for (int j = 0; j < dimsf1[1]; j++)
				for (int k = 0; k < dimsf1[2]; k++)
					data1[m][j][k] = AlltoPlasmaBoundary_.begin()->E_PartoWall[m][j][k];
		DataSet *dataset2 = new DataSet(
			file->createDataSet("T", PredType::NATIVE_DOUBLE, dataspace1));
		dataset2->write(data1, PredType::NATIVE_DOUBLE);
		delete dataset2;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoPlasmaBoundary_.begin()->All_PartoWall[m];
		DataSet *dataset3 = new DataSet(
			file->createDataSet("n_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset3->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset3;

		for (int m = 0; m < dimsf1[0]; m++)
			data2[m] = AlltoPlasmaBoundary_.begin()->All_E_PartoWall[m];
		DataSet *dataset4 = new DataSet(
			file->createDataSet("T_all", PredType::NATIVE_DOUBLE, dataspace2));
		dataset4->write(data2, PredType::NATIVE_DOUBLE);
		delete dataset4;
	}
	delete file;
}

double Particle::CalAngle(int num_wall)
{
	// std::cout << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
	// std::cout << num_wall << '\t' << cosang[num_wall] << '\t' << sinang[num_wall] << endl;
	// std::cout << 180 / pi * acos(abs(V_[0] * sinang[num_wall] - V_[1] * cosang[num_wall]) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)))) << endl;
	if (MeshMode == 1)
	{
		if (InterscePoint[0][4] == 11)
			return 180 / pi * acos(abs(V_[0] * sinang[num_wall] - V_[1] * cosang[num_wall]) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2))));
		else if (InterscePoint[0][4] == 1)
			return 180 / pi * acos(abs(V_[0] * Grid1.Sin_target(num_wall) - V_[1] * Grid1.Cos_target(num_wall)) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2))));
		else
			throw std::logic_error("angle calculation of reflect have some problem in Calangle()");
	}
	else if (MeshMode == 3)
	{
		if (InterscePoint[0][4] == 11)
			return 180 / pi * acos(abs(V_[0] * Grid4.Sin_trimesh(Tri_Index_, num_wall) - V_[1] * Grid4.Cos_trimesh(Tri_Index_, num_wall)) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2))));
		else if (InterscePoint[0][4] == 1)
			return 180 / pi * acos(abs(V_[0] * Grid4.Sin_trimesh(Tri_Index_, num_wall) - V_[1] * Grid4.Cos_trimesh(Tri_Index_, num_wall)) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2))));
		else
			throw std::logic_error("angle calculation of reflect have some problem in Calangle()");
	}
	else
		throw std::logic_error("angle calculation of reflect have some problem in Calangle()");
}

double Particle::CalAngle(double Sin, double Cos)
{
	// std::cout << V_[0] << '\t' << V_[1] << '\t' << V_[2] << endl;
	// std::cout << num_wall << '\t' << cosang[num_wall] << '\t' << sinang[num_wall] << endl;
	// std::cout << 180 / pi * acos(abs(V_[0] * sinang[num_wall] - V_[1] * cosang[num_wall]) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)))) << endl;
	return 180 / pi * acos(abs(V_[0] * Sin - V_[1] * Cos) / sqrt((pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2))));
}

int Particle::MaxCharge() { return MaxCharge_; }

void Particle::SetCharge(int i) { Charge_ = i; }

double Particle::mass() { return mass_; }

double Particle::V_Charge(int i) { return V_Charge_[i]; }

double Particle::V_Grid(int i, int j, int k, int m)
{
	return V_Grid_[i][j][k][m];
}
double Particle::V_D_1(int i, int j, int k, int m)
{
	return V_D_1_[i][j][k][m];
}
double Particle::V_Grid_CX_Ion_Af(int i, int j, int k, int m)
{
	return V_Grid_CX_Ion_Af_[i][j][k][m];
}
double Particle::V_Grid_CX_Ion_Be(int i, int j, int k, int m)
{
	return V_Grid_CX_Ion_Be_[i][j][k][m];
}
double Particle::V_Grid_CX_Neu_Af(int i, int j, int k, int m)
{
	return V_Grid_CX_Neu_Af_[i][j][k][m];
}
double Particle::V_Grid_CX_Neu_Be(int i, int j, int k, int m)
{
	return V_Grid_CX_Neu_Be_[i][j][k][m];
}

int Particle::Zone() { return Zone_; }

/*int Particle::Zone_new()
{
	return Zone_new_;
}*/

int Particle::Charge() { return Charge_; }

double Particle::Tn() { return Tn_; }

int Particle::IfColl()
{
	return IfColl_;
}

int Particle::IfFlightOut()
{
	return IfFlightOut_;
}

void Particle::SetIfFlightOut(int i)
{
	IfFlightOut_ = i;
}

double Particle::n(int i, int j, int k) { return n_[i][j][k]; }

double Particle::T_n(int i, int j, int k)
{
	return T_[i][j][k];
}

double Particle::Gamma()
{
	return Weight_ * NumPar_now;
}

double Particle::Weight_Target(int i) { return Weight_Target_[i]; }

double Particle::NumPar_Target(int i) { return NumPar_Target_[i]; }

double Particle::Weight() { return Weight_; }

int Particle::ChargeTag() { return ChargeTag_; }

void Particle::SetWeight(double w) { Weight_ = w; }

void Particle::SetsourceWall(double *sourceWall)
{
	for (int i = 0; i < 4; i++)
		sourceWall_[i] = sourceWall[i];
}

void Particle::SetTn(double T) { Tn_ = T; }

void Particle::SetZone(int zone)
{
	Zone_ = zone;
}

void Particle::SetX(int i, double X) { X_[i] = X; }

void Particle::SetV(int i, double V) { V_[i] = V; }

void Particle::SetX_new(int i, double X) { X_new_[i] = X; }

void Particle::setname(string particle_name) { name_ = particle_name; }

void Particle::SetChargeTag(int i) { ChargeTag_ = i; }

void Particle::setPar(string particle_name, double mass, int Charge)
{
	name_ = particle_name;
	mass_ = mass;
	Charge_ = Charge;
	Zone_ = 0;
}

void Particle::Dump_Flux()
{
	ofstream out;
	string pathout;

	/*pathout = Outputpath + name_ + "_" + "Flux_Grid.txt";
	out.open(pathout);
	for (int i = 0; i < num_GridBoundry; i++)
		out << Flux_GridBoundry[i][0] << '\t' << Flux_GridBoundry[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "Flux_PFR.txt";
	out.open(pathout);
	for (int i = 0; i < num_PFRBoundry; i++)
		out << Flux_PFRBoundry[i][0] << '\t' << Flux_PFRBoundry[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "Flux_Core.txt";
	out.open(pathout);
	for (int i = 0; i < num_CoreBoundry; i++)
		out << Flux_CoreBoundry[i][0] << '\t' << Flux_CoreBoundry[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "Flux_Target_l.txt";
	out.open(pathout);
	for (int i = 0; i < 36; i++)
		out << Flux_Target_l[i][0] << '\t' << Flux_Target_l[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "Flux_Target_r.txt";
	out.open(pathout);
	for (int i = 0; i < 36; i++)
		out << Flux_Target_r[i][0] << '\t' << Flux_Target_r[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "EneFlux_Grid.txt";
	out.open(pathout);
	for (int i = 0; i < num_GridBoundry; i++)
		out << EneFlux_GridBoundry[i][0] << '\t' << EneFlux_GridBoundry[i][1]
			<< endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "EneFlux_PFR.txt";
	out.open(pathout);
	for (int i = 0; i < num_PFRBoundry; i++)
		out << EneFlux_PFRBoundry[i][0] << '\t' << EneFlux_PFRBoundry[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "EneFlux_Core.txt";
	out.open(pathout);
	for (int i = 0; i < num_CoreBoundry; i++)
		out << EneFlux_CoreBoundry[i][0] << '\t' << EneFlux_CoreBoundry[i][1]
			<< endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "EneFlux_Target_l.txt";
	out.open(pathout);
	for (int i = 0; i < 36; i++)
		out << EneFlux_Target_l[i][0] << '\t' << EneFlux_Target_l[i][1] << endl;
	out.close();

	pathout = Outputpath + name_ + "_" + "EneFlux_Target_r.txt";
	out.open(pathout);
	for (int i = 0; i < 36; i++)
		out << EneFlux_Target_r[i][0] << '\t' << EneFlux_Target_r[i][1] << endl;
	out.close();*/
	FluxOutput();
	FT_.write_H5(Outputpath + name_ + "_FluxZone.h5");
}

double Particle::sourceWall(int i) { return sourceWall_[i]; }

void Particle::divimp_Ti()
{
	if (Charge_ == 0)
	{
		return;
	}
	double b;
	static double u, u5_2, u3_2, u2;
	u = mass_ / (mass_ + Dmass);
	/*if (K_back == 1)
		u = mass_ / (mass_ + Hmass);
	if (K_back == 2)
		u = mass_ / (mass_ + Dmass);
	if (K_back == 3)
		u = mass_ / (mass_ + Tmass);*/
	u2 = pow(u, 2.0);
	u5_2 = pow(u, 2.5);
	u3_2 = pow(u, 1.5);

	int i, j, q;
	double dv;
	i = XY_[0];
	j = XY_[1];
	q = Charge_;

	b = -3.0 * (1.0 - u - 5.0 * 1.414 * q * q * (1.1 * u5_2 - 0.35 * u3_2)) /
		(2.6 - 2.0 * u + 5.4 * u2);
	dv = b * dTip[i][j] * qe / mass_ * dt * fabs(Bpol[i][j] / Btor[i][j]);
	Vlog << dv << endl;
	V_Charge_[3] += dv;
}

void Particle::divimp_F()
{
	if (Charge_ == 0)
	{
		return;
	}
	int i, j, q;
	i = XY_[0];
	j = XY_[1];
	q = Charge_;
	double ts1, lnA;
	lnA = 15.0;
	ts1 = mass_ * 6.02e26 * Ti[i][j] * sqrt(Ti[i][j] / (Dmass * 6.02e26)) /
		  (6.8 * 1e4 * (1.0 + Dmass / mass_) * (n_D_1[i][j] + n_T_1[i][j]) * 1e-18 * lnA);
	/*if (K_back == 1)
		ts1 = mass_ * 6.02e26 * Ti[i][j] * sqrt(Ti[i][j] / (Hmass * 6.02e26)) /
			  (6.8 * 1e4 * (1.0 + Hmass / mass_) * ni[i][j] * 1e-18 * lnA);
	if (K_back == 2)
		ts1 = mass_ * 6.02e26 * Ti[i][j] * sqrt(Ti[i][j] / (Dmass * 6.02e26)) /
			  (6.8 * 1e4 * (1.0 + Dmass / mass_) * ni[i][j] * 1e-18 * lnA);
	if (K_back == 3)
		ts1 = mass_ * 6.02e26 * Ti[i][j] * sqrt(Ti[i][j] / (Tmass * 6.02e26)) /
			  (6.8 * 1e4 * (1.0 + Tmass / mass_) * ni[i][j] * 1e-18 * lnA);*/
	double ts2;
	ts2 = ts1 / q / q;
	double dv;
	// dv = (vp[i][j] - p->vx) / ts2;
	dv = (ua_D_1[i][j] * 1 - V_Charge_[3]) / ts2 * dt;
	Vlog << dv << ' ';
	V_Charge_[3] += dv;
}

void Particle::divimp_E()
{
	int i, j, q;
	i = XY_[0];
	j = XY_[1];
	q = Charge_;
	V_Charge_[3] +=
		Charge_ * qe * Epol[i][j] / mass_ * dt * Bpol[i][j] / Btor[i][j];
}

void Particle::divimp_drift_E()
{
	int i, j;
	double vdp;
	double vdr;
	if (Charge_ == 0)
	{
		return;
	}
	i = XY_[0];
	j = XY_[1];

	// 电势问题
	// 取-1.0
	vdp = vdEp[i][j];
	vdr = vdEr[i][j];

	/*pull_mode*/
	/*
	0:epx[i][j]
	1:p->epx
	*/

	X_[0] += vdp * dt * epx[i][j];
	X_[1] += vdp * dt * epy[i][j];
	X_[0] += vdr * dt * erx[i][j];
	X_[1] += vdr * dt * ery[i][j];
}

void Particle::divimp_anomalous_diffusion()
{
	if (Charge_ == 0)
	{
		return;
	}
	int i, j;
	int Nt;
	double d;
	double dx, dy;
	double dx31, dx42, dy31, dy42;
	double ex, ey;
	double Dan;
	if (K_dn == 2)
	{
		if (XY_[0] < 25 || XY_[0] > 72)
			Dan = 0.3;
		else
			Dan = Dn[XY_[1]];
	}
	else if (K_dn == 1)
		Dan = 0.3;
	double rand_g;
	i = XY_[0];
	j = XY_[1];
	Nt = N_poloidal;
	if (i == 1 || i == Nt - 2)
	{
		return;
	}

	d = sqrt(2.0 * Dan * dt);
	rand_g = Random();
	if (rand_g > 0.5)
	{
		rand_g = 1.0;
	}
	else
	{
		rand_g = -1.0;
	}
	// rand_g = 1.0;
	// if (i == 0 || i == 1 || i == Nt - 2 || i == Nt - 1)
	// if (i == 0 || i == Nt - 1)
	//{
	//	dx31 = crx3[i][j] - crx1[i][j];
	//	dy31 = cry3[i][j] - cry1[i][j];
	//	dx42 = crx4[i][j] - crx2[i][j];
	//	dy42 = cry4[i][j] - cry2[i][j];
	//	dx = (dx31 + dx42) / 2.0;
	//	dy = (dy31 + dy42) / 2.0;
	//	ex = dx / sqrt(dx*dx + dy*dy);
	//	ey = dy / sqrt(dx*dx + dy*dy);
	//	dx = d*rand_g*ex;
	//	dy = d*rand_g*ey;
	//	p->x = p->x + dx;
	//	p->y = p->y + dy;
	// }
	// else
	//{
	//	dx = d*rand_g*erx[i][j];
	//	dy = d*rand_g*ery[i][j];
	//	p->x = p->x + dx;
	//	p->y = p->y + dy;
	// }
	dx = d * rand_g * erx[i][j];
	dy = d * rand_g * ery[i][j];
	X_new_[0] += dx;
	X_new_[1] += dy;
}

void Particle::Dump_array_2D(string type, string coll, int Charge)
{
	ofstream out;
	string pathout;
	for (int k = 0; k < Num_array_ + 1; k++)
	{
		if (coll != "n" && coll != "T")
			pathout = Outputpath + type + "_" + coll + "_" + name_ + "_" + to_string(Charge) + "_" + to_string(k);
		else
			pathout = Outputpath + type + "_" + name_ + "_" + to_string(Charge) + "_" + to_string(k);

		out.open(pathout);
		if (type == "n")
		{
			for (int i = 0; i < N_poloidal; i++)
			{
				for (int j = 0; j < N_radial; j++)
					out << Sum_n_array_[i][j][Charge][k] << "\t";
				out << endl;
			}
			out.close();
		}
	}
}

void Particle::Dump_2D(string type, string coll, int Charge)
{
	ofstream out;
	string pathout;
	if (coll != "n" && coll != "T")
		pathout = Outputpath + type + "_" + coll + "_" + name_ + "_" + to_string(Charge);
	else
		pathout = Outputpath + type + "_" + name_ + "_" + to_string(Charge);
	out.open(pathout);

	if (type == "n")
	{
		for (int i = 0; i < N_poloidal; i++)
		{
			for (int j = 0; j < N_radial; j++)
				out << n_[i][j][Charge] << "\t";
			out << endl;
		}
		out.close();
		return;
	}
	if (type == "T")
	{
		for (int i = 0; i < N_poloidal; i++)
		{
			for (int j = 0; j < N_radial; j++)
				out << T_[i][j][Charge] << "\t";
			out << endl;
		}
		out.close();
		return;
	}
	if (type == "E")
	{
		for (int i = 0; i < N_poloidal; i++)
		{
			for (int j = 0; j < N_radial; j++)
				out << E_[i][j][Charge] << "\t";
			out << endl;
		}
		out.close();
		return;
	}
	if (type == "Sn")
	{
		if (coll == "Ion")
		{
			Ion_[Charge].OutSn(out);
			return;
		}
		if (coll == "Rec")
		{
			Rec_[Charge].OutSn(out);
			return;
		}
		if (coll == "CX")
		{
			CX_[Charge].OutSn(out);
			return;
		}
		if (coll == "CXDT")
		{
			CX_DT_[Charge].OutSn(out);
			return;
		}
		if (coll == "Ela")
		{
			Ela_[Charge].OutSn(out);
			return;
		}
		if (coll == "Mar")
		{
			MAR_[Charge].OutSn(out);
			return;
		}
		if (coll == "Diss1")
		{
			Diss1_[Charge].OutSn(out);
			return;
		}
		if (coll == "Diss2")
		{
			Diss2_[Charge].OutSn(out);
			return;
		}
		if (coll == "Diss3")
		{
			Diss3_[Charge].OutSn(out);
			return;
		}
		if (coll == "DS1")
		{
			DS_[Charge][0].OutSn(out);
			return;
		}
		if (coll == "DS2")
		{
			DS_[Charge][1].OutSn(out);
			return;
		}
		if (coll == "DS3")
		{
			DS_[Charge][2].OutSn(out);
			return;
		}
		if (coll == "DS4")
		{
			DS_[Charge][3].OutSn(out);
			return;
		}
		out.close();
		return;
	}
	if (type == "Smu")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSmu(out);
		if (coll == "Rec")
			Rec_[Charge].OutSmu(out);
		if (coll == "CX")
			CX_[Charge].OutSmu(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSmu(out);
		if (coll == "Ela")
			Ela_[Charge].OutSmu(out);
	}
	if (type == "Smu1")
	{
		if (coll == "CX")
			CX_[Charge].OutSmu1(out);
	}
	if (type == "Smu2")
	{
		if (coll == "CX")
			CX_[Charge].OutSmu2(out);
	}
	if (type == "SE")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSE(out);
		if (coll == "Rec")
			Rec_[Charge].OutSE(out);
		if (coll == "CX")
			CX_[Charge].OutSE(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSE(out);
		if (coll == "Ela")
			Ela_[Charge].OutSE(out);
	}
	if (type == "SE_n")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSE_n(out);
		if (coll == "Rec")
			Rec_[Charge].OutSE_n(out);
		if (coll == "CX")
			CX_[Charge].OutSE_n(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSE_n(out);
		if (coll == "Ela")
			Ela_[Charge].OutSE_n(out);
	}
	if (type == "Pra")
	{
		if (coll == "Ion")
			Ion_[Charge].OutPra(out);
		if (coll == "Rec")
			Rec_[Charge].OutPra(out);
		if (coll == "CX")
			CX_[Charge].OutPra(out);
		if (coll == "DS3")
			DS_[Charge][2].OutPra(out);
	}
	/*for (int i = 0; i < N_poloidal; i++)
	{
			for (int j = 0; j < N_radial; j++)
					out << VarOut[i][j] << " ";
			out << endl;
	}*/
	out.close();
}

void Particle::Dump_Tri(string type, string coll, int Charge)
{
	ofstream out;
	string pathout;
	if (coll != "n" && coll != "T")
		pathout = Outputpath + type + "_" + coll + "_" + name_ + "_" + to_string(Charge) + "_Tri";
	else
		pathout = Outputpath + type + "_" + name_ + "_" + to_string(Charge) + "_Tri";
	out.open(pathout);

	if (type == "n")
	{
		for (int i = 0; i < Grid4.num_tris(); i++)
		{
			out << Tri_n_[i][Charge] << endl;
		}
		out.close();
		return;
	}
	if (type == "T")
	{
		for (int i = 0; i < Grid4.num_tris(); i++)
		{
			out << Tri_T_[i][Charge] << endl;
		}
		out.close();
		return;
	}
	if (type == "E")
	{
		for (int i = 0; i < Grid4.num_tris(); i++)
		{
			out << Tri_E_[i][Charge] << endl;
		}
		out.close();
		return;
	}
	if (type == "Sn")
	{
		if (coll == "Ion")
		{
			Ion_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Rec")
		{
			Rec_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "CX")
		{
			CX_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "CXDT")
		{
			CX_DT_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Ela")
		{
			Ela_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Mar")
		{
			MAR_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Diss1")
		{
			Diss1_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Diss2")
		{
			Diss2_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "Diss3")
		{
			Diss3_[Charge].OutSn_Tri(out);
			return;
		}
		if (coll == "DS1")
		{
			DS_[Charge][0].OutSn_Tri(out);
			return;
		}
		if (coll == "DS2")
		{
			DS_[Charge][1].OutSn_Tri(out);
			return;
		}
		if (coll == "DS3")
		{
			DS_[Charge][2].OutSn_Tri(out);
			return;
		}
		if (coll == "DS4")
		{
			DS_[Charge][3].OutSn_Tri(out);
			return;
		}
		out.close();
		return;
	}
	if (type == "Smu")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSmu_Tri(out);
		if (coll == "Rec")
			Rec_[Charge].OutSmu_Tri(out);
		if (coll == "CX")
			CX_[Charge].OutSmu_Tri(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSmu_Tri(out);
		if (coll == "Ela")
			Ela_[Charge].OutSmu_Tri(out);
	}
	if (type == "Smu1")
	{
		if (coll == "CX")
			CX_[Charge].OutSmu1(out);
	}
	if (type == "Smu2")
	{
		if (coll == "CX")
			CX_[Charge].OutSmu2(out);
	}
	if (type == "SE")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSE_Tri(out);
		if (coll == "Rec")
			Rec_[Charge].OutSE_Tri(out);
		if (coll == "CX")
			CX_[Charge].OutSE_Tri(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSE_Tri(out);
		if (coll == "Ela")
			Ela_[Charge].OutSE_Tri(out);
	}
	if (type == "SE_n")
	{
		if (coll == "Ion")
			Ion_[Charge].OutSE_n(out);
		if (coll == "Rec")
			Rec_[Charge].OutSE_n(out);
		if (coll == "CX")
			CX_[Charge].OutSE_n(out);
		if (coll == "CXDT")
			CX_DT_[Charge].OutSE_n(out);
		if (coll == "Ela")
			Ela_[Charge].OutSE_n(out);
	}
	if (type == "Pra")
	{
		if (coll == "Ion")
			Ion_[Charge].OutPra_Tri(out);
		if (coll == "Rec")
			Rec_[Charge].OutPra_Tri(out);
		if (coll == "CX")
			CX_[Charge].OutPra_Tri(out);
		if (coll == "DS3")
			DS_[Charge][2].OutPra_Tri(out);
	}
	/*for (int i = 0; i < N_poloidal; i++)
	{
			for (int j = 0; j < N_radial; j++)
					out << VarOut[i][j] << " ";
			out << endl;
	}*/
	out.close();
}

/*void Particle::FluxCal_Grid()
{
	double i_x, i_y, res;
	for (int i = 0; i < num_CoreBoundry; i++)
		if (get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1],
								  CoreBoundry[i][0], CoreBoundry[i][1],
								  CoreBoundry[i + 1][0], CoreBoundry[i + 1][1],
								  &i_x, &i_y))
		{
			res = (CoreBoundry[i + 1][0] - CoreBoundry[i][0]) *
					  (X_[1] - CoreBoundry[i][1]) -
				  (CoreBoundry[i + 1][1] - CoreBoundry[i][1]) *
					  (X_[0] - CoreBoundry[i][0]);
			if (res < 0)
			{
				Flux_CoreBoundry[i][0] += Weight_ * NumPar_now;
				EneFlux_CoreBoundry[i][0] += Weight_ * NumPar_now * Tn_;
			}
			else if (res > 0)
			{
				Flux_CoreBoundry[i][1] += Weight_ * NumPar_now;
				EneFlux_CoreBoundry[i][1] += Weight_ * NumPar_now * Tn_;
			}
			else
			{
				res = (CoreBoundry[i + 1][0] - CoreBoundry[i][0]) *
						  (X_new_[1] - CoreBoundry[i][1]) -
					  (CoreBoundry[i + 1][1] - CoreBoundry[i][1]) *
						  (X_new_[0] - CoreBoundry[i][0]);
				if (res > 0)
				{
					Flux_CoreBoundry[i][0] += Weight_ * NumPar_now;
					EneFlux_CoreBoundry[i][0] += Weight_ * NumPar_now * Tn_;
				}
				else if (res < 0)
				{
					Flux_CoreBoundry[i][1] += Weight_ * NumPar_now;
					EneFlux_CoreBoundry[i][1] += Weight_ * NumPar_now * Tn_;
				}
			}
		}

	for (int i = 0; i < num_GridBoundry; i++)
		if (get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1],
								  GridBoundry[i][0], GridBoundry[i][1],
								  GridBoundry[i + 1][0], GridBoundry[i + 1][1],
								  &i_x, &i_y))
		{
			res = (GridBoundry[i + 1][0] - GridBoundry[i][0]) *
					  (X_[1] - GridBoundry[i][1]) -
				  (GridBoundry[i + 1][1] - GridBoundry[i][1]) *
					  (X_[0] - GridBoundry[i][0]);
			if (res < 0)
			{
				Flux_GridBoundry[i][1] += Weight_ * NumPar_now;
				EneFlux_GridBoundry[i][1] += Weight_ * NumPar_now * Tn_;
			}
			else if (res > 0)
			{
				Flux_GridBoundry[i][0] += Weight_ * NumPar_now;
				EneFlux_GridBoundry[i][0] += Weight_ * NumPar_now * Tn_;
			}
			else
			{
				res = (GridBoundry[i + 1][0] - GridBoundry[i][0]) *
						  (X_new_[1] - GridBoundry[i][1]) -
					  (GridBoundry[i + 1][1] - GridBoundry[i][1]) *
						  (X_new_[0] - GridBoundry[i][0]);
				if (res > 0)
				{
					Flux_GridBoundry[i][1] += Weight_ * NumPar_now;
					EneFlux_GridBoundry[i][1] += Weight_ * NumPar_now * Tn_;
				}
				else if (res < 0)
				{
					Flux_GridBoundry[i][0] += Weight_ * NumPar_now;
					EneFlux_GridBoundry[i][0] += Weight_ * NumPar_now * Tn_;
				}
			}
		}

	for (int i = 0; i < num_PFRBoundry; i++)
		if (get_line_intersection(X_[0], X_[1], X_new_[0], X_new_[1],
								  PFRBoundry[i][0], PFRBoundry[i][1],
								  PFRBoundry[i + 1][0], PFRBoundry[i + 1][1], &i_x,
								  &i_y))
		{
			res = (PFRBoundry[i + 1][0] - PFRBoundry[i][0]) *
					  (X_[1] - PFRBoundry[i][1]) -
				  (PFRBoundry[i + 1][1] - PFRBoundry[i][1]) *
					  (X_[0] - PFRBoundry[i][0]);
			if (res < 0)
			{
				Flux_PFRBoundry[i][0] += Weight_ * NumPar_now;
				EneFlux_PFRBoundry[i][0] += Weight_ * NumPar_now * Tn_;
			}
			else if (res > 0)
			{
				Flux_PFRBoundry[i][1] += Weight_ * NumPar_now;
				EneFlux_PFRBoundry[i][1] += Weight_ * NumPar_now * Tn_;
			}
			else
			{
				res = (PFRBoundry[i + 1][0] - PFRBoundry[i][0]) *
						  (X_new_[1] - PFRBoundry[i][1]) -
					  (PFRBoundry[i + 1][1] - PFRBoundry[i][1]) *
						  (X_new_[0] - PFRBoundry[i][0]);
				if (res > 0)
				{
					Flux_PFRBoundry[i][0] += Weight_ * NumPar_now;
					EneFlux_PFRBoundry[i][0] += Weight_ * NumPar_now * Tn_;
				}
				else if (res < 0)
				{
					Flux_PFRBoundry[i][1] += Weight_ * NumPar_now;
					EneFlux_PFRBoundry[i][1] += Weight_ * NumPar_now * Tn_;
				}
			}
		}
}*/

void Particle::FluxCal_Target()
{
}

void Particle::FluxOutput()
{
	ofstream out;
	string pathout;
	pathout = Outputpath + name_ + "_" + "n_Flux_Grid.txt";
	out.open(pathout);
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			out << i << "\t" << j << "\t";
			for (int k = 0; k < 4; k++)
				out << n_Flux_Grid_[i][j][k] << "\t";
			out << endl;
		}
	}
	out.close();

	pathout = Outputpath + name_ + "_" + "T_Flux_Grid.txt";
	out.open(pathout);
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			out << i << "\t" << j << "\t";
			for (int k = 0; k < 4; k++)
				out << T_Flux_Grid_[i][j][k] << "\t";
			out << endl;
		}
	}
	out.close();
}

void Particle::regression()
{
	X_[0] = GridIndex_old_ / N_radial;
	X_[0] = GridIndex_old_ % N_radial;
	Zone_ = Zone_old_;
}

void Particle::Find()
{
	GridIndex_old_ = XY_[0] * N_radial + XY_[1];
	int GridIndex_GRID;
	Zone_old_ = Zone_;
	double Z[250], Y[250];
	Point Start, End;
	int k = 0;
	int ifOK = 0;

	if (XY_[0] > 0 && XY_[0] < 97 && XY_[1] > 0 && XY_[1] < 37)
	{
		if (K_GRID)
		{
			if (Grid1.IfinIndex(XY_[0], XY_[1], X_new_[0], X_new_[1]))
			{
				GridIndex_ = XY_[0] * N_radial + XY_[1];
				GridIndex_GRID = GridIndex_;
				// ifOK = 1;
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				// std::cout << 1;
				return;
			}
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				Z[i] = Grid[XY_[0]][XY_[1]][i];
				Y[i] = Grid[XY_[0]][XY_[1]][i + 4];
			}
			Z[4] = Grid[XY_[0]][XY_[1]][0];
			Y[4] = Grid[XY_[0]][XY_[1]][0 + 4];
			if (IfinEdge(Z, Y, 4, X_new_[0], X_new_[1]))
			{
				// std::cout << "1";
				GridIndex_ = XY_[0] * N_radial + XY_[1];
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				return;
			}
		}
	}
	if (XY_[0] > 1 && XY_[0] < 96 && XY_[1] > 1 && XY_[1] < 36)
	{
		if (K_GRID)
		{
			int ifOK_GRID = 0;
			if (Grid1.IfinIndex(XY_[0], XY_[1], X_new_[0], X_new_[1]))
			{
				GridIndex_ = XY_[0] * N_radial + XY_[1];
				GridIndex_GRID = GridIndex_;
				ifOK = 1;
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] + 1, XY_[1], X_new_[0], X_new_[1]))
				{
					XY_[0] += 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 2;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] - 1, XY_[1], X_new_[0], X_new_[1]))
				{
					XY_[0] -= 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 3;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0], XY_[1] + 1, X_new_[0], X_new_[1]))
				{
					XY_[1] += 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 4;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0], XY_[1] - 1, X_new_[0], X_new_[1]))
				{
					XY_[1] -= 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 5;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] + 1, XY_[1] + 1, X_new_[0], X_new_[1]))
				{
					XY_[0] += 1;
					XY_[1] += 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 2;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] + 1, XY_[1] - 1, X_new_[0], X_new_[1]))
				{
					XY_[0] += 1;
					XY_[1] -= 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 2;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] - 1, XY_[1] + 1, X_new_[0], X_new_[1]))
				{
					XY_[0] -= 1;
					XY_[1] += 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 2;
				}
			}

			if (!ifOK)
			{
				if (Grid1.IfinIndex(XY_[0] - 1, XY_[1] - 1, X_new_[0], X_new_[1]))
				{
					XY_[0] -= 1;
					XY_[1] -= 1;
					GridIndex_ = XY_[0] * 38 + XY_[1];
					if (Charge_ > 0)
						Vchargefix();
					ifOK_GRID = 1;
					ifOK = 1;
					// std::cout << 2;
				}
			}

			if (ifOK == 1)
			{
				// std::cout << GridIndex_ << '\t';
				GridIndex_GRID = GridIndex_;
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				ifOK_GRID = 0;
				return;
			}
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				Z[i] = Grid[XY_[0] + 1][XY_[1]][i];
				Y[i] = Grid[XY_[0] + 1][XY_[1]][i + 4];
			}
			Z[4] = Grid[XY_[0] + 1][XY_[1]][0];
			Y[4] = Grid[XY_[0] + 1][XY_[1]][0 + 4];
			if (IfinEdge(Z, Y, 4, X_new_[0], X_new_[1]))
			{
				// std::cout << "2" << endl;
				XY_[0] += 1;
				GridIndex_ = XY_[0] * 38 + XY_[1];
				if (Charge_ > 0)
					Vchargefix();
				if (GridIndex_ != GridIndex_GRID)
				{
					if (K_GRID)
						std::cerr << "error1\t" << GridIndex_ << GridIndex_GRID << endl;
				}
				// std::cout << GridIndex_ << endl;
				ifOK = 1;
			}

			if (ifOK)
			{
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				return;
			}
			for (int i = 0; i < 4; i++)
			{
				Z[i] = Grid[XY_[0] - 1][XY_[1]][i];
				Y[i] = Grid[XY_[0] - 1][XY_[1]][i + 4];
			}
			Z[4] = Grid[XY_[0] - 1][XY_[1]][0];
			Y[4] = Grid[XY_[0] - 1][XY_[1]][0 + 4];
			if (IfinEdge(Z, Y, 4, X_new_[0], X_new_[1]))
			{
				// std::cout << "3" << endl;
				XY_[0] -= 1;
				GridIndex_ = XY_[0] * 38 + XY_[1];
				if (Charge_ > 0)
					Vchargefix();
				if (GridIndex_ != GridIndex_GRID)
					if (K_GRID)
						std::cerr << "error2\t" << GridIndex_ << GridIndex_GRID << endl;
				// std::cout << GridIndex_ << endl;
				ifOK = 1;
			}
			if (ifOK)
			{
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				return;
			}
			for (int i = 0; i < 4; i++)
			{
				Z[i] = Grid[XY_[0]][XY_[1] + 1][i];
				Y[i] = Grid[XY_[0]][XY_[1] + 1][i + 4];
			}
			Z[4] = Grid[XY_[0]][XY_[1] + 1][0];
			Y[4] = Grid[XY_[0]][XY_[1] + 1][0 + 4];
			if (IfinEdge(Z, Y, 4, X_new_[0], X_new_[1]))
			{
				// std::cout << "4" << endl;
				XY_[1] += 1;
				GridIndex_ = XY_[0] * 38 + XY_[1];
				if (Charge_ > 0)
					Vchargefix();
				if (GridIndex_ != GridIndex_GRID)
					if (K_GRID)
						std::cerr << "error3\t" << GridIndex_ << GridIndex_GRID << endl;
				// std::cout << GridIndex_ << endl;
				ifOK = 1;
			}
			if (ifOK)
			{
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				return;
			}
			for (int i = 0; i < 4; i++)
			{
				Z[i] = Grid[XY_[0]][XY_[1] - 1][i];
				Y[i] = Grid[XY_[0]][XY_[1] - 1][i + 4];
			}
			Z[4] = Grid[XY_[0]][XY_[1] - 1][0];
			Y[4] = Grid[XY_[0]][XY_[1] - 1][0 + 4];
			if (IfinEdge(Z, Y, 4, X_new_[0], X_new_[1]))
			{
				// std::cout << "5" << endl;
				XY_[1] -= 1;
				GridIndex_ = XY_[0] * 38 + XY_[1];
				if (Charge_ > 0)
					Vchargefix();
				if (GridIndex_ != GridIndex_GRID)
					if (K_GRID)
						std::cerr << "error4\t" << GridIndex_ << GridIndex_GRID << endl;
				// std::cout << GridIndex_ << endl;
				ifOK = 1;
			}

			if (ifOK)
			{
				if (XY_[1] >= 19 && XY_[1] <= 36)
					Zone_ = 3;
				else if (XY_[0] <= 24)
					Zone_ = 4;
				else if (XY_[0] <= 72)
					Zone_ = 2;
				else
					Zone_ = 5;
				return;
			}
		}
	}

	///
	for (int i = 25; i < 73; i++)
	{
		Z[k] = Grid[i][1][0];
		Y[k++] = Grid[i][1][4];
	}
	Z[k] = Z[0];
	Y[k] = Y[0];

	if (IfinEdge(Z, Y, k, X_new_[0], X_new_[1]))
	{
		Zone_ = 1;
		XY_[0] = -1;
		XY_[1] = -1;
		return;
	}

	k = 0;
	Start.PointIndex(25, 1);
	End.PointIndex(72, 18);
	Edge(Start.Index(), End.Index(), Z, Y, &k);
	if (IfinEdge(Z, Y, k, X_new_[0], X_new_[1]))
	{
		Zone_ = 2;
		if (ifOK)
		{
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (GridIndex_ != GridIndex_GRID)
				if (K_GRID)
					std::cerr << "error5\t" << GridIndex_ << GridIndex_GRID << endl;
			// std::cout << GridIndex_ << endl;
			if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
				Vchargefix();
			return;
		}
		XY_[0] = Xfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		Start.PointIndex(XY_[0], Start.Index()[1]);
		End.PointIndex(XY_[0], End.Index()[1]);
		XY_[1] = Yfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		GridIndex_ = XY_[0] * N_radial + XY_[1];
		// std::cout << GridIndex_ << endl;
		if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
			Vchargefix();
		return;
	}

	k = 0;
	Start.PointIndex(1, 19);
	End.PointIndex(96, 36);
	Edge(Start.Index(), End.Index(), Z, Y, &k);
	if (IfinEdge(Z, Y, k, X_new_[0], X_new_[1]))
	{
		Zone_ = 3;
		if (ifOK)
		{
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (GridIndex_ != GridIndex_GRID)
				if (K_GRID)
					std::cerr << "error6\t" << GridIndex_ << GridIndex_GRID << endl;
			// std::cout << GridIndex_ << endl;
			if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
				Vchargefix();
			return;
		}
		XY_[0] = Xfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		Start.PointIndex(XY_[0], Start.Index()[1]);
		End.PointIndex(XY_[0], End.Index()[1]);
		XY_[1] = Yfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		GridIndex_ = XY_[0] * N_radial + XY_[1];
		// std::cout << GridIndex_ << endl;
		if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
			Vchargefix();
		return;
	}

	k = 0;
	Start.PointIndex(1, 1);
	End.PointIndex(24, 18);
	Edge(Start.Index(), End.Index(), Z, Y, &k);
	if (IfinEdge(Z, Y, k, X_new_[0], X_new_[1]))
	{
		Zone_ = 4;
		if (ifOK)
		{
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (GridIndex_ != GridIndex_GRID)
				if (K_GRID)
					std::cerr << "error7\t" << GridIndex_ / 38 << '\t' << GridIndex_GRID / 38 << '\t' << GridIndex_ % 38 << '\t' << GridIndex_GRID % 38 << endl;
			// std::cout << GridIndex_ << endl;
			if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
				Vchargefix();
			return;
		}
		XY_[0] = Xfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		Start.PointIndex(XY_[0], Start.Index()[1]);
		End.PointIndex(XY_[0], End.Index()[1]);
		XY_[1] = Yfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		GridIndex_ = XY_[0] * N_radial + XY_[1];
		// std::cout << GridIndex_ << endl;
		if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
			Vchargefix();
		return;
	}

	k = 0;
	Start.PointIndex(73, 1);
	End.PointIndex(96, 18);
	Edge(Start.Index(), End.Index(), Z, Y, &k);
	if (IfinEdge(Z, Y, k, X_new_[0], X_new_[1]))
	{
		Zone_ = 5;
		if (ifOK)
		{
			GridIndex_ = XY_[0] * N_radial + XY_[1];
			if (GridIndex_ != GridIndex_GRID)
				if (K_GRID)
					std::cerr << "error8\t" << GridIndex_ << GridIndex_GRID << endl;
			// std::cout << GridIndex_ << endl;
			if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
				Vchargefix();
			return;
		}
		XY_[0] = Xfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		Start.PointIndex(XY_[0], Start.Index()[1]);
		End.PointIndex(XY_[0], End.Index()[1]);
		XY_[1] = Yfind(Start.Index(), End.Index(), X_new_[0], X_new_[1]);
		GridIndex_ = XY_[0] * N_radial + XY_[1];
		// std::cout << GridIndex_ << endl;
		if (Charge_ > 0 && GridIndex_ != GridIndex_old_)
			Vchargefix();
		return;
	}

	/// judge if the particle is inside of the wall
	for (int i = 0; i <= Grid1.Wall_num(); i++)
	{
		Z[i] = Grid1.Wall(i, 0);
		Y[i] = Grid1.Wall(i, 1);
	}
	if (IfinEdge(Z, Y, Grid1.Wall_num(), X_new_[0], X_new_[1]))
	{
		Zone_ = 6;
		XY_[0] = -1;
		XY_[1] = -1;
		GridIndex_ = XY_[0] * N_radial + XY_[1];
	}
	else
	{
		Zone_ = 7;
		XY_[0] = -1;
		XY_[1] = -1;
		GridIndex_ = XY_[0] * N_radial + XY_[1];
		return;
	}
}

void Particle::setmass(double m) { mass_ = m; }

void Particle::setfate(int i, int FATE_temp, int sourcePar)
{
	fate_[i] = FATE_temp;
	sourcePar_[i] = sourcePar;
	if (i == 0)
		for (int j = 0; j < 2; j++)
			sourceGrid_[j] = XY_[j];
}

void Particle::CalWeight1(int num)
{
	if (MeshMode == 1)
		for (int i = 0; i < N_poloidal; i++)
		{
			for (int j = 0; j < N_radial; j++)
			{
				if (K_DT)
				{
					if (this == &D)
					{
						NumPar_Grid_[i][j] = n_D_1[i][j] * Rec_[1].cs(i, j) * Volume[i][j] / (double)num;
						// std::cout << NumPar_Grid_[i][j] << endl;
						if (NumPar_Grid_[i][j] != 0)
							Weight_Grid_[i][j] = 1.;
						else
							Weight_Grid_[i][j] = 0.;
						// NumPar_sum_Grid_ += NumPar[i][j];
					}
					if (this == &T)
					{
						NumPar_Grid_[i][j] = n_T_1[i][j] * Rec_[1].cs(i, j) * Volume[i][j] / (double)num;
						// std::cout << NumPar_Grid_[i][j] << endl;
						if (NumPar_Grid_[i][j] != 0)
							Weight_Grid_[i][j] = 1.;
						else
							Weight_Grid_[i][j] = 0.;
						// NumPar_sum_Grid_ += NumPar[i][j];
					}
				}
				else
				{
					if (MeshMode == 1)
					{
						if (this == &D)
							NumPar_Grid_[i][j] = n_D_1[i][j] * Rec_[1].cs(i, j) * Volume[i][j] / (double)num;
						if (this == &T)
							NumPar_Grid_[i][j] = n_T_1[i][j] * Rec_[1].cs(i, j) * Volume[i][j] / (double)num;
						// std::cout << NumPar_Grid_[i][j] << endl;
						if (NumPar_Grid_[i][j] != 0)
							Weight_Grid_[i][j] = 1.;
						else
							Weight_Grid_[i][j] = 0.;
						// NumPar_sum_Grid_ += NumPar[i][j];
					}
				}
			}
		}
	else if (MeshMode == 3)
	{
		for (int i = 0; i < Grid4.num_tris(); i++)
		{
			if (Grid4.if_in_plasmagrid(i))
			{

				if (this == &D)
				{
					Tri_NumPar_Grid_[i] = n_D_1[Tri_B2_[i][0]][Tri_B2_[i][1]] * Rec_[1].cs(Tri_B2_[i][0], Tri_B2_[i][1]) * Grid4.triVolume(i) / (double)num;
				}
				if (this == &H)
				{
					Tri_NumPar_Grid_[i] = n_H_1[Tri_B2_[i][0]][Tri_B2_[i][1]] * Rec_[1].cs(Tri_B2_[i][0], Tri_B2_[i][1]) * Grid4.triVolume(i) / (double)num;
				}
				if (this == &T)
				{
					Tri_NumPar_Grid_[i] = n_T_1[Tri_B2_[i][0]][Tri_B2_[i][1]] * Rec_[1].cs(Tri_B2_[i][0], Tri_B2_[i][1]) * Grid4.triVolume(i) / (double)num;
				}
				if (Tri_NumPar_Grid_[i] != 0)
					Weight_Grid_[Tri_B2_[i][0]][Tri_B2_[i][1]] = 1.;
				else
					Weight_Grid_[Tri_B2_[i][0]][Tri_B2_[i][1]] = 0.;
			}
		}
	}
}

void Particle::CalWeight2(double NumPar[76], int num)
{
	NumPar_sum_Target_ = 0;
	for (int i = 0; i < 2 * N_radial; i++)
	{
		// std::cout << NumPar[i] << endl;
		NumPar_sum_Target_ += NumPar[i];
	}
	for (int i = 0; i < 2 * N_radial; i++)
	{
		Num_Target_[i] = (double)num / NumPar_sum_Target_ * NumPar[i];
		// std::cout << Num_Target_[i] << endl;

		if (Num_Target_[i] < 50 && Num_Target_[i] != 0 && NumPar[i] != 0)
		{
			Num_Target_[i] = 50;
		}

		if (Num_Target_[i] != 0)
			NumPar_Target_[i] = NumPar[i] / (double)Num_Target_[i];
		// else
		//	  NumPar_Target_[i] = 0;
		if (Num_Target_[i] != 0)
			Weight_Target_[i] = 1.; // NumPar[i] / NumPar_sum_Target_ / num;
	}
}

int Particle::Num_Target(int i) { return Num_Target_[i]; }

void Edge(int Start[2], int End[2], double a[], double b[], int *k)
{
	a[*k] = Grid[Start[0]][Start[1]][0];
	b[(*k)++] = Grid[Start[0]][Start[1]][4];
	for (int i = Start[0]; i <= End[0]; i++)
	{
		a[*k] = Grid[i][Start[1]][1];
		b[(*k)++] = Grid[i][Start[1]][5];
	}
	for (int j = Start[1]; j <= End[1]; j++)
	{
		a[*k] = Grid[End[0]][j][2];
		b[(*k)++] = Grid[End[0]][j][6];
	}
	for (int i = End[0]; i >= Start[0]; i--)
	{
		a[*k] = Grid[i][End[1]][3];
		b[(*k)++] = Grid[i][End[1]][7];
	}
	for (int j = End[1]; j >= Start[1]; j--)
	{
		a[*k] = Grid[Start[0]][j][0];
		b[(*k)++] = Grid[Start[0]][j][4];
	}
	(*k)--;
}

bool IfinEdge(double a[], double b[], int k, double x, double y)
{
	int count = 0;
	double secx, secy;
	for (int i = 0; i < k; i++)
	{
		if (min(b[i], b[i + 1]) < y && max(b[i], b[i + 1]) > y)
		{
			if (max(a[i], a[i + 1]) < x)
			{
				count += 1;
			}
			else if (a[i] < x || a[i + 1] < x)
			{
				if (get_line_intersection(x, y, -1, y, a[i], b[i], a[i + 1], b[i + 1],
										  &secx, &secy))
					count += 1;
			}
		}
	}
	if (count % 2 == 1)
		return true;
	else
		return false;
}

int Xfind(int Start[], int End[], double x, double y)
{
	if (End[0] == Start[0])
		return Start[0];
	double a[250], b[250];
	int k = 0, Xmid[2] = {(End[0] + Start[0]) / 2, End[1]};
	Edge(Start, Xmid, a, b, &k);
	if (IfinEdge(a, b, k, x, y))
		End[0] = Xmid[0];
	else
		Start[0] = Xmid[0] + 1;
	return Xfind(Start, End, x, y);
}

int Yfind(int Start[], int End[], double x, double y)
{
	if (End[1] == Start[1])
		return Start[1];
	double a[250], b[250];
	int k = 0, Ymid[2] = {End[0], (End[1] + Start[1]) / 2};
	Edge(Start, Ymid, a, b, &k);
	if (IfinEdge(a, b, k, x, y))
		End[1] = Ymid[1];
	else
		Start[1] = Ymid[1] + 1;
	return Yfind(Start, End, x, y);
}

void Vinit_cosine(std::vector<double> &Vt, double Rcos, double Rsin)
{
	do
	{
		double ztheta = 6.283185307 * Random();
		double zcost = cos(ztheta);
		double zsint = sin(ztheta);
		double zpsi = Random() * 3.1415927;
		double zcosp = abs(cos(zpsi));
		double zsinp = sin(zpsi);
		double Vx = zsint * zsinp;
		double Vy = zcosp;
		double Vz = zsinp * zcost;
		// std::cout << pow(Vx, 2) << ", " << pow(Vy, 2) << ", " << pow(Vz, 2) << "  ";

		Vt[0] = Rcos * Vx - Rsin * Vy;
		Vt[1] = Rsin * Vx + Rcos * Vy;
		Vt[2] = Vz;
		// std::cout << pow(Vt[0], 2) << "," << pow(Vt[1], 2) << ", " << pow(Vt[2], 2) << endl;
	} while (Vt[2] == 1 || Vt[0] < -3e8);
}

int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
						  double p2_x, double p2_y, double p3_x, double p3_y,
						  double *i_x, double *i_y)
{
	double s02_x, s02_y, s10_x, s10_y, s32_x, s32_y, s_numer, t_numer, denom, t;
	s10_x = p1_x - p0_x;
	s10_y = p1_y - p0_y;
	s32_x = p3_x - p2_x;
	s32_y = p3_y - p2_y;

	denom = s10_x * s32_y - s32_x * s10_y;
	if (denom == 0)
		return 0; // Collinear
	int denomPositive;
	if (denom > 0)
	{
		denomPositive = 1;
	}
	else
	{
		denomPositive = 0;
	}

	s02_x = p0_x - p2_x;
	s02_y = p0_y - p2_y;
	s_numer = s10_x * s02_y - s10_y * s02_x;
	t_numer = s32_x * s02_y - s32_y * s02_x;
	t = t_numer / denom;
	if (i_x != NULL)
		*i_x = p0_x + (t * s10_x);
	if (i_y != NULL)
		*i_y = p0_y + (t * s10_y);

	if ((s_numer < 0) == denomPositive)
		return 0; // No collision

	if ((t_numer < 0) == denomPositive)
		return 0; // No collision

	if (fabs(s_numer) > fabs(denom) || fabs(t_numer) > fabs(denom))
		return 0; // No collision
				  // Collision detected

	return 1;
}

PartoPar::PartoPar() { ifstore_ = 0; }

bool PartoPar::ifstore() { return ifstore_; }

bool PartoPar::ifChange() { return ifChange_; }

void PartoPar::store(Particle *PP, double factor, Particle *QQ, std::vector<double> &V,
					 int SOURCE, int toCollPar, int Charge)
{
	Particle PQ_temp;
	toV_.push(V);
	toCollPar_.push(toCollPar);
	toSource_.push(SOURCE);
	PQ_temp.Particlefrom(PP, factor, Charge);
	PQ_.push(PQ_temp);
	Q_.push(QQ);
	// Q_[Q_.size() - 1]->Particlefrom(PP, factor, Charge);
	// Q_[Q_.size() - 1]->setfate(0, SOURCE, toCollPar);
	if (K_test1)
	{
		std::cout << "CX Particle: " << Q_.back()->X(0) << ',' << Q_.back()->X(1) << '\t';
		std::cout << Q_.back()->X_new(0) << ',' << Q_.back()->X_new(1) << endl;
	}
	ifstore_ += 1;
}

void PartoPar::ParChange()
{
	if (ifstore_)
	{
		// std::cout << "Q.size()= " << Q_.size() << endl;
		while (!Q_.empty())
		{
			// std::cout << "Q_.size: " << Q_.size() << endl;
			Q_.front()->Particlefrom(&PQ_.front(), 1, PQ_.front().Charge());
			for (int i = 0; i < 3; i++)
				Q_.front()->SetV(i, toV_.front()[i]);
			// std::cout << toSource_.front() << '\t' << toCollPar_.front() << endl;
			Q_.front()->setfate(0, toSource_.front(), toCollPar_.front());
			if (toSource_.front() == 1 || toSource_.front() == 9)
				Q_.front()->CalTn();
			else
				Q_.front()->Init(6);

			if (K_test1)
			{
				std::cout << "before flight: " << Q_.front()->name() << ' ' << Q_.front()->X(0) << ' ' << Q_.front()->X(1) << ' ';
				std::cout << Q_.front()->X_new(0) << ' ' << Q_.front()->X_new(1) << ' ' << "Charge: " << Q_.front()->Charge() << endl;
			}

			Push(Q_.front());
			/// delete the Particle

			PQ_.pop();
			Q_.pop();
			toV_.pop();
			toCollPar_.pop();
			toSource_.pop();

			/*PQ_.shrink_to_fit();
			Q_.shrink_to_fit();
			toV_.shrink_to_fit();
			toCollPar_.shrink_to_fit();
			toSource_.shrink_to_fit();*/
		}
	}
	ifstore_ = 0;
}

void ParCollCar::initialize(int e_or_i, int num_trimesh)
{
	e_or_i_ = e_or_i;
	for (int i = 0; i < 98; i++)
		for (int j = 0; j < 38; j++)
		{
			cs_[i][j] = 1.e-74;
		}
	Cor_cs_ = 0.;
	num_trimesh_ = num_trimesh;
	for (int i = 0; i < num_trimesh_; i++)
	{
		Tri_cs_.push_back(1.e-74);
		Tri_Sn_.push_back(0.);
		Tri_Smu_.push_back(0.);
		Tri_Smu_0_.push_back(0.);
		Tri_Smu_1_.push_back(0.);
		Tri_Smu_2_.push_back(0.);
		Tri_Smu1_.push_back(0.);
		Tri_Smu2_.push_back(0.);
		Tri_SE_.push_back(0.);

		Tri_Num_n_.push_back(0.);
		Tri_Num_mu_.push_back(0.);
		Tri_Num_mu_0_.push_back(0.);
		Tri_Num_mu_1_.push_back(0.);
		Tri_Num_mu_2_.push_back(0.);
		Tri_Num_mu1_.push_back(0.);
		Tri_Num_mu2_.push_back(0.);
		Tri_Num_E_.push_back(0.);
		Tri_Pra_.push_back(0.);
		Tri_crossSection.push_back(0.);
		Tri_SE_n_.push_back(0.);
		Tri_Num_E_n_.push_back(0.);
	}
	V_relative_.resize(3);
}

/// @brief Calculate the cross-section of Charge exchange collision
/// @param ParColl data from ADAS database
/// @param Charge Charge of the particle who loss the electron
/// @param Par Particle who get the electron
/// @param cross_CX option of the cross-CX
void ParCollCar::ADASInput(ADAS *ParColl, int Charge, int Par, int cross_CX)
{
	if (ParColl == NULL)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				cs_[i][j] = 0;
			}
		return;
	}
	if (e_or_i_ == 1)
	{
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				if (Par == 2)
				{
					cs_[i][j] = ne[i][j] * ParColl->cal(Te[i][j], ne[i][j], Charge);
				}
				else if (Par == 3)
				{
					cs_[i][j] = ne[i][j] * ParColl->cal(Te[i][j], ne[i][j], Charge);
				}
				else
					cs_[i][j] = ne[i][j] * ParColl->cal(Te[i][j], ne[i][j], Charge);
				// cs_[i][j] = ne[i][j] * ParColl->cal(Te[i][j], ne[i][j], Charge);
			}
	}
	if (e_or_i_ == 2)
	{
		for (int i = 0; i < N_poloidal; i++)
		{
			for (int j = 0; j < N_radial; j++)
			{
				if (Par == 1)
				{
					cs_[i][j] = n_H_1[i][j] * ParColl->cal(Ti[i][j], n_H_1[i][j], Charge);
				}
				else if (Par == 2)
				{
					if (cross_CX == 1)
						cs_[i][j] = n_D_1[i][j] * ParColl->cal(Ti[i][j] / 2., n_D_1[i][j], Charge);
					else
						cs_[i][j] = n_D_1[i][j] * ParColl->cal(Ti[i][j], n_D_1[i][j], Charge);
				}
				else if (Par == 3)
				{
					if (cross_CX == 1)
						cs_[i][j] = n_T_1[i][j] * ParColl->cal(Ti[i][j] / 3., n_T_1[i][j], Charge);
					else
						cs_[i][j] = n_T_1[i][j] * ParColl->cal(Ti[i][j], n_T_1[i][j], Charge);
				}
				else
					std::cerr << "something is wrong here!" << endl;

				/*if (cross_CX == 1)
					cs_[i][j] /= Ratio_DT_Coll;*/
				// std::cout << cs_[i][j] << '\t';
				//  cs_[i][j] = ni[i][j] * ParColl->cal(Te[i][j], ne[i][j], Charge);
			}
			// std::cout << endl;
		}
	}
	Cor_cs_ = ne_core * ParColl->cal(Te_core, ne_core, Charge);
}

void ParCollCar::EIRENEInput(EIRENE *ParColl)
{
	std::cerr << "ni in this function should be set." << endl;
	/*if (e_or_i_ == 1)
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				cs_[i][j] = ne[i][j] * ParColl->cal(ne[i][j], Te[i][j]);
			}
	if (e_or_i_ == 2)
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				cs_[i][j] = ni[i][j] * ParColl->cal(ni[i][j], Te[i][j]);
			}*/
}

void ParCollCar::Setcs_now(double cs)
{
	cs_now_ = cs;
}

void ParCollCar::SetTrics_cs(std::vector<std::vector<int>> &Tri_b2)
{
	for (int i = 0; i < Tri_b2.size(); i++)
	{
		if (Tri_b2[i][0] != -1)
		{
			Tri_cs_[i] = cs_[Tri_b2[i][0]][Tri_b2[i][1]];
		}
	}
}

double ParCollCar::cs(int XY[3], double Zone)
{
	// std::cout << XY[0] <<
	if (Zone < 6 && Zone > 1)
	{
		return cs_[XY[0]][XY[1]];
	}
	if (Zone >= 6)
		return 1.e-74;
	if (Zone == 1)
		return Cor_cs_;
	else
		throw std::logic_error("Invalid Zone_ in cs()");
}

double ParCollCar::Cor_cs()
{
	return Cor_cs_;
}

double ParCollCar::cs(int i, int j)
{
	return cs_[i][j];
}

double ParCollCar::cs(int i)
{
	return Tri_cs_[i];
}

void ParCollCar::SetCor_cs(double Cor_cs)
{
	Cor_cs_ = Cor_cs;
}

double ParCollCar::cs_now()
{
	return cs_now_;
}

void ParCollCar::Clear()
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			Sn_[i][j] = 0.;
			Smu_[i][j] = 0.;
			SE_[i][j] = 0.;
			Pra_[i][j] = 0.;
			// std::cout << i << ", " << j << endl;
			Num_n_[i][j] = 0.;
			Num_mu_[i][j] = 0.;
			Num_E_[i][j] = 0.;
		}
	}

	for (int i = 0; i < num_trimesh_; i++)
	{
		Tri_Sn_[i] = 0.;
		Tri_Smu_[i] = 0.;
		Tri_SE_[i] = 0.;
		Tri_Pra_[i] = 0.;
		// std::cout << i << ", " << j << endl;
		Tri_Num_n_[i] = 0.;
		Tri_Num_mu_[i] = 0.;
		Tri_Num_E_[i] = 0.;
	}
}

void ParCollCar::Set_K(int i)
{
	K_ - i;
}

int ParCollCar::K()
{
	return K_;
}

void ParCollCar::Stat(int Charge, double *n[98][38], ADAS *PraADAS, const std::vector<std::vector<double>> &ni)
{
	if (MeshMode == 1)
	{
		for (int i = 1; i < N_poloidal - 1; i++)
			for (int j = 1; j < N_radial - 1; j++)
			{
				Sn_[i][j] = Num_n_[i][j];
				if (Num_n_[i][j] >= 1e-8)
				{
					Smu_[i][j] = Num_mu_[i][j] / Volume[i][j];
					SE_[i][j] = Num_E_[i][j] / Volume[i][j];
					SE_n_[i][j] = Num_E_n_[i][j] / Volume[i][j];
				}
			}
		if (PraADAS == NULL)
			return;
		for (int i = 0; i < N_poloidal; i++)
			for (int j = 0; j < N_radial; j++)
			{
				if (e_or_i_ == 1)
					Pra_[i][j] = n[i][j][Charge] * ne[i][j] * PraADAS->cal(Te[i][j], ne[i][j], Charge);
				if (e_or_i_ == 2)
					Pra_[i][j] = n[i][j][Charge] * ni[i][j] * PraADAS->cal(Ti[i][j], ni[i][j], Charge);
			}
	}
	else if (MeshMode == 3)
	{
	}
}

void ParCollCar::Stat_Tri(int Charge, const std::vector<std::vector<double>> &n, ADAS *PraADAS, const std::vector<std::vector<double>> &ni)
{
	for (int i = 0; i < num_trimesh_; i++)
	{
		Tri_Sn_[i] = Tri_Num_n_[i];
		if (Tri_Num_n_[i] >= 1e-8)
		{
			Tri_Smu_[i] = Tri_Num_mu_[i] / Grid4.triVolume(i);
			Tri_SE_[i] = Tri_Num_E_[i] / Grid4.triVolume(i);
			Tri_SE_n_[i] = Tri_Num_E_n_[i] / Grid4.triVolume(i);
		}
	}
	if (PraADAS == NULL)
		return;
	for (int i = 0; i < num_trimesh_; i++)
	{
		if (Grid4.if_in_plasmagrid(i))
		{
			if (e_or_i_ == 1)
				Tri_Pra_[i] = n[i][Charge] * ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * PraADAS->cal(Te[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Charge);
			else if (e_or_i_ == 2)
				Tri_Pra_[i] = n[i][Charge] * ni[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * PraADAS->cal(Ti[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], ni[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Charge);
		}
	}
}

void ParCollCar::Stat1(int Charge, double *n[98][38], ADAS *PraADAS, const std::vector<std::vector<double>> &ni)
{
	for (int i = 1; i < N_poloidal - 1; i++)
		for (int j = 1; j < N_radial - 1; j++)
		{
			if (Num_n_[i][j] >= 1e-8)
			{
				Smu1_[i][j] = Num_mu1_[i][j] / Volume[i][j];
				Smu2_[i][j] = Num_mu2_[i][j] / Volume[i][j];
			}
		}
	if (PraADAS == NULL)
		return;
}

void ParCollCar::StatEIRENE(int Charge, double *n[98][38], EIRENE *PraEIRENE, double E, double **ni)
{
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
		{
			Sn_[i][j] = Num_n_[i][j];
			Smu_[i][j] = Num_mu_[i][j] / Volume[i][j];
			SE_[i][j] = Num_E_[i][j] / Volume[i][j];
			if (e_or_i_ == 1)
				Pra_[i][j] = n[i][j][Charge] * ne[i][j] * PraEIRENE->cal(ne[i][j], Te[i][j]) * qe - Num_n_[i][j] * E * 1.6e-19;
			if (e_or_i_ == 2)
				Pra_[i][j] = n[i][j][Charge] * ni[i][j] * PraEIRENE->cal(ni[i][j], Ti[i][j]) * qe;
		}
}

void ParCollCar::StatEIRENE_Tri(int Charge, std::vector<std::vector<double>> &n, EIRENE *PraEIRENE, double E, const std::vector<std::vector<double>> &ni)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		Tri_Sn_[i] = Tri_Num_n_[i];
		Tri_Smu_[i] = Tri_Num_mu_[i] / Grid4.triVolume(i);
		Tri_SE_[i] = Tri_Num_E_[i] / Grid4.triVolume(i);
		if (Grid4.if_in_plasmagrid(i))
		{
			if (e_or_i_ == 1)
				Tri_Pra_[i] = n[i][Charge] * ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * PraEIRENE->cal(ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Te[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)]) * qe - Tri_Num_n_[i] * E * 1.6e-19;
			if (e_or_i_ == 2)
				Tri_Pra_[i] = n[i][Charge] * ni[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * PraEIRENE->cal(ni[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Ti[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)]) * qe;
		}
	}
}

void ParCollCar::RecStat(int Charge, double *n[98][38], ADAS *PraADAS)
{
	for (int i = 1; i < N_poloidal - 1; i++)
		for (int j = 1; j < N_radial - 1; j++)
		{
			Sn_[i][j] = Num_n_[i][j];
			Smu_[i][j] = Num_mu_[i][j] / Volume[i][j];
			SE_[i][j] = Num_E_[i][j] / Volume[i][j];
			Pra_[i][j] = n[i][j][Charge] * ne[i][j] * PraADAS->cal(Te[i][j], ne[i][j], Charge - 1);
		}
}

void ParCollCar::RecStat_Tri(int Charge, const std::vector<std::vector<double>> &n, ADAS *PraADAS)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		Tri_Sn_[i] = Tri_Num_n_[i];
		Tri_Smu_[i] = Tri_Num_mu_[i] / Grid4.triVolume(i);
		Tri_SE_[i] = Tri_Num_E_[i] / Grid4.triVolume(i);
		Tri_Pra_[i] = n[i][Charge] * ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)] * PraADAS->cal(Te[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], ne[Grid4.b2_index(i, 0)][Grid4.b2_index(i, 1)], Charge - 1);
	}
}

void ParCollCar::OutSn(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << Sn_[i][j] << "\t";
		in << endl;
	}
}

void ParCollCar::OutSn_Tri(std::ofstream &in)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		in << Tri_Sn_[i] << endl;
	}
}

void ParCollCar::OutSmu(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << Smu_[i][j] << "\t";
		in << endl;
	}
}

void ParCollCar::OutSmu_Tri(std::ofstream &in)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		in << Tri_Smu_[i] << endl;
	}
}

void ParCollCar::OutSmu1(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << Smu1_[i][j] << "\t";
		in << endl;
	}
}
void ParCollCar::OutSmu2(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << Smu2_[i][j] << "\t";
		in << endl;
	}
}

void ParCollCar::OutSE(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << SE_[i][j] << "\t";
		in << endl;
	}
}

void ParCollCar::OutSE_Tri(std::ofstream &in)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		in << Tri_SE_[i] << endl;
	}
}

void ParCollCar::OutSE_n(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << SE_n_[i][j] << "\t";
		in << endl;
	}
}

void ParCollCar::OutPra(std::ofstream &in)
{
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
			in << Pra_[i][j] << " ";
		in << endl;
	}
}

void ParCollCar::OutPra_Tri(std::ofstream &in)
{
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		in << Tri_Pra_[i] << endl;
	}
}

double ParCollCar::Sn(int X, int Y)
{
	return Sn_[X][Y];
}

double ParCollCar::Smu(int X, int Y)
{
	return Smu_[X][Y];
}

double ParCollCar::SE(int X, int Y)
{
	return SE_[X][Y];
}

void ParCollCar::Setcs(int i, double cs)
{
	Tri_cs_[i] = cs;
}

void ParCollCar::Setcs(int i, int j, double cs)
{
	cs_[i][j] = cs;
}

void ParCollCar::n_Add(int XY[], double Var)
{
	Num_n_[XY[0]][XY[1]] += Var;
}

void ParCollCar::n_Add(int XY, double Var)
{
	Tri_Num_n_[XY] += Var;
}

void ParCollCar::Mu_Add(int XY[], double Var)
{
	Num_mu_[XY[0]][XY[1]] += Var;
}

void ParCollCar::Mu_Add(int XY, double Var0, double Var1, double Var2)
{
	Tri_Num_mu_0_[XY] += Var0;
	Tri_Num_mu_1_[XY] += Var1;
	Tri_Num_mu_2_[XY] += Var2;
}

void ParCollCar::Mu_Add(int XY, double Var)
{
	Tri_Num_mu_[XY] += Var;
}

void ParCollCar::Mu_Add1(int XY[], double Var)
{
	Num_mu1_[XY[0]][XY[1]] += Var;
}

void ParCollCar::Mu_Add2(int XY[], double Var)
{
	Num_mu2_[XY[0]][XY[1]] += Var;
}

void ParCollCar::E_Add(int XY[], double Var)
{
	Num_E_[XY[0]][XY[1]] += Var;
}

void ParCollCar::E_Add(int XY, double Var)
{
	Tri_Num_E_[XY] += Var;
}

void ParCollCar::E_n_Add(int XY[], double Var)
{
	Num_E_n_[XY[0]][XY[1]] += Var;
}

void ParCollCar::E_n_Add(int XY, double Var)
{
	Tri_Num_E_n_[XY] += Var;
}

void ParCollCar::SetPra(int i, int j, double Var)
{
	Pra_[i][j] = Var;
}

void ParCollCar::SetPra_Tri(int i, double Var)
{
	Tri_Pra_[i] = Var;
}

void ParCollCar::Pra_Add(int XY[], double Var)
{
	Pra_[XY[0]][XY[1]] += Var;
}

void ParCollCar::Pra_Add(int XY, double Var)
{
	Tri_Pra_[XY] += Var;
}

int ZonefromXY(int X, int Y)
{
	if (Y >= 19 && Y <= 36)
		return 3;
	else if (X <= 24)
		return 4;
	else if (X <= 72)
		return 2;
	else
		return 5;
}

int Particle::IfParticleOut(int i)
{
	double X1, Y1;
	if (i == 1 || i == 2 || i == 3)
	{
		X1 = Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i - 1, 0) - Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i, 0);
		Y1 = Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i - 1, 1) - Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i, 1);
	}
	else if (i == 4)
	{
		X1 = Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i - 1, 0) - Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 0);
		Y1 = Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(i - 1, 1) - Grid1.Plasma_Grid(XY_[0], XY_[1]).Grid_Point(0, 1);
	}
	if (V_[0] * Y1 - V_[1] * X1 > 0)
		return 0;
	else
		return 1;
}

int PointandLine(double line_x1, double line_y1, double line_x2, double line_y2, double point_x, double point_y)
{
	double x1, x2, y1, y2;
	x1 = line_x2 - line_x1;
	y1 = line_y2 - line_y1;
	x2 = point_x - line_x1;
	y2 = point_y - line_y1;
	if (x1 * y2 - x2 * y1 >= 0)
		return 1;
	else
		return 0;
}

double randomSign()
{
	if (Random() < 0.5)
	{
		return 1.;
	}
	else
	{
		return -1.;
	}
}

double Xrecycling(double rand_, int a, int b)
{
	if (a < 38)
	{
		if (b == 0)
			return (Grid[1][a][3] - Grid[1][a][0]) * rand_ + Grid[1][a][0];
		else
			return (Grid[1][a][7] - Grid[1][a][4]) * rand_ + Grid[1][a][4];
	}
	else
	{
		if (b == 0)
			return (Grid[96][a - 38][1] * Grid[96][a - 38][2]) * rand_ + Grid[96][a - 38][2];
		else
			return (Grid[96][a - 38][5] * Grid[96][a - 38][6]) * rand_ + Grid[96][a - 38][6];
	}
}

string fatename(int i)
{
	string Fate;
	if (i == -1)
		Fate = "All";
	else if (i == 1)
		Fate = "CX";
	else if (i == 2)
		Fate = "recycling";
	else if (i == 3)
		Fate = "puff";
	else if (i == 4)
		Fate = "rec";
	else if (i == 5)
		Fate = "ion";
	else if (i == 6)
		Fate = "mar";
	else if (i == 7)
		Fate = "Diss";
	else if (i == 8)
		Fate = "Diss2";
	else if (i == 9)
		Fate = "elastic";
	else if (i == 10)
		Fate = "DS0";
	else if (i == 11)
		Fate = "DS1";
	else if (i == 12)
		Fate = "DS2";
	else if (i == 13)
		Fate = "DS3";
	else if (i == 14)
		Fate = "Diss1";
	else if (i == 15)
		Fate = "reflect";
	else if (i == 16)
		Fate = "thermalrelease";
	return Fate;
}

string sourcename(int i)
{
	string Source;
	if (i == 0)
		Source = "no";
	else if (i == 1)
		Source = "target";
	else if (i == 2)
		Source = "puffport";
	else if (i == 3)
		Source = "D";
	else if (i == 4)
		Source = "D2";
	else if (i == 5)
		Source = "CD4";
	else if (i == 6)
		Source = "CD3";
	else if (i == 7)
		Source = "CD2";
	else if (i == 8)
		Source = "CD1";
	else if (i == 9)
		Source = "C";
	else if (i == 10)
		Source = "Ar";
	else if (i == 11)
		Source = "wall";
	else if (i == 12)
		Source = "SOL";
	else if (i == 13)
		Source = "H";
	else if (i == 14)
		Source = "H2";
	else if (i == 15)
		Source = "T";
	else if (i == 16)
		Source = "T2";
	else if (i == 17)
		Source = "e";
	else if (i == 18)
		Source = "Core";
	else if (i == 19)
		Source = "PlasmaBoundry";
	else if (i == 20)
		Source = "target_1";
	else
		Source = "not set";
	return Source;
}

/*void LXCatRead(double *Te_DATA, double *Cross_DATA, int Num_Array)
{
	ifstream input_temp;
	input_temp.open("path of file");
	for (int i = 0; i < Num_Array; i++)
	{
		input_temp >> Te_DATA[i] >> Cross_DATA[i];
	}
}

BoutReal interpolation_1D(BoutReal Te, double *Te_DATA, double *Cross_DATA, int Num_Array)
{
	int Num_Key;
	if (Te < Te_DATA[0])
	{
		std::cerr << "Te for this collison is too low" << endl;
		return Cross_DATA[0];
	}
	else if (Te > Cross_DATA[Num_Array])
	{
		std::cerr << "Te for this collison is too large" << endl;
		return Cross_DATA[Num_Array];
	}
	else
		Num_Key = BinarySearch(Te_DATA, 0, Num_Array, Te);
	return (Te * (Cross_DATA[Num_Key + 1] - Cross_DATA[Num_Key]) + Te_DATA[Num_Key] * Cross_DATA[Num_Key] - Te_DATA[Num_Key] * Cross_DATA[Num_Key + 1]) / (Te_DATA[Num_Key + 1] - Te_DATA[Num_Key]);
}*/
FluxTracker::FluxTracker() {}

FluxTracker::FluxTracker(int nRegion, int nBin, double eMin, double eMax,
						 bool trackMom, bool trackPow)
	: NREG_(nRegion), NBIN_(nBin),
	  E_MIN_(eMin), E_MAX_(eMax),
	  trackMom_(trackMom), trackPow_(trackPow),
	  flux_(nRegion, std::vector<double>(nRegion, 0.0)),
	  hist_(nRegion, std::vector<std::vector<double>>(nRegion, std::vector<double>(nBin, 0.0)))
{
	if (trackMom_)
		momFlux_.assign(NREG_,
						std::vector<std::array<double, 3>>(NREG_, {0.0, 0.0, 0.0}));
	if (trackPow_)
		powFlux_.assign(NREG_, std::vector<double>(NREG_, 0.0));
}

bool FluxTracker::FluxTraInit(int nRegion, int nBin,
							  double eMin, double eMax,
							  int nx, int ny,
							  bool trackMom, bool trackPow)
{
	// 防止重复初始化；如需支持重复，可先清空数据
	if (inited_)
		return false;

	// 1) 存参数
	NREG_ = nRegion;
	NBIN_ = nBin;
	E_MIN_ = eMin;
	E_MAX_ = eMax;
	trackMom_ = trackMom;
	trackPow_ = trackPow;

	// 2) 分配并置零
	flux_.assign(NREG_, std::vector<double>(NREG_, 0.0));

	hist_.assign(NREG_,
				 std::vector<std::vector<double>>(
					 NREG_, std::vector<double>(NBIN_, 0.0)));

	if (trackMom_)
		momFlux_.assign(NREG_,
						std::vector<std::array<double, 3>>(NREG_, {0.0, 0.0, 0.0}));

	if (trackPow_)
		powFlux_.assign(NREG_, std::vector<double>(NREG_, 0.0));

	NX_ = nx;
	NY_ = ny;
	gridEnabled_ = true;

	gridHist_.assign(NX_, std::vector<std::vector<double>>(NY_, std::vector<double>(NBIN_, 0.0)));

	if (trackMom_)
		gridMom_.assign(NX_, std::vector<std::array<double, 3>>(NY_, {0.0, 0.0, 0.0}));

	if (trackPow_)
		gridPow_.assign(NX_, std::vector<double>(NY_, 0.0));

	inited_ = true;
	return true;
}

void FluxTracker::accumulate(int from, int to, double E,
							 double vx, double vy, double vz, double w)
{
	assert(from >= 0 && from < NREG_ && to >= 0 && to < NREG_);
	flux_[from][to] += w;

	int bin = int((E - E_MIN_) / (E_MAX_ - E_MIN_) * NBIN_);
	if (bin < 0)
		bin = 0;
	else if (bin >= NBIN_)
		bin = NBIN_ - 1;
	hist_[from][to][bin] += w;

	if (trackMom_)
	{
		momFlux_[from][to][0] += w * vx;
		momFlux_[from][to][1] += w * vy;
		momFlux_[from][to][2] += w * vz;
	}
	if (trackPow_)
		powFlux_[from][to] += w * E; // eV·权重
}

void FluxTracker::write_H5(const std::string &fn, const std::string &grp) const
{ // 在最外层加一个 try-catch，捕获所有可能的 HDF5 写入错误（磁盘满、权限等）
	try
	{
		// 1. 创建全新的 HDF5 文件
		// 使用 H5F_ACC_TRUNC: 如果文件已存在，则截断清空它；如果不存在，则新建。
		H5File file(fn, H5F_ACC_TRUNC);

		// 2. 直接创建 Group
		// 因为文件是刚建的空文件，绝不会发生 "name already exists" 的错误
		Group g = file.createGroup("MyGroup");

		// ---- 粒子通量 (2D: NREG_ x NREG_) ----
		hsize_t d2[2]{(hsize_t)NREG_, (hsize_t)NREG_};
		DataSpace s2(2, d2);
		DataSet dsFlux = g.createDataSet("flux", PredType::NATIVE_DOUBLE, s2);

		std::vector<double> flatFlux;
		flatFlux.reserve(NREG_ * NREG_); // 【优化】：一次性分配好所需内存
		for (const auto &r : flux_)		 // 【优化】：使用 const auto& 避免拷贝
			flatFlux.insert(flatFlux.end(), r.begin(), r.end());
		dsFlux.write(flatFlux.data(), PredType::NATIVE_DOUBLE);

		// ---- 能谱 (3D: NREG_ x NREG_ x NBIN_) ----
		hsize_t d3[3]{(hsize_t)NREG_, (hsize_t)NREG_, (hsize_t)NBIN_};
		DataSpace s3(3, d3);
		DataSet dsHist = g.createDataSet("hist", PredType::NATIVE_DOUBLE, s3);

		std::vector<double> flatHist;
		flatHist.reserve(NREG_ * NREG_ * NBIN_);
		for (const auto &m : hist_)
			for (const auto &r : m)
				flatHist.insert(flatHist.end(), r.begin(), r.end());
		dsHist.write(flatHist.data(), PredType::NATIVE_DOUBLE);

		// ---- 动量 ----
		if (trackMom_)
		{
			DataSet dsMom = g.createDataSet("momFlux", PredType::NATIVE_DOUBLE, s3);
			std::vector<double> flatMom;
			flatMom.reserve(NREG_ * NREG_ * NBIN_);
			for (const auto &m : momFlux_)
				for (const auto &r : m)
					flatMom.insert(flatMom.end(), r.begin(), r.end());
			dsMom.write(flatMom.data(), PredType::NATIVE_DOUBLE);
		}

		// ---- 能量 ----
		if (trackPow_)
		{
			DataSet dsPow = g.createDataSet("powFlux", PredType::NATIVE_DOUBLE, s2);
			std::vector<double> flatPow;
			flatPow.reserve(NREG_ * NREG_);
			for (const auto &r : powFlux_)
				flatPow.insert(flatPow.end(), r.begin(), r.end());
			dsPow.write(flatPow.data(), PredType::NATIVE_DOUBLE);
		}

		// ---- 元数据 (Attributes) ----
		// 记录能量边界等标量信息
		g.createAttribute("E_min", PredType::NATIVE_DOUBLE, DataSpace())
			.write(PredType::NATIVE_DOUBLE, &E_MIN_);
		g.createAttribute("E_max", PredType::NATIVE_DOUBLE, DataSpace())
			.write(PredType::NATIVE_DOUBLE, &E_MAX_);

		// ---- 空间网格分布数据 ----
		if (gridEnabled_)
		{
			// 网格能谱 (3D: NX_ x NY_ x NBIN_)
			hsize_t dims3[3]{(hsize_t)NX_, (hsize_t)NY_, (hsize_t)NBIN_};
			DataSpace space3(3, dims3);
			DataSet dsGridHist = g.createDataSet("gridHist", PredType::NATIVE_DOUBLE, space3);

			std::vector<double> flatGridHist;
			flatGridHist.reserve(NX_ * NY_ * NBIN_);
			for (const auto &row : gridHist_)
				for (const auto &cell : row)
					flatGridHist.insert(flatGridHist.end(), cell.begin(), cell.end());
			dsGridHist.write(flatGridHist.data(), PredType::NATIVE_DOUBLE);

			// 网格动量分布
			if (trackMom_)
			{
				hsize_t dimsMom[3]{(hsize_t)NX_, (hsize_t)NY_, 3};
				DataSpace spaceMom(3, dimsMom);
				DataSet dsGridMom = g.createDataSet("gridMom", PredType::NATIVE_DOUBLE, spaceMom);

				std::vector<double> flatGridMom;
				flatGridMom.reserve(NX_ * NY_ * 3);
				for (const auto &row : gridMom_)
					for (const auto &cell : row)
						flatGridMom.insert(flatGridMom.end(), cell.begin(), cell.end());
				dsGridMom.write(flatGridMom.data(), PredType::NATIVE_DOUBLE);
			}

			// 网格能量分布
			if (trackPow_)
			{
				hsize_t dimsPow[2]{(hsize_t)NX_, (hsize_t)NY_};
				DataSpace spacePow(2, dimsPow);
				DataSet dsGridPow = g.createDataSet("gridPow", PredType::NATIVE_DOUBLE, spacePow);

				std::vector<double> flatGridPow;
				flatGridPow.reserve(NX_ * NY_);
				for (const auto &row : gridPow_)
					flatGridPow.insert(flatGridPow.end(), row.begin(), row.end());
				dsGridPow.write(flatGridPow.data(), PredType::NATIVE_DOUBLE);
			}
		}

		std::cout << ">>> 数据已成功保存至 HDF5 文件: " << fn << std::endl;
	}
	catch (const H5::Exception &ee)
	{
		// 捕获并打印出任何 HDF5 层面的错误
		std::cerr << "!!! HDF5 写入失败 !!!" << std::endl;
		ee.printErrorStack();
	}
}

void FluxTracker::accumulateGrid(int i, int j,
								 double E,
								 double vx, double vy, double vz,
								 double w)
{
	if (!gridEnabled_ || i < 0 || i >= NX_ || j < 0 || j >= NY_)
		return;

	// 能谱 bin 编号
	int bin = int((E - E_MIN_) / (E_MAX_ - E_MIN_) * NBIN_);
	if (bin < 0)
		bin = 0;
	else if (bin >= NBIN_)
		bin = NBIN_ - 1;

	// 能谱直方图统计
	gridHist_[i][j][bin] += w;

	// 动量通量
	if (trackMom_)
	{
		gridMom_[i][j][0] += w * vx;
		gridMom_[i][j][1] += w * vy;
		gridMom_[i][j][2] += w * vz;
	}

	// 功率通量
	if (trackPow_)
		gridPow_[i][j] += w * E; // eV·权重
}

void FluxTracker::normalizeByVolume()
{
	assert(gridEnabled_);

	// --- 能谱 ---
	for (int i = 0; i < NX_; ++i)
	{
		for (int j = 0; j < NY_; ++j)
		{
			double V = Volume[i][j];
			if (V <= 0.0)
				continue;

			for (int b = 0; b < NBIN_; ++b)
				gridHist_[i][j][b] /= V;

			// 动量
			if (trackMom_)
			{
				gridMom_[i][j][0] /= V;
				gridMom_[i][j][1] /= V;
				gridMom_[i][j][2] /= V;
			}
			// 功率
			if (trackPow_)
			{
				gridPow_[i][j] /= V;
			}
		}
	}
}

/// @brief Set and calculate the relative velocity
/// @param Vi0 V_x of ion
/// @param Vi1 V_y of ion
/// @param Vi2 V_z of ion
/// @param Vn0 V_x of neutral particle
/// @param Vn1 V_y of neutral particle
/// @param Vn2 V_z of neutral particle
void ParCollCar::Set_V_relative(double Vi0, double Vi1, double Vi2, double Vn0, double Vn1, double Vn2)
{
	V_relative_[0] = Vi0 - Vn0;
	V_relative_[1] = Vi1 - Vn1;
	V_relative_[2] = Vi2 - Vn2;
	V_relative_2_ = pow(V_relative_[0], 2) + pow(V_relative_[0], 2) + pow(V_relative_[0], 2);
	V_2_relative_ = pow(Vi0, 2) + pow(Vi1, 2) + pow(Vi2, 2) - (pow(Vn0, 2) + pow(Vn1, 2) + pow(Vn2, 2));
}

/// @brief Difference of speed between ion and neutral particle
double ParCollCar::V_relative(int i)
{
	return V_relative_[i];
}

/// @brief Square of speed difference
double ParCollCar::V_relative_2()
{
	return V_relative_2_;
}

/// @brief Difference of the square of speed
double ParCollCar::V_2_relative()
{
	return V_2_relative_;
}

// --- 3. 发射函数 ---
/**
 * @param v0, v1, v2 三角形的三个顶点 (逆时针顺序)
 * @param mode 1 for cos distribution, 2 for mirror reflection
 * @param incidentDir 入射方向 (仅 Specular 模式需要，建议归一化)
 */
void Particle::emit_particle(const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, int mode, Vector3 incidentDir)
{
	// --- A. 计算几何属性 ---
	Vector3 edge1 = v1 - v0;
	Vector3 edge2 = v2 - v0;
	Vector3 normal = edge1.cross(edge2).normalize();
	// std::cout << "normal: " << normal.x << "\t" << normal.y << "\t" << normal.z << endl;

	// --- B. 表面随机采样 (重心坐标) ---
	double r1 = Random();
	double r2 = Random();

	// 如果在四边形外侧，折叠回三角形内
	if (r1 + r2 > 1.0)
	{
		r1 = 1.0 - r1;
		r2 = 1.0 - r2;
	}
	// 计算位置 P = v0 + r1*edge1 + r2*edge2
	Vector3 pos = v0 + (edge1 * r1) + (edge2 * r2);

	// [关键] 沿法线偏移，防止粒子陷入表面 (Shadow Acne)
	double epsilon = 1e-5;
	pos = pos + (normal * epsilon);

	// --- C. 计算发射方向 ---
	Vector3 dir;

	if (mode == 2)
	{
		// 镜面反射: R = I - 2(N·I)N
		// 确保入射光是归一化的
		incidentDir = incidentDir.normalize();
		double dotNI = normal.dot(incidentDir);

		// 公式计算
		dir = incidentDir - (normal * (2.0 * dotNI));
	}
	else
	{
		// 余弦分布 (Cosine Weighted Lambertian)

		// 1. 生成局部球坐标
		double u1 = Random(); // 随机数1
		double u2 = Random(); // 随机数2

		// 半球采样公式
		double theta_r = 2.0 * 3.1415926535 * u1; // 方位角 phi
		double cos_phi = std::sqrt(1.0 - u2);	  // 极角 (z分量)
		double sin_phi = std::sqrt(u2);			  // 水平分量半径

		// 局部坐标向量 (假设法线是 Z 轴)
		double x_local = std::cos(theta_r) * sin_phi;
		double y_local = std::sin(theta_r) * sin_phi;
		double z_local = cos_phi;

		// 2. 构建局部坐标系 (TBN 矩阵)
		// 我们需要构建两个垂直于 Normal 的向量 (Tangent, Bitangent)
		Vector3 helper = (std::abs(normal.x) > 0.9) ? Vector3(0, 1, 0) : Vector3(1, 0, 0);

		Vector3 tangent = normal.cross(helper).normalize();
		Vector3 bitangent = normal.cross(tangent); // 已经垂直且归一化

		// 3. 转换到世界坐标: dir = T*x + B*y + N*z
		dir = (tangent * x_local) + (bitangent * y_local) + (normal * z_local);
	}

	// --- D. 数据回填到粒子结构体 (Vector3 -> double[]) ---
	dir = dir.normalize(); // 确保最终方向归一化

	// X_[0] = pos.x;
	// X_[1] = pos.y;
	// X_[2] = 0;
	V_[0] = dir.x;
	V_[1] = dir.y;
	V_[2] = dir.z;
}

Vector3 calculate_dissociation_velocity(Vector3 v_mol, double E_diss, double m_h)
{
	// 物理常数
	const double eV_to_Joule = 1.60218e-19;

	// 1. 计算质心系下的各向同性速度 u
	// 能量守恒: E_diss = 0.5 * m * u^2 + 0.5 * m * u^2 = m * u^2
	double u_mag = std::sqrt((E_diss * eV_to_Joule) / m_h);

	// 2. 蒙特卡洛采样随机方向 (球坐标)
	double xi1 = Random(); // 随机数 0-1
	double xi2 = Random(); // 随机数 0-1

	double phi = 2.0 * M_PI * xi1;
	double cos_theta = 2.0 * xi2 - 1.0;
	double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

	// 3. 构建局部坐标系 (Basis)
	// 我们需要一个基底 (n, e1, e2)，其中 n 平行于 v_mol
	Vector3 n = v_mol.normalize();

	// 处理分子静止的边缘情况
	if (v_mol.norm() == 0)
	{
		n = Vector3(0, 0, 1);
	}

	// 寻找垂直向量 e1
	// 技巧：找一个不与 n 平行的任意轴 aux，做叉乘
	Vector3 aux;
	if (std::abs(n.z) < 0.9)
	{
		aux = Vector3(0, 0, 1);
	}
	else
	{
		aux = Vector3(0, 1, 0);
	}

	Vector3 e1 = n.cross(aux).normalize();
	Vector3 e2 = n.cross(e1).normalize();

	// 4. 合成速度增量向量 vec_u
	// u_perp (垂直分量) = u * sin_theta * (cos_phi * e1 + sin_phi * e2)
	// u_para (平行分量) = u * cos_theta * n

	Vector3 u_vec = n * (u_mag * cos_theta) +
					e1 * (u_mag * sin_theta * std::cos(phi)) +
					e2 * (u_mag * sin_theta * std::sin(phi));

	// 5. 返回实验室系速度
	return v_mol + u_vec;
}