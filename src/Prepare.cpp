#include "Global.h"
#include "DBoundarySource.h"
#include "Particle.h"
#include "RecyclingFlightAllocation.h"
#include "TargetIncident.h"
// #include "sys/timeb.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <stdexcept>

void get_unit_vector_xy();
void get_center();
void get_dTi();
void CX_DT_Fix();

namespace
{
	double NormalProjectionFromAngle(double angle_degree)
	{
		const double kPi = 3.14159265358979323846;
		return std::abs(std::cos(angle_degree * kPi / 180.0));
	}

	double RadicalInverse(unsigned int index, unsigned int base)
	{
		double value = 0.;
		double factor = 1. / static_cast<double>(base);
		while (index > 0)
		{
			value += factor * static_cast<double>(index % base);
			index /= base;
			factor /= static_cast<double>(base);
		}
		return value;
	}

	void BuildDTargetIncidentAverages()
	{
		const int target_count = 2 * N_radial;
		DTargetFastProbability.assign(target_count, 0.);
		DTargetMeanIncidentEnergy.assign(target_count, 0.);
		DTargetMeanReflectedEnergy.assign(target_count, 0.);
		if (DTargetIncidentModel == 0 || K_DWTrimReflection != 1 || K_Wallelement != 1)
			return;

		double weighted_energy = 0.;
		double weighted_angle = 0.;
		double weighted_fast_probability = 0.;
		double weighted_sheath_factor = 0.;
		double source_weight = 0.;
		double inner_fast_probability = 0.;
		double inner_source_weight = 0.;
		double outer_fast_probability = 0.;
		double outer_source_weight = 0.;
		double minimum_sheath_factor = 0.;
		double maximum_sheath_factor = 0.;
		for (int target = 0; target < target_count; ++target)
		{
			const int grid_i = target < N_radial ? 0 : N_poloidal - 1;
			const int grid_j = target < N_radial ? target : target - N_radial;
			const double sheath_factor = DTargetSheathFactorAtCell(grid_i, grid_j);
			if (target == 0)
				minimum_sheath_factor = maximum_sheath_factor = sheath_factor;
			else
			{
				minimum_sheath_factor = std::min(minimum_sheath_factor, sheath_factor);
				maximum_sheath_factor = std::max(maximum_sheath_factor, sheath_factor);
			}
			if (DTargetIncidentModel == 2)
			{
				Tools::IncidentFluxSample incident;
				incident.energy_eV = EireneTargetIncidentEnergy(target);
				incident.angle_deg = 0.;
				Ei_Dion[target] = incident.energy_eV;
				DTargetMeanIncidentEnergy[target] = incident.energy_eV;
				DTargetIncidentAngle[target] = incident.angle_deg;
				DTargetFastProbability[target] =
					DTargetFastReflectionProbability(incident);
				if (DTargetFastProbability[target] > 0.)
					DTargetMeanReflectedEnergy[target] =
						D_W_Trim.MeanReflectedEnergy(
							incident.energy_eV, incident.angle_deg);
			}
			else
			{
				double energy_sum = 0.;
				double angle_sum = 0.;
				double fast_probability_sum = 0.;
				double reflected_energy_sum = 0.;
				for (int sample_index = 1;
					 sample_index <= DTargetIncidentSamples; ++sample_index)
				{
					auto incident = SampleDTargetIncidentFlux(
						target,
						RadicalInverse(sample_index, 2),
						RadicalInverse(sample_index, 3),
						RadicalInverse(sample_index, 5));
					const double fast_probability =
						DTargetFastReflectionProbability(incident);
					energy_sum += incident.energy_eV;
					angle_sum += incident.angle_deg;
					fast_probability_sum += fast_probability;
					if (fast_probability > 0.)
						reflected_energy_sum += fast_probability *
							D_W_Trim.MeanReflectedEnergy(
								incident.energy_eV, incident.angle_deg);
				}

				const double inverse_samples = 1. / DTargetIncidentSamples;
				DTargetMeanIncidentEnergy[target] = energy_sum * inverse_samples;
				DTargetIncidentAngle[target] = angle_sum * inverse_samples;
				DTargetFastProbability[target] =
					fast_probability_sum * inverse_samples;
				if (fast_probability_sum > 0.)
					DTargetMeanReflectedEnergy[target] =
						reflected_energy_sum / fast_probability_sum;
				Ei_Dion[target] = DTargetMeanIncidentEnergy[target];
			}

			const double target_weight = std::max(0., NumPar_wall_D[target]);
			weighted_energy += target_weight * DTargetMeanIncidentEnergy[target];
			weighted_angle += target_weight * DTargetIncidentAngle[target];
			weighted_fast_probability +=
				target_weight * DTargetFastProbability[target];
			if (target < N_radial)
			{
				inner_fast_probability +=
					target_weight * DTargetFastProbability[target];
				inner_source_weight += target_weight;
			}
			else
			{
				outer_fast_probability +=
					target_weight * DTargetFastProbability[target];
				outer_source_weight += target_weight;
			}
			weighted_sheath_factor += target_weight * sheath_factor;
			source_weight += target_weight;
		}

		if (source_weight > 0.)
		{
			std::cout << (DTargetIncidentModel == 1
				? "D target EIRENE NEMODS=7 incident averages: E="
				: "D target normal-monoenergetic sensitivity averages: E=")
					  << weighted_energy / source_weight
					  << " eV, angle=" << weighted_angle / source_weight
					  << " deg, fast reflection="
					  << weighted_fast_probability / source_weight
					  << " (IT="
					  << (inner_source_weight > 0.
							  ? inner_fast_probability / inner_source_weight
							  : 0.)
					  << ", OT="
					  << (outer_source_weight > 0.
							  ? outer_fast_probability / outer_source_weight
							  : 0.)
					  << ")"
					  << ", sheath factor="
					  << weighted_sheath_factor / source_weight
					  << " [" << minimum_sheath_factor << ", "
					  << maximum_sheath_factor << "]";
			std::cout << std::endl;
		}
	}

	void ApplyDTargetIonFluxProfile(const std::vector<double> &angle_B_with_target)
	{
		double inner_total = 0.0;
		double outer_total = 0.0;
		const int outer_poloidal = poloidalLastIndex();

		for (int i = 0; i < N_radial; ++i)
		{
			const double inner_density = std::max(0.0, n_D_1[1][i]);
			const double outer_density = std::max(0.0, n_D_1[outer_poloidal][i]);
			const double inner_source =
				inner_density * std::abs(ua_D_1[1][i]) *
				NormalProjectionFromAngle(angle_B_with_target[i]) * Grid4.Vol_Target(i);
			const double outer_source =
				outer_density * std::abs(ua_D_1[outer_poloidal][i]) *
				NormalProjectionFromAngle(angle_B_with_target[i + N_radial]) * Grid4.Vol_Target(i + N_radial);
			NumPar_wall_D[i] = inner_source;
			NumPar_wall_D[i + N_radial] = outer_source;
			inner_total += inner_source;
			outer_total += outer_source;
		}
		if (inner_total <= 0.0 || outer_total <= 0.0)
			throw std::runtime_error("D target ion-flux profile is zero; check nDion, ua_Dion and target geometry");
		std::cout << "D target source mode: SOLPS projected ion flux profile"
				  << " (IT D+ flux=" << inner_total
				  << ", OT D+ flux=" << outer_total << ")" << std::endl;
	}

	void BindEireneCollisionRates()
	{
		if (K_database == 2)
		{
			if (K_H)
				H.Ion_[0].EIRENEInput(&R2_1_5_H4, EireneDensitySource::Electron);
			if (K_D)
			{
				D.Ion_[0].EIRENEInput(&R2_1_5_H4, EireneDensitySource::Electron);
				D.Rec_[1].EIRENEInput(&R2_1_8_H4, EireneDensitySource::Electron);
			}
			if (K_T)
			{
				T.Ion_[0].EIRENEInput(&R2_1_5_H4, EireneDensitySource::Electron);
				T.Rec_[1].EIRENEInput(&R2_1_8_H4, EireneDensitySource::Electron);
			}
		}

		if (K_MarColl == 1)
		{
			if (K_H)
			{
				H2.MAR_[0].EIRENEInput(nullptr);
				H2.DS_[1][2].EIRENEInput(&R2_2_14_H4, EireneDensitySource::Electron);
			}
			if (K_D)
			{
				D2.MAR_[0].EIRENEInput(nullptr);
				D2.DS_[1][2].EIRENEInput(&R2_2_14_H4, EireneDensitySource::Electron);
			}
			if (K_T)
			{
				T2.MAR_[0].EIRENEInput(nullptr);
				T2.DS_[1][2].EIRENEInput(&R2_2_14_H4, EireneDensitySource::Electron);
			}
		}
		else if (K_MarColl == 2)
		{
			if (K_H)
			{
				H2.MAR_[0].EIRENEInput(&R3_2_3r_H4, EireneDensitySource::Electron);
				H2.DS_[1][2].EIRENEInput(nullptr);
			}
			if (K_D)
			{
				D2.MAR_[0].EIRENEInput(&R3_2_3r_H4, EireneDensitySource::Electron);
				D2.DS_[1][2].EIRENEInput(nullptr);
			}
			if (K_T)
			{
				T2.MAR_[0].EIRENEInput(&R3_2_3r_H4, EireneDensitySource::Electron);
				T2.DS_[1][2].EIRENEInput(nullptr);
			}
		}

		if (K_H)
		{
			H2.CX_[0].EIRENEInput(&R3_2_3_H2, EireneDensitySource::HIon);
			H2.Ion_[0].EIRENEInput(&R2_2_9_H4, EireneDensitySource::Electron);
			H2.Diss1_[0].EIRENEInput(&R2_2_5g_H4, EireneDensitySource::Electron);
			H2.Diss2_[0].EIRENEInput(&R2_2_10_H4, EireneDensitySource::Electron);
			H2.DS_[1][0].EIRENEInput(&R2_2_11_H4, EireneDensitySource::Electron);
			H2.DS_[1][1].EIRENEInput(&R2_2_12_H4, EireneDensitySource::Electron);
			H.R_with_H_[0].EIRENEInput(&R_H_H, EireneDensitySource::HNeutralTri, EireneArgumentMode::SameDensityTemperature);
			H.R_with_H2_[0].EIRENEInput(&R_H_H2, EireneDensitySource::H2NeutralTri, EireneArgumentMode::SameDensityTemperature);
			H2.R_with_H_[0].EIRENEInput(&R_H2_H, EireneDensitySource::HNeutralTri, EireneArgumentMode::SameDensityTemperature);
			H2.R_with_H2_[0].EIRENEInput(&R_H2_H2, EireneDensitySource::H2NeutralTri, EireneArgumentMode::SameDensityTemperature);
		}
		if (K_D)
		{
			D2.CX_[0].EIRENEInput(&R3_2_3_H2, EireneDensitySource::DIon);
			if (K_DT && K_CX_DT)
				D2.CX_DT_[0].EIRENEInput(&R3_2_3_H2, EireneDensitySource::TIon);
			D2.Ion_[0].EIRENEInput(&R2_2_9_H4, EireneDensitySource::Electron);
			D2.Diss1_[0].EIRENEInput(&R2_2_5g_H4, EireneDensitySource::Electron);
			D2.Diss2_[0].EIRENEInput(&R2_2_10_H4, EireneDensitySource::Electron);
			D2.DS_[1][0].EIRENEInput(&R2_2_11_H4, EireneDensitySource::Electron);
			D2.DS_[1][1].EIRENEInput(&R2_2_12_H4, EireneDensitySource::Electron);
			D.R_with_H_[0].EIRENEInput(&R_H_H, EireneDensitySource::DNeutralTri, EireneArgumentMode::SameDensityTemperature);
			D.R_with_H2_[0].EIRENEInput(&R_H_H2, EireneDensitySource::D2NeutralTri, EireneArgumentMode::SameDensityTemperature);
			D2.R_with_H_[0].EIRENEInput(&R_H2_H, EireneDensitySource::DNeutralTri, EireneArgumentMode::SameDensityTemperature);
			D2.R_with_H2_[0].EIRENEInput(&R_H2_H2, EireneDensitySource::D2NeutralTri, EireneArgumentMode::SameDensityTemperature);
		}
		if (K_T)
		{
			T2.CX_[0].EIRENEInput(&R3_2_3_H2, EireneDensitySource::TIon);
			if (K_DT && K_CX_DT)
				T2.CX_DT_[0].EIRENEInput(&R3_2_3_H2, EireneDensitySource::DIon);
			T2.Ion_[0].EIRENEInput(&R2_2_9_H4, EireneDensitySource::Electron);
			T2.Diss1_[0].EIRENEInput(&R2_2_5g_H4, EireneDensitySource::Electron);
			T2.Diss2_[0].EIRENEInput(&R2_2_10_H4, EireneDensitySource::Electron);
			T2.DS_[1][0].EIRENEInput(&R2_2_11_H4, EireneDensitySource::Electron);
			T2.DS_[1][1].EIRENEInput(&R2_2_12_H4, EireneDensitySource::Electron);
			T.R_with_H_[0].EIRENEInput(&R_H_H, EireneDensitySource::TNeutralTri, EireneArgumentMode::SameDensityTemperature);
			T.R_with_H2_[0].EIRENEInput(&R_H_H2, EireneDensitySource::T2NeutralTri, EireneArgumentMode::SameDensityTemperature);
			T2.R_with_H_[0].EIRENEInput(&R_H2_H, EireneDensitySource::TNeutralTri, EireneArgumentMode::SameDensityTemperature);
			T2.R_with_H2_[0].EIRENEInput(&R_H2_H2, EireneDensitySource::T2NeutralTri, EireneArgumentMode::SameDensityTemperature);
		}
	}
}

void Prepare()
{
	/// fix the CX cross-section
	// CX_DT_Fix();
	/*for (double i = 0.1; i < 200; i += 0.1)
	{
		std::cout << i << '\t' << CCD96_H.cal(i, 1e20) << '\t' << CCD96_D.cal(i, 1e20) << '\t' << CCD96_T.cal(i, 1e20) << '\t'
				  << CCD96_H.cal(i, 1e20) * Ratio_D_Coll << '\t' << CCD96_H.cal(i, 1e20) * Ratio_T_Coll << '\t' << CCD96_H.cal(i, 1e20) * Ratio_DT_Coll << endl;
	}*/
	/*double Ti_He[98][38], Ta_He[98][38];
	ifstream He_read;
	He_read.open("doc/He/Ti_96.txt");
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
			He_read >> Ti_He[i][j];
	He_read.close();
	He_read.open("doc/He/Ta_96.txt");
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
			He_read >> Ta_He[i][j];
	out.open("doc/He/Cross_CX_He_96.txt");
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			out << R5_3_1_H3.cal(1.5 * Ta_He[i][j], Ti_He[i][j]) << "\t";
		}
		out << endl;
	}
	out.close();
	std::cout << Ti_He[88][0] << "\t" << Ta_He[88][0] << endl;
	std::cout << R5_3_1_H3.cal(1.5 * 0.0960353288986987, 4.19076540653769) << endl;*/

	/*out.open("doc/Cross_Ela_D_H2_1.txt");
	for (double i = 0.01; i < 200; i += 0.01)
	{
		out << i << "\t";
		out << R0_3D_H3.cal(0.01, i) << "\t";
		out << R0_3D_H3.cal(0.1, i) << "\t";
		out << R0_3D_H3.cal(1, i) << "\t";
		out << R0_3D_H3.cal(10, i) << endl;
	}
	out.close();*/
	double deltaX, deltaY, deltaL;
	std::vector<Point> polygon;
	polygon.resize(4);
	srand((unsigned)time(NULL));

	std::cout << setprecision(8);
	std::cerr << setprecision(8);
	Xlog << setprecision(8);
	Vlog << setprecision(8);

	std::ofstream out;

	/// Calculate the Reduced mass
	DTmass = Dmass * Tmass / (Dmass + Tmass);
	D2Tmass = D2mass * Tmass / (D2mass + Tmass);
	DT2mass = Dmass * T2mass / (Dmass + T2mass);
	HH2mass = Hmass * H2mass / (Hmass + H2mass);
	DD2mass = Dmass * D2mass / (Dmass + D2mass);
	TT2mass = Tmass * T2mass / (Tmass + T2mass);
	H2H2mass = H2mass * H2mass / (H2mass + H2mass);
	D2D2mass = D2mass * D2mass / (D2mass + D2mass);
	T2T2mass = T2mass * T2mass / (T2mass + T2mass);

	// prepare for Reflection coeffition
	H_W.SetReflect(H_W_n_a, H_W_E_a, H_W_n_b, H_W_E_b, H_W_epsilon_L);
	D_W.SetReflect(D_W_n_a, D_W_E_a, D_W_n_b, D_W_E_b, D_W_epsilon_L);
	T_W.SetReflect(T_W_n_a, T_W_E_a, T_W_n_b, T_W_E_b, T_W_epsilon_L);
	H_C.SetReflect(H_C_n_a, H_C_E_a, H_C_n_b, H_C_E_b, H_C_epsilon_L);
	D_C.SetReflect(D_C_n_a, D_C_E_a, D_C_n_b, D_C_E_b, D_C_epsilon_L);
	T_C.SetReflect(T_C_n_a, T_C_E_a, T_C_n_b, T_C_E_b, T_C_epsilon_L);
	H_Be.SetReflect(H_Be_n_a, H_Be_E_a, H_Be_n_b, H_Be_E_b, H_Be_epsilon_L);
	D_Be.SetReflect(D_Be_n_a, D_Be_E_a, D_Be_n_b, D_Be_E_b, D_Be_epsilon_L);
	T_Be.SetReflect(T_Be_n_a, T_Be_E_a, T_Be_n_b, T_Be_E_b, T_Be_epsilon_L);

	const int target_stats_count = MeshMode == 3 ? Grid4.Num_Target() : 0;
	D_W_sputtered.init(Grid1.Wall_num() + target_stats_count, 200, 0.0, 1500, 50, 0.0, 1.0);
	D_W_sputtered.setSurfaceCounts(Grid1.Wall_num(), target_stats_count);
	D_W_sputtered.loadYieldTable("Inputfile/database/D_W_sputtered.data");
	/*for (int j = 0; j < 9; j++)
	{
		std::cout << "\t" << D_W_sputtered.Y_A(j);
	}
	std::cout << endl;
	for (int i = 0; i < 22; i++)
	{
		std::cout << D_W_sputtered.Y_E(i) << endl;
		for (int j = 0; j < 9; j++)
		{
			std::cout << D_W_sputtered.Y_table(i, j) << "\t";
		}
		std::cout << endl;
	}*/

	/// Calculate the cos and sin of Wall
	// std::cout << numWall << endl;
	/*for (int i = 0; i < Grid1.Wall_num(); i++)
	{
		// std::cout << i << '\t' << Wall[i][0] << '\t' << Wall[i][1] << endl;
		deltaX = Grid1.Wall(i + 1, 0) - Grid1.Wall(i, 0);
		deltaY = Grid1.Wall(i + 1, 1) - Grid1.Wall(i, 1);
		deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
		cosang[i] = -deltaX / deltaL;
		sinang[i] = -deltaY / deltaL;
	}
	for (int i = 25; i < 73; i++)
	{
		deltaX = Grid[i][0][1] - Grid[i][0][0];
		deltaY = Grid[i][0][5] - Grid[i][0][4];
		deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
		cosCore[i - 25] = deltaX / deltaL;
		sinCore[i - 25] = deltaY / deltaL;
	}*/
	Grid1.CalAngleTarget();
	/// calculate the area of Grid
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				polygon[(3 - k)].SetX(Grid[i][j][k]);
				polygon[(3 - k)].SetY(Grid[i][j][k + 4]);
			}
			Area[i][j] = Tools::polygonarea(polygon, 4);
		}

	/// calculate the epx, epy, erx, ery
	get_unit_vector_xy();
	get_center();
	get_dTi();

	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			B[i][j][0] = (Bpol[i][j] * epx[i][j] + Brad[i][j] * erx[i][j]) / Btot[i][j];
			B[i][j][1] = (Bpol[i][j] * epy[i][j] + Brad[i][j] * ery[i][j]) / Btot[i][j];
			B[i][j][2] = Btor[i][j] / Btot[i][j];
			vdEp[i][j] = Erad[i][j] * Btor[i][j] - Etor[i][j] * Brad[i][j];
			vdEr[i][j] = Etor[i][j] * Bpol[i][j] - Epol[i][j] * Btor[i][j];

			/// calculate the Vi
			Vi_D[i][j] = sqrt(3. * qe * Ti[i][j] / Dmass / 3.);
			if (K_T)
				Vi_T[i][j] = sqrt(3. * qe * Ti[i][j] / Tmass / 3.);
			// std::cout << Tools::sqr(B[i][j][0]) + Tools::sqr(B[i][j][1]) << " ";
			//  std::cout << erx[i][j] << '\t';
		}
		// std::cout << endl;
	}
	vector<int> K_collision(15, 0);
	for (int i = 0; i < K_collision.size(); i++)
		K_collision[i] = 0;

	std::vector<std::vector<int>> convert;
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		convert.push_back(std::vector<int>());
		for (int j = 0; j < 11; j++)
		{
			convert[i].push_back(Grid4.tris_[i].neigh[j]);
		}
	}

	/*std::ofstream out_temp("doc/convert.txt");
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		for (int j = 0; j < 11; j++)
		{
			out_temp << convert[i][j] << "\t";
		}
		out_temp << endl;
	}
	out_temp.close();*/
	if (K_H)
	{
		H.ParInit(K_collision, convert);
		H2.ParInit(K_collision, convert);
	}
	if (K_D)
	{
		D.ParInit(K_collision, convert);
		D2.ParInit(K_collision, convert);
	}
	if (K_T)
	{
		T.ParInit(K_collision, convert);
		T2.ParInit(K_collision, convert);
	}
	if (K_Methane)
	{
		CD4.ParInit(K_collision, convert);
		CD3.ParInit(K_collision, convert);
		CD2.ParInit(K_collision, convert);
		CD1.ParInit(K_collision, convert);
	}
	if (K_Methane || K_C)
	{
		C.ParInit(K_collision, convert);
	}
	if (K_Ar)
	{
		Ar.ParInit(K_collision, convert);
	}

	/// calculate the collision rate
	if (K_database == 1)
	{
		if (K_H)
		{
			H.calIonRateADAS(&SCD12_H);
			H.calCXRateADAS(&CCD96_H);
			H.calRecRateADAS(&ACD12_H);
		}

		if (K_D)
		{
			D.calIonRateADAS(&SCD12_H);
			D.calCXRateADAS(&CCD96_D);
			D.calRecRateADAS(&ACD12_H);
			if (K_CX_DT)
				D.calCXRateADAS(&CCD96_H, 1);
		}
		if (K_T)
		{
			T.calIonRateADAS(&SCD96_T);
			T.calCXRateADAS(&CCD96_T);
			T.calRecRateADAS(&ACD96_T);
			if (K_CX_DT)
				T.calCXRateADAS(&CCD96_H, 1);
		}
	}
	BindEireneCollisionRates();
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			if (false && K_database == 2)
			{
				if (K_H)
				{
					H.Ion_[0].Setcs(i, j, ne[i][j] * R2_1_5_H4.cal(ne[i][j], Te[i][j]));
				}

				if (K_D)
				{
					D.Ion_[0].Setcs(i, j, ne[i][j] * R2_1_5_H4.cal(ne[i][j], Te[i][j]));
					D.Rec_[1].Setcs(i, j, ne[i][j] * R2_1_8_H4.cal(ne[i][j], Te[i][j]));
				}
				if (K_T)
				{
					T.Ion_[0].Setcs(i, j, ne[i][j] * R2_1_5_H4.cal(ne[i][j], Te[i][j]));
					T.Rec_[1].Setcs(i, j, ne[i][j] * R2_1_8_H4.cal(ne[i][j], Te[i][j]));
				}
			}
			if (false && K_MarColl == 1)
			{
				if (K_H)
				{
					H2.MAR_[0].Setcs(i, j, 0.);
					H2.DS_[1][2].Setcs(i, j, ne[i][j] * R2_2_14_H4.cal(ne[i][j], Te[i][j]));
					H2.CX_[0].Setcs(i, j, n_H_1[i][j] * R3_2_3_H3.cal(3. / 2. * T_N, Te[i][j]));
				}
				if (K_D)
				{
					D2.MAR_[0].Setcs(i, j, 0.);
					D2.DS_[1][2].Setcs(i, j, ne[i][j] * R2_2_14_H4.cal(ne[i][j], Te[i][j] / 1.));
				}
				if (K_T)
				{
					T2.MAR_[0].Setcs(i, j, 0.);
					T2.DS_[1][2].Setcs(i, j, ne[i][j] * R2_2_14_H4.cal(ne[i][j], Te[i][j]));
				}
			}
			if (false && K_MarColl == 2)
			{
				if (K_H)
				{
					H2.MAR_[0].Setcs(i, j, ne[i][j] * R3_2_3r_H4.cal(ne[i][j], Te[i][j]));
					H2.DS_[1][2].Setcs(i, j, 0);
				}
				if (K_D)
				{
					D2.MAR_[0].Setcs(i, j, ne[i][j] * R3_2_3r_H4.cal(ne[i][j], Te[i][j]));
					D2.DS_[1][2].Setcs(i, j, 0);
				}
				if (K_T)
				{
					T2.MAR_[0].Setcs(i, j, ne[i][j] * R3_2_3r_H4.cal(ne[i][j], Te[i][j]));
					T2.DS_[1][2].Setcs(i, j, 0);
				}
			}
			if (false && K_H)
			{
				H2.Ion_[0].Setcs(i, j, ne[i][j] * R2_2_9_H4.cal(ne[i][j], Te[i][j]));
				H2.Diss1_[0].Setcs(i, j, ne[i][j] * R2_2_5g_H4.cal(ne[i][j], Te[i][j]));
				H2.Diss2_[0].Setcs(i, j, ne[i][j] * R2_2_10_H4.cal(ne[i][j], Te[i][j]));
				H2.DS_[1][0].Setcs(i, j, ne[i][j] * R2_2_11_H4.cal(ne[i][j], Te[i][j]));
				H2.DS_[1][1].Setcs(i, j, ne[i][j] * R2_2_12_H4.cal(ne[i][j], Te[i][j]));
			}

			if (false && K_D)
			{
				D2.Diss1_[0].Setcs(i, j, ne[i][j] * R2_2_5g_H4.cal(ne[i][j], Te[i][j] / 1.));
				D2.Ion_[0].Setcs(i, j, ne[i][j] * R2_2_9_H4.cal(ne[i][j], Te[i][j] / 1.));
				D2.Diss2_[0].Setcs(i, j, ne[i][j] * R2_2_10_H4.cal(ne[i][j], Te[i][j] / 1.));
				D2.DS_[1][0].Setcs(i, j, ne[i][j] * R2_2_11_H4.cal(ne[i][j], Te[i][j] / 1.));
				D2.DS_[1][1].Setcs(i, j, ne[i][j] * R2_2_12_H4.cal(ne[i][j], Te[i][j] / 1.));
			}

			if (false && K_T)
			{
				T2.Diss1_[0].Setcs(i, j, ne[i][j] * R2_2_5g_H4.cal(ne[i][j], Te[i][j]));
				T2.Ion_[0].Setcs(i, j, ne[i][j] * R2_2_9_H4.cal(ne[i][j], Te[i][j]));
				T2.Diss2_[0].Setcs(i, j, ne[i][j] * R2_2_10_H4.cal(ne[i][j], Te[i][j]));
				T2.DS_[1][0].Setcs(i, j, ne[i][j] * R2_2_11_H4.cal(ne[i][j], Te[i][j]));
				T2.DS_[1][1].Setcs(i, j, ne[i][j] * R2_2_12_H4.cal(ne[i][j], Te[i][j]));
			}

			if (K_Methane)
			{
				/// methane collision cross section calculation
				CD4.Ion_[0].Setcs(i, j, ne[i][j] * RM_2_1.cal(ne[i][j], Te[i][j]));
				CD4.CX_[0].Setcs(i, j, n_D_1[i][j] * RM_3_1_1.cal(T_D2_0[i][j], Ti[i][j]));
				CD4.Diss1_[0].Setcs(i, j, ne[i][j] * RM_2_11.cal(ne[i][j], Te[i][j]));
				CD4.Diss2_[0].Setcs(i, j, ne[i][j] * RM_2_5.cal(ne[i][j], Te[i][j]));
				CD4.DS_[1][0].Setcs(i, j, ne[i][j] * RM_2_15_1.cal(ne[i][j], Te[i][j]));
				CD4.DS_[1][1].Setcs(i, j, ne[i][j] * RM_2_15_2.cal(ne[i][j], Te[i][j]));
				CD4.DS_[1][2].Setcs(i, j, ne[i][j] * RM_2_19_1.cal(ne[i][j], Te[i][j]));
				CD4.DS_[1][3].Setcs(i, j, ne[i][j] * RM_2_19_2.cal(ne[i][j], Te[i][j]));
				// std::cout << D2.Ion_rate[i][j] << '\t';
				CD3.Ion_[0].Setcs(i, j, ne[i][j] * RM_2_2.cal(ne[i][j], Te[i][j]));
				CD3.CX_[0].Setcs(i, j, n_D_1[i][j] * RM_3_1_2.cal(T_D2_0[i][j], Ti[i][j]));
				CD3.Diss1_[0].Setcs(i, j, ne[i][j] * RM_2_12.cal(ne[i][j], Te[i][j]));
				CD3.Diss2_[0].Setcs(i, j, ne[i][j] * RM_2_8.cal(ne[i][j], Te[i][j]));
				CD3.DS_[1][0].Setcs(i, j, ne[i][j] * RM_2_16_1.cal(ne[i][j], Te[i][j]));
				CD3.DS_[1][1].Setcs(i, j, ne[i][j] * RM_2_16_2.cal(ne[i][j], Te[i][j]));
				CD3.DS_[1][2].Setcs(i, j, ne[i][j] * RM_2_20.cal(ne[i][j], Te[i][j]));

				CD2.Ion_[0].Setcs(i, j, ne[i][j] * RM_2_3.cal(ne[i][j], Te[i][j]));
				CD2.CX_[0].Setcs(i, j, n_D_1[i][j] * RM_3_1_3.cal(T_D2_0[i][j], Ti[i][j]));
				CD2.Diss1_[0].Setcs(i, j, ne[i][j] * RM_2_13.cal(ne[i][j], Te[i][j]));
				CD2.Diss2_[0].Setcs(i, j, ne[i][j] * RM_2_9.cal(ne[i][j], Te[i][j]));
				CD2.DS_[1][0].Setcs(i, j, ne[i][j] * RM_2_17_1.cal(ne[i][j], Te[i][j]));
				CD2.DS_[1][1].Setcs(i, j, ne[i][j] * RM_2_17_2.cal(ne[i][j], Te[i][j]));
				CD2.DS_[1][2].Setcs(i, j, ne[i][j] * RM_2_21.cal(ne[i][j], Te[i][j]));

				CD1.Ion_[0].Setcs(i, j, ne[i][j] * RM_2_4.cal(ne[i][j], Te[i][j]));
				CD1.CX_[0].Setcs(i, j, n_D_1[i][j] * RM_3_1_4.cal(T_D2_0[i][j], Ti[i][j]));
				CD1.Diss1_[0].Setcs(i, j, ne[i][j] * RM_2_14.cal(ne[i][j], Te[i][j]));
				CD1.Diss2_[0].Setcs(i, j, ne[i][j] * RM_2_10_1.cal(ne[i][j], Te[i][j]));
				CD1.Diss3_[0].Setcs(i, j, ne[i][j] * RM_2_10_2.cal(ne[i][j], Te[i][j]));
				CD1.DS_[1][0].Setcs(i, j, ne[i][j] * RM_2_18_1.cal(ne[i][j], Te[i][j]));
				CD1.DS_[1][1].Setcs(i, j, ne[i][j] * RM_2_18_2.cal(ne[i][j], Te[i][j]));
				CD1.DS_[1][2].Setcs(i, j, ne[i][j] * RM_2_22.cal(ne[i][j], Te[i][j]));
			}
		}
		// std::cout << endl;
	}
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		if (false && K_H)
		{
			if (n_H_0_Tri[i] > 1e8)
			{
				H.R_with_H_[0].Setcs(i, n_H_0_Tri[i] * R_H_H.cal(n_H_0_Tri[i], T_H_0_Tri[i]));
				H2.R_with_H_[0].Setcs(i, n_H_0_Tri[i] * R_H2_H.cal(n_H_0_Tri[i], T_H_0_Tri[i]));
			}
			if (n_H2_0_Tri[i] > 1e8)
			{
				H.R_with_H2_[0].Setcs(i, n_H2_0_Tri[i] * R_H_H2.cal(n_H2_0_Tri[i], T_H2_0_Tri[i]));
				H2.R_with_H2_[0].Setcs(i, n_H_0_Tri[i] * R_H2_H2.cal(n_H2_0_Tri[i], T_H2_0_Tri[i]));
			}
		}
		if (false && K_D)
		{
			if (n_D_0_Tri[i] > 1e8)
			{
				D.R_with_H_[0].Setcs(i, n_D_0_Tri[i] * R_H_H.cal(n_D_0_Tri[i], T_D_0_Tri[i]));
				D2.R_with_H_[0].Setcs(i, n_D_0_Tri[i] * R_H2_H.cal(n_D_0_Tri[i], T_D_0_Tri[i]));
			}
			if (n_D2_0_Tri[i] > 1e8)
			{
				D.R_with_H2_[0].Setcs(i, n_D2_0_Tri[i] * R_H_H2.cal(n_D2_0_Tri[i], T_D2_0_Tri[i]));
				D2.R_with_H2_[0].Setcs(i, n_D2_0_Tri[i] * R_H2_H2.cal(n_D2_0_Tri[i], T_D2_0_Tri[i]));
			}
		}
		if (false && K_T)
		{
			if (n_T_0_Tri[i] > 1e8)
			{
				T.R_with_H_[0].Setcs(i, n_T_0_Tri[i] * R_H_H.cal(n_T_0_Tri[i], T_T_0_Tri[i]));
				T2.R_with_H_[0].Setcs(i, n_T_0_Tri[i] * R_H2_H.cal(n_T_0_Tri[i], T_T_0_Tri[i]));
			}
			if (n_T2_0_Tri[i] > 1e8)
			{
				T.R_with_H2_[0].Setcs(i, n_T2_0_Tri[i] * R_H_H2.cal(n_T2_0_Tri[i], T_T2_0_Tri[i]));
				T2.R_with_H2_[0].Setcs(i, n_T2_0_Tri[i] * R_H2_H2.cal(n_T2_0_Tri[i], T_T2_0_Tri[i]));
			}
		}
	}
	/*ofstream out_temp("doc/cs_n-n.txt");
	for (int i = 0; i < Grid4.num_tris(); i++)
	{
		out_temp << i << "\t" << D.R_with_H_[0].cs(i) << "\t" << D.R_with_H2_[0].cs(i) << "\t" << D2.R_with_H_[0].cs(i) << "\t" << D2.R_with_H2_[0].cs(i) << endl;
	}
	out_temp.close();*/

	if (K_C || K_Methane)
	{
		C.calIonRateADAS(&SCD96_C);
		C.calRecRateADAS(&ACD96_C);
		if (K_CX_impurity)
			C.calCXRateADAS(&CCD96_C);
	}
	if (K_Ar)
	{
		Ar.calIonRateADAS(&SCD89_Ar);
		Ar.calRecRateADAS(&ACD89_Ar);
		if (K_CX_impurity)
			Ar.calCXRateADAS(&CCD89_Ar);
	}

	if (K_database == 1)
	{
		if (K_H)
		{
			H.Ion_[0].SetCor_cs(ne_core * SCD12_H.cal(Te_core, ne_core));
			H.CX_[0].SetCor_cs(ne_core * CCD96_D.cal(Te_core, ne_core));
		}
		if (K_D)
		{
			D.Ion_[0].SetCor_cs(ne_core * SCD12_H.cal(Te_core, ne_core));
			D.CX_[0].SetCor_cs(ne_core * CCD96_D.cal(Te_core, ne_core));
		}
		if (K_T)
		{
			T.Ion_[0].SetCor_cs(ne_core * SCD12_H.cal(Te_core, ne_core));
			T.CX_[0].SetCor_cs(ne_core * CCD96_D.cal(Te_core, ne_core));
		}
	}
	else if (K_database == 2)
	{
		if (K_H)
		{
			H.Ion_[0].SetCor_cs(ne_core * R2_1_5_H4.cal(ne_core, Te_core));
			H.CX_[0].SetCor_cs(ne_core * R3_1_8_H3.cal(3. / 2. * T_N, Te_core));
		}
		if (K_D)
		{
			D.Ion_[0].SetCor_cs(ne_core * R2_1_5_H4.cal(ne_core, Te_core));
			D.CX_[0].SetCor_cs(ne_core * R3_1_8_H3.cal(3. / 2. * T_N, Te_core));
		}
		if (K_T)
		{
			T.Ion_[0].SetCor_cs(ne_core * R2_1_5_H4.cal(ne_core, Te_core));
			T.CX_[0].SetCor_cs(ne_core * R3_1_8_H3.cal(3. / 2. * T_N, Te_core));
		}
	}
	if (K_MarColl == 1)
	{
		if (K_H)
		{
			H2.CX_[0].SetCor_cs(ne_core * R3_2_3_H2.cal(ne_core, Te_core));
			H2.MAR_[0].SetCor_cs(0.);
		}
		if (K_D)
		{
			D2.CX_[0].SetCor_cs(ne_core * R3_2_3_H2.cal(ne_core, Te_core));
			D2.MAR_[0].SetCor_cs(0.);
		}
		if (K_T)
		{
			T2.CX_[0].SetCor_cs(ne_core * R3_2_3_H2.cal(ne_core, Te_core));
			T2.MAR_[0].SetCor_cs(0.);
		}
	}
	else if (K_MarColl == 2)
	{
		if (K_H)
		{
			H2.CX_[0].SetCor_cs(0.);
			H2.MAR_[0].SetCor_cs(ne_core * R3_2_3r_H4.cal(ne_core, Te_core));
		}
		if (K_D)
		{
			D2.CX_[0].SetCor_cs(0.);
			D2.MAR_[0].SetCor_cs(ne_core * R3_2_3r_H4.cal(ne_core, Te_core));
		}
		if (K_T)
		{
			T2.CX_[0].SetCor_cs(0.);
			T2.MAR_[0].SetCor_cs(ne_core * R3_2_3r_H4.cal(ne_core, Te_core));
		}
	}
	if (K_H2_elastic == 1 && D2ElasticModel > 0)
	{
		if (K_H)
		{
			H2.Ela_[0].SetCor_cs(ne_core * R0_3T_H3.cal(3. / 2. * T_N, Te_core));
		}
		if (K_D)
		{
			D2.Ela_[0].SetCor_cs(ne_core * R0_3T_H3.cal(3. / 2. * T_N / 2., Te_core));
			if (K_DT)
			{
				D2.Ela_DT_[0].SetCor_cs(ne_core * R0_3T_H3.cal(3. / 2. * T_N / 2., Te_core));
			}
			else
			{
				D2.Ela_DT_[0].SetCor_cs(0.);
			}
		}
		if (K_T)
		{
			T2.Ela_[0].SetCor_cs(ne_core * R0_3T_H3.cal(3. / 2. * T_N / 3., Te_core));
			if (K_DT)
			{
				T2.Ela_DT_[0].SetCor_cs(ne_core * R0_3T_H3.cal(3. / 2. * T_N / 3., Te_core));
			}
			else
			{
				T2.Ela_DT_[0].SetCor_cs(0.);
			}
		}
	}
	else
	{
		if (K_H)
		{
			H2.Ela_[0].SetCor_cs(0.);
		}
		if (K_D)
		{
			D2.Ela_[0].SetCor_cs(0.);
			if (K_DT)
			{
				D2.Ela_DT_[0].SetCor_cs(0.);
			}
		}
		if (K_T)
		{
			T2.Ela_[0].SetCor_cs(0.);
			if (K_DT)
			{
				T2.Ela_DT_[0].SetCor_cs(0.);
			}
		}
	}
	if (K_H)
	{
		H2.Diss1_[0].SetCor_cs(ne_core * R2_2_5_H4.cal(ne_core, Te_core));
		H2.Ion_[0].SetCor_cs(ne_core * R2_2_9_H4.cal(ne_core, Te_core));
		H2.Diss2_[0].SetCor_cs(ne_core * R2_2_10_H4.cal(ne_core, Te_core));
		H2.DS_[1][0].SetCor_cs(ne_core * R2_2_11_H4.cal(ne_core, Te_core));
		H2.DS_[1][1].SetCor_cs(ne_core * R2_2_12_H4.cal(ne_core, Te_core));
		H2.DS_[1][2].SetCor_cs(ne_core * R2_2_14_H4.cal(ne_core, Te_core));
	}

	if (K_D)
	{
		D2.Diss1_[0].SetCor_cs(ne_core * R2_2_5_H4.cal(ne_core, Te_core));
		D2.Ion_[0].SetCor_cs(ne_core * R2_2_9_H4.cal(ne_core, Te_core));
		D2.Diss2_[0].SetCor_cs(ne_core * R2_2_10_H4.cal(ne_core, Te_core));
		D2.DS_[1][0].SetCor_cs(ne_core * R2_2_11_H4.cal(ne_core, Te_core));
		D2.DS_[1][1].SetCor_cs(ne_core * R2_2_12_H4.cal(ne_core, Te_core));
		D2.DS_[1][2].SetCor_cs(ne_core * R2_2_14_H4.cal(ne_core, Te_core));
	}

	if (K_T)
	{
		T2.Diss1_[0].SetCor_cs(ne_core * R2_2_5_H4.cal(ne_core, Te_core));
		T2.Ion_[0].SetCor_cs(ne_core * R2_2_9_H4.cal(ne_core, Te_core));
		T2.Diss2_[0].SetCor_cs(ne_core * R2_2_10_H4.cal(ne_core, Te_core));
		T2.DS_[1][0].SetCor_cs(ne_core * R2_2_11_H4.cal(ne_core, Te_core));
		T2.DS_[1][1].SetCor_cs(ne_core * R2_2_12_H4.cal(ne_core, Te_core));
		T2.DS_[1][2].SetCor_cs(ne_core * R2_2_14_H4.cal(ne_core, Te_core));
	}

	if (K_Methane)
	{
		CD4.Ion_[0].SetCor_cs(ne_core * RM_2_1.cal(ne_core, Te_core));
		CD4.CX_[0].SetCor_cs(ne_core * RM_3_1_1.cal(Te_core, Te_core));
		CD4.Diss1_[0].SetCor_cs(ne_core * RM_2_11.cal(ne_core, Te_core));
		CD4.Diss2_[0].SetCor_cs(ne_core * RM_2_5.cal(ne_core, Te_core));
		CD4.DS_[1][0].SetCor_cs(ne_core * RM_2_15_1.cal(ne_core, Te_core));
		CD4.DS_[1][1].SetCor_cs(ne_core * RM_2_15_2.cal(ne_core, Te_core));
		CD4.DS_[1][2].SetCor_cs(ne_core * RM_2_19_1.cal(ne_core, Te_core));
		CD4.DS_[1][3].SetCor_cs(ne_core * RM_2_19_2.cal(ne_core, Te_core));

		CD3.Ion_[0].SetCor_cs(ne_core * RM_2_2.cal(ne_core, Te_core));
		CD3.CX_[0].SetCor_cs(ne_core * RM_3_1_2.cal(Te_core, Te_core));
		CD3.Diss1_[0].SetCor_cs(ne_core * RM_2_12.cal(ne_core, Te_core));
		CD3.Diss2_[0].SetCor_cs(ne_core * RM_2_8.cal(ne_core, Te_core));
		CD3.DS_[1][0].SetCor_cs(ne_core * RM_2_16_1.cal(ne_core, Te_core));
		CD3.DS_[1][1].SetCor_cs(ne_core * RM_2_16_2.cal(ne_core, Te_core));
		CD3.DS_[1][2].SetCor_cs(ne_core * RM_2_20.cal(ne_core, Te_core));

		CD2.Ion_[0].SetCor_cs(ne_core * RM_2_3.cal(ne_core, Te_core));
		CD2.CX_[0].SetCor_cs(ne_core * RM_3_1_3.cal(Tn_core, Te_core));
		CD2.Diss1_[0].SetCor_cs(ne_core * RM_2_13.cal(ne_core, Te_core));
		CD2.Diss2_[0].SetCor_cs(ne_core * RM_2_9.cal(ne_core, Te_core));
		CD2.DS_[1][0].SetCor_cs(ne_core * RM_2_17_1.cal(ne_core, Te_core));
		CD2.DS_[1][1].SetCor_cs(ne_core * RM_2_17_2.cal(ne_core, Te_core));
		CD2.DS_[1][2].SetCor_cs(ne_core * RM_2_21.cal(ne_core, Te_core));

		CD1.Ion_[0].SetCor_cs(ne_core * RM_2_4.cal(ne_core, Te_core));
		CD1.CX_[0].SetCor_cs(ne_core * RM_3_1_4.cal(Tn_core, Te_core));
		CD1.Diss1_[0].SetCor_cs(ne_core * RM_2_14.cal(ne_core, Te_core));
		CD1.Diss2_[0].SetCor_cs(ne_core * RM_2_10_1.cal(ne_core, Te_core));
		CD1.Diss3_[0].SetCor_cs(ne_core * RM_2_10_2.cal(ne_core, Te_core));
		CD1.DS_[1][0].SetCor_cs(ne_core * RM_2_18_1.cal(ne_core, Te_core));
		CD1.DS_[1][1].SetCor_cs(ne_core * RM_2_18_2.cal(ne_core, Te_core));
		CD1.DS_[1][2].SetCor_cs(ne_core * RM_2_22.cal(ne_core, Te_core));
	}

	/// false
	// std::ofstream Ti_out;
	// Ti_out.open("doc/Ti_D_.txt");
	/*for (int i = 1; i < N_poloidal - 1; i++)
	{
		for (int j = 1; j < N_radial - 1; j++)
		{
			if (K_D)
				Ti_D_thermal[i][j] = Ti[i][j] - 0.5 * Dmass * Tools::sqr(ua_D_1[i][j]) / qe;
			if (K_T)
				Ti_T_thermal[i][j] = Ti[i][j] - 0.5 * Tmass * Tools::sqr(ua_T_1[i][j]) / qe;
			if (Ti_D_thermal[i][j] < 0)
				std::cerr << "error of Ti_D_thermal[i][j]" << endl;
			// Ti_out << Ti_D_thermal[i][j] << '\t';
		}
		// Ti_out << endl;
	}*/

	// Ti_out.close();
	/*double V_temp[3], sum_1 = 0, sum_2 = 0;
	for (int i = 0; i < 1000; i++)
	{
		Tools::Vinit_cosine(V_temp, 0, -1);
		std::cout << V_temp[0] << ", " << V_temp[1] << ", " << V_temp[2] << endl;
		if (V_temp[1] < 0)
		{
			sum_1 += 1;
		}
		if (V_temp[2] < 0)
		{
			sum_2 += 1;
		}
	}
	std::cout << sum_1 << ", " << sum_2 << endl;*/

	if (K_H)
	{
		H.CalProb();
		H2.CalProb();
	}
	if (K_D)
	{
		D.CalProb();
		D2.CalProb();
		// std::cout << "D.lambda_min: " << D.lambda_min(0) << " " << D2.lambda_min(0) << endl;
		/*out.open("doc/D2_lambda.txt");
		for (double i = 0; i < 98; i++)
		{
			for (int j = 0; j < N_radial; j++)
			{
				out << D.lambda(i, j, 0) << ' ';
			}
			out << endl;
		}
		out.close();*/
	}
	if (K_T)
	{
		T.CalProb();
		T2.CalProb();
	}
	int xy[2];
	/*for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			xy[0] = i;
			xy[1] = j;
			// std::cout << D2.Ion_[0].cs(xy, 4) << ' ';
			std::cout << D.Ion_[0].cs_[i][j] << ' ';
		}
		std::cout << endl;
	}*/
	/*out.open("doc/D2_Ela_cs_1.txt");
	for (int i = 0; i < N_poloidal; i++)
	{
		for (int j = 0; j < N_radial; j++)
		{
			out << D2.Ela_[0].cs_[i][j] << ' ';
		}
		out << endl;
	}
	out.close();*/
	if (K_Methane)
	{
		CD4.CalProb();
		CD3.CalProb();
		CD2.CalProb();
		CD1.CalProb();
	}
	if (K_C || K_Methane)
	{
		C.CalProb();
	}
	if (K_Ar)
	{
		Ar.CalProb();
	}

	/*// Dump_2D_Global(D.Ion_rate, "Ion_rate_D");
	ofstream out;
	out.open("doc/D2_collision.txt");
	for (double i = 0.1; i < 200; i += 0.1)
	{
		out << i << '\t' << CCD96_H.cal(i, 1e20) << '\t' << CCD96_D.cal(i, 1e20) << '\t' << CCD96_T.cal(i, 1e20) << '\t' << CCD96_H.cal(i, 1e20) * Ratio_DT_Coll << endl;
	}
	out.close();*/

	/// prepare for particle recycle
	for (int i = 0; i < N_radial; i++)
	{
		if (K_Ei == 2)
		{
			Ei_Dion[i] = 2. * Ti[1][i] + 3. * Te[1][i];
			Ei_Dion[i + N_radial] =
				2. * Ti[poloidalLastIndex()][i] +
				3. * Te[poloidalLastIndex()][i];
		}
	}

	/*for (double E_input = 0.1; E_input < 10; E_input = E_input + 0.1)
	{
		std::cout << E_input << '\t';
		std::cout << D_W.E_RefCoeff(1, E_input) << '\t';
		for (double angle_input = 0; angle_input <= 90; angle_input = angle_input + 5)
			std::cout << D_W.E_RefCoeff(K_Reflect, E_input, angle_input) << '\t';
		std::cout << endl;
	}
	for (double E_input = 10; E_input <= 1000; E_input = E_input + 10)
	{
		std::cout << E_input << '\t';
		std::cout << D_W.E_RefCoeff(1, E_input) << '\t';
		for (double angle_input = 0; angle_input <= 90; angle_input = angle_input + 5)
			std::cout << D_W.E_RefCoeff(K_Reflect, E_input, angle_input) << '\t';
		std::cout << endl;
	}*/

	std::vector<double> angle_B_with_target(N_radial * 2);
	if (K_DWTargetActualAngle == 1)
		for (int i = 0; i < N_radial; i++)
		{
			angle_B_with_target[i] = Tools::CalBFieldToWallNormalAngle(B[1][i][0], B[1][i][1], B[1][i][2], Grid4.Cos_Target(i), Grid4.Sin_Target(i));
			angle_B_with_target[i + N_radial] = Tools::CalBFieldToWallNormalAngle(B[poloidalLastIndex()][i][0], B[poloidalLastIndex()][i][1], B[poloidalLastIndex()][i][2], Grid4.Cos_Target(i + N_radial), Grid4.Sin_Target(i + N_radial));
		}
	else
	{
		for (int i = 0; i < N_radial * 2; i++)
		{
			angle_B_with_target[i] = 60;
		}
	}
	DTargetIncidentAngle = angle_B_with_target;
	/// get the area of the target before target source construction
	Grid4.get_sx();
	/*ofstream Out_temp("doc/angel.txt");
	for (int i = 0; i < N_radial * 2; i++)
	{
		Out_temp << i << "\t" << angle_B_with_target[i] << endl;
	}
	Out_temp.close();*/
	// H recycling
	if (K_H)
	{
		if (K_Wallelement == 1)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				if (coeff_recyc_target > H_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]))
				{
					NumPar_H_recyc[i] = NumPar_wall_H[i] * H_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]);
					Tn_H_recyc[i] = H_W.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / H_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) * Ei_Dion[i] / 1.5;
					NumPar_H2_recyc[i] = (coeff_recyc_target - H_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i])) * NumPar_wall_H[i] / 2.;
				}
				else
				{
					NumPar_H_recyc[i] = NumPar_wall_H[i] * coeff_recyc_target;
					Tn_H_recyc[i] = H_W.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / coeff_recyc_target * Ei_Dion[i] / 1.5;
					NumPar_H2_recyc[i] = 0.;
				}
			}
		}
		else if (K_Wallelement == 2)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				if (coeff_recyc_target > H_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]))
				{
					NumPar_H_recyc[i] = NumPar_wall_H[i] * H_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]);
					Tn_H_recyc[i] = H_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / H_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) * Ei_Dion[i] / 1.5;
					NumPar_H2_recyc[i] = (coeff_recyc_target - H_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i])) *
										 NumPar_wall_H[i] / 2.;
				}
				else
				{
					NumPar_H_recyc[i] = NumPar_wall_H[i] * coeff_recyc_target;
					Tn_H_recyc[i] = H_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / coeff_recyc_target * Ei_Dion[i] / 1.5;
					NumPar_H2_recyc[i] = 0.;
				}
			}
		}
		H.CalWeight1(numPar_flight); // numPar_flight is Test flight particle of each grid
		H.RecycledCal(NumPar_H_recyc, Tn_H_recyc);
		H2.RecycledCal(NumPar_H2_recyc);
		const auto h_recycling_flights = AllocateRecyclingFlights(
			H.RecycledSourceSum(), H2.RecycledSourceSum(), numPar_flight_Target);
		H.CalWeight2(NumPar_H_recyc, h_recycling_flights.first);
		H2.CalWeight2(NumPar_H2_recyc, h_recycling_flights.second);
	}
	// D recycling
	if (K_D)
	{
		if (K_DTargetSourceMode == 2)
			ApplyDTargetIonFluxProfile(angle_B_with_target);
		else if (K_DTargetSourceMode != 1)
			throw std::runtime_error("Unsupported K_DTargetSourceMode; use 1 or 2");
		BuildDTargetIncidentAverages();
		// out.open("doc/recycling_D.txt");
		if (K_Wallelement == 1)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				const bool prepared_incident =
					DTargetIncidentModel != 0 && K_DWTrimReflection == 1;
				const double reflection_probability = prepared_incident
					? DTargetFastProbability[i]
					: (K_DWTrimReflection == 1
						   ? (Ei_Dion[i] >= DWTrimERMIN
								  ? D_W_Trim.ReflectionProbability(
										Ei_Dion[i], angle_B_with_target[i])
								  : 0.0)
						   : D_W.n_RefCoeff(
								 K_Reflect, Ei_Dion[i], angle_B_with_target[i]));
				const double fast_probability = prepared_incident
					? reflection_probability
					: std::min(coeff_recyc_target, reflection_probability);
				NumPar_D_recyc[i] = NumPar_wall_D[i] * fast_probability;
				NumPar_D2_recyc[i] =
					(coeff_recyc_target - fast_probability) * NumPar_wall_D[i] / 2.0;
				if (prepared_incident)
				{
					Tn_D_recyc[i] = DTargetMeanReflectedEnergy[i] / 1.5;
				}
				else if (K_DWTrimReflection == 1)
				{
					Tn_D_recyc[i] = D_W_Trim.MeanReflectedEnergy(Ei_Dion[i], angle_B_with_target[i]) / 1.5;
				}
				else if (fast_probability > 0.0)
					Tn_D_recyc[i] =
						D_W.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) /
						fast_probability * Ei_Dion[i] / 1.5;
				else
					Tn_D_recyc[i] = 0.0;
				if (!prepared_incident)
				{
					DTargetMeanIncidentEnergy[i] = Ei_Dion[i];
					DTargetFastProbability[i] = fast_probability;
					DTargetMeanReflectedEnergy[i] = 1.5 * Tn_D_recyc[i];
				}
				// out << NumPar_D_recyc[i] << '\t' << NumPar_D2_recyc[i] << endl;
			}
		}
		else if (K_Wallelement == 2)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				if (coeff_recyc_target > D_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]))
				{
					NumPar_D_recyc[i] = NumPar_wall_D[i] * D_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]);
					Tn_D_recyc[i] = D_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / D_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) * Ei_Dion[i] / 1.5;
					NumPar_D2_recyc[i] = (coeff_recyc_target - D_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i])) * NumPar_wall_D[i] / 2.;
				}
				else
				{
					NumPar_D_recyc[i] = NumPar_wall_D[i] * coeff_recyc_target;
					Tn_D_recyc[i] = D_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / coeff_recyc_target * Ei_Dion[i] / 1.5;
					NumPar_D2_recyc[i] = 0.;
				}
				// out << NumPar_D_recyc[i] << '\t' << NumPar_D2_recyc[i] << endl;
			}
		}
		// out.close();
		D.CalWeight1(numPar_flight); // numPar_flight is Test flight particle of each grid
		D.RecycledCal(NumPar_D_recyc, Tn_D_recyc);
		D2.RecycledCal(NumPar_D2_recyc);
		const auto d_recycling_flights = AllocateRecyclingFlights(
			D.RecycledSourceSum(), D2.RecycledSourceSum(), numPar_flight_Target);
		D.CalWeight2(NumPar_D_recyc, d_recycling_flights.first);
		D2.CalWeight2(NumPar_D2_recyc, d_recycling_flights.second);
		if (K_DBoundarySource)
			D_BoundarySource.Prepare();
	}
	// T recycling
	if (K_T)
	{
		// out.open("doc/recycling_T.txt");
		if (K_Wallelement == 1)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				if (coeff_recyc_target > T_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]))
				{
					NumPar_T_recyc[i] = NumPar_wall_T[i] * T_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]);
					Tn_T_recyc[i] = T_W.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / T_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) * Ei_Dion[i] / 1.5;
					NumPar_T2_recyc[i] = (coeff_recyc_target - T_W.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i])) * NumPar_wall_T[i] / 2.;
				}
				else
				{
					NumPar_T_recyc[i] = NumPar_wall_T[i] * coeff_recyc_target;
					Tn_T_recyc[i] = T_W.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / coeff_recyc_target * Ei_Dion[i] / 1.5;
					NumPar_T2_recyc[i] = 0.;
				}
			}
		}
		else if (K_Wallelement == 2)
		{
			for (int i = 0; i < N_radial * 2; i++)
			{
				if (coeff_recyc_target > T_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]))
				{
					NumPar_T_recyc[i] = NumPar_wall_T[i] * T_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]);
					Tn_T_recyc[i] = T_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / T_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) * Ei_Dion[i] / 1.5;
					NumPar_T2_recyc[i] = (coeff_recyc_target - T_C.n_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i])) * NumPar_wall_T[i] / 2.;
				}
				else
				{
					NumPar_T_recyc[i] = NumPar_wall_T[i] * coeff_recyc_target;
					Tn_T_recyc[i] = T_C.E_RefCoeff(K_Reflect, Ei_Dion[i], angle_B_with_target[i]) / coeff_recyc_target * Ei_Dion[i] / 1.5;
					NumPar_T2_recyc[i] = 0.;
				}
			}
			// out << NumPar_T_recyc[i] << '\t'	<< NumPar_T2_recyc[i] << endl;
		}
		// out.close();
		T.CalWeight1(numPar_flight); // numPar_flight is Test flight particle of each grid
		T.RecycledCal(NumPar_T_recyc, Tn_T_recyc);
		T2.RecycledCal(NumPar_T2_recyc);
		const auto t_recycling_flights = AllocateRecyclingFlights(
			T.RecycledSourceSum(), T2.RecycledSourceSum(), numPar_flight_Target);
		T.CalWeight2(NumPar_T_recyc, t_recycling_flights.first);
		T2.CalWeight2(NumPar_T2_recyc, t_recycling_flights.second);
	}

	for (int i = 0; i <= num_CoreBoundry / 2; i++)
	{
		CoreBoundry[i][0] = Grid[i + 25][1][0];
		CoreBoundry[i][1] = Grid[i + 25][1][4];
		CoreBoundry[num_CoreBoundry - i][0] = Grid[i + 25][1][1];
		CoreBoundry[num_CoreBoundry - i][1] = Grid[i + 25][1][5];
	}

	for (int i = 0; i <= num_GridBoundry / 2; i++)
	{
		GridBoundry[num_GridBoundry - i][0] = Grid[i][radialLastIndex()][3];
		GridBoundry[num_GridBoundry - i][1] = Grid[i][radialLastIndex()][7];
		GridBoundry[i][0] = Grid[poloidalLastIndex() - i][radialLastIndex()][2];
		GridBoundry[i][1] = Grid[poloidalLastIndex() - i][radialLastIndex()][6];
	}

	for (int i = 0; i <= num_PFRBoundry / 2; i++)
	{
		PFRBoundry[num_PFRBoundry - i][0] = Grid[i][radialLastIndex()][3];
		PFRBoundry[num_PFRBoundry - i][1] = Grid[i][radialLastIndex()][7];
		PFRBoundry[i][0] = Grid[poloidalLastIndex() - i][1][1];
		PFRBoundry[i][1] = Grid[poloidalLastIndex() - i][1][6];
	}

	/*out.open("doc/plasmaboundary.txt");
	for (int i = 0; i < Grid1.PLasma_Grid_Boundry_num(); i++)
	{
		out << i << " " << Grid1.PLasma_Grid_Boundry(i, 0) << "\t" << Grid1.PLasma_Grid_Boundry(i, 1) << endl;
	}
	out.close();*/

	/*for (int i = 0; i < N_radial; i++)
	{
		D.SetX_new(0, Grid[poloidalLastIndex()][i][2]);
		D.SetX_new(1, Grid[poloidalLastIndex()][i][6]);
		D.Find();
		std::cout << i << '\t' << D.Zone() << endl;
	}*/
	/*Particle Par_temp;
	Par_temp.SetX_new(0, 2.683969437001213);
	Par_temp.SetX_new(1, 0.363076757001213);
	Par_temp.Find();
	std::cout << Par_temp.Zone() << endl;*/
}

/// @brief Calculate the epx,epy,erx,ery from Magnetic Fields
void get_unit_vector_xy()
{
	int Nt, Nr;
	Nt = N_poloidal;
	Nr = N_radial;
	double d12, d34, d13, d24;
	double dpx, dpy, drx, dry;
	double epx_tmp, epy_tmp, erx_tmp, ery_tmp;
	double eper, erer;

	epx.assign(Nt, std::vector<double>(Nr, 0.0));
	epy.assign(Nt, std::vector<double>(Nr, 0.0));
	erx.assign(Nt, std::vector<double>(Nr, 0.0));
	ery.assign(Nt, std::vector<double>(Nr, 0.0));

	for (int i = 0; i <= Nt - 1; i = i + 1)
	{
		for (int j = 0; j <= Nr - 1; j = j + 1)
		{
			d12 = (Grid[i][j][1] - Grid[i][j][0]) * (Grid[i][j][1] - Grid[i][j][0]) +
				  (Grid[i][j][5] - Grid[i][j][4]) * (Grid[i][j][5] - Grid[i][j][4]);
			d12 = sqrt(d12);

			d34 = (Grid[i][j][2] - Grid[i][j][3]) * (Grid[i][j][2] - Grid[i][j][3]) +
				  (Grid[i][j][6] - Grid[i][j][7]) * (Grid[i][j][6] - Grid[i][j][7]);
			d34 = sqrt(d34);

			d13 = (Grid[i][j][3] - Grid[i][j][0]) * (Grid[i][j][3] - Grid[i][j][0]) +
				  (Grid[i][j][7] - Grid[i][j][4]) * (Grid[i][j][7] - Grid[i][j][4]);
			d13 = sqrt(d13);

			d24 = (Grid[i][j][2] - Grid[i][j][1]) * (Grid[i][j][2] - Grid[i][j][1]) +
				  (Grid[i][j][6] - Grid[i][j][5]) * (Grid[i][j][6] - Grid[i][j][5]);
			d24 = sqrt(d24);

			dpx = (Grid[i][j][1] - Grid[i][j][0]) / d12 +
				  (Grid[i][j][2] - Grid[i][j][3]) / d34;
			dpy = (Grid[i][j][5] - Grid[i][j][4]) / d12 +
				  (Grid[i][j][6] - Grid[i][j][7]) / d34;
			drx = (Grid[i][j][3] - Grid[i][j][0]) / d13 +
				  (Grid[i][j][2] - Grid[i][j][1]) / d24;
			dry = (Grid[i][j][7] - Grid[i][j][4]) / d13 +
				  (Grid[i][j][6] - Grid[i][j][5]) / d24;

			epx_tmp = dpx / sqrt(dpx * dpx + dpy * dpy);
			epy_tmp = dpy / sqrt(dpx * dpx + dpy * dpy);
			erx_tmp = drx / sqrt(drx * drx + dry * dry);
			ery_tmp = dry / sqrt(drx * drx + dry * dry);

			eper = epx_tmp * erx_tmp + epy_tmp * ery_tmp;
			erx_tmp = erx_tmp - eper * epx_tmp;
			ery_tmp = ery_tmp - eper * epy_tmp;

			erer = sqrt(erx_tmp * erx_tmp + ery_tmp * ery_tmp);

			erx_tmp = erx_tmp / erer;
			ery_tmp = ery_tmp / erer;

			epx[i][j] = epx_tmp;
			epy[i][j] = epy_tmp;
			erx[i][j] = erx_tmp;
			ery[i][j] = ery_tmp;
		}
	}
	for (int i = 0; i < N_poloidal; i++)
		for (int j = 0; j < N_radial; j++)
		{
			if (i == 0 || i == 97 || j == 0 || j == 37)
			{
				int x_temp, y_temp;
				x_temp = (i == 0) * (1) + (i == 97) * (-1) + i;
				y_temp = (j == 0) * (1) + (j == 37) * (-1) + j;
				epx[i][j] = epx[x_temp][y_temp];
				epy[i][j] = epy[x_temp][y_temp];
				erx[i][j] = erx[x_temp][y_temp];
				ery[i][j] = ery[x_temp][y_temp];
			}
		}
	/*ofstream f1, f2, f3, f4;
	f1.open(Outputpath + "epx.txt");
	f2.open(Outputpath + "epy.txt");
	f3.open(Outputpath + "erx.txt");
	f4.open(Outputpath + "ery.txt");
	for (int i = 0; i <= Nt - 1; i = i + 1)
	{
			for (int j = 0; j <= Nr - 1; j = j + 1)
			{
					f1 << epx[i][j] << '\t';
					f2 << epy[i][j] << '\t';
					f3 << erx[i][j] << '\t';
					f4 << ery[i][j] << '\t';
			}
			f1 << endl;
			f2 << endl;
			f3 << endl;
			f4 << endl;
	}
	f1.close();
	f2.close();
	f3.close();
	f4.close();*/
}

void get_center()
{
	int Nt, Nr;
	Nt = N_poloidal;
	Nr = N_radial;
	crx.assign(Nt, std::vector<double>(Nr, 0.0));
	cry.assign(Nt, std::vector<double>(Nr, 0.0));

	for (int i = 0; i <= Nt - 1; i = i + 1)
	{
		for (int j = 0; j <= Nr - 1; j = j + 1)
		{
			Tools::get_line_intersection(Grid[i][j][0], Grid[i][j][4], Grid[i][j][2],
										 Grid[i][j][6], Grid[i][j][1], Grid[i][j][5],
										 Grid[i][j][3], Grid[i][j][7], &crx[i][j],
										 &cry[i][j]);
		}
	}
	/*ofstream f1, f2;
	f1.open(Outputpath + "crx.txt");
	f2.open(Outputpath + "cry.txt");
	for (int i = 0; i <= Nt - 1; i = i + 1)
	{
			for (int j = 0; j <= Nr - 1; j = j + 1)
			{
					f1 << crx[i][j] << '\t';
					f2 << cry[i][j] << '\t';
			}
			f1 << endl;
			f2 << endl;
	}
	f1.close();
	f2.close();*/
}
void get_dTi()
{
	int Nt, Nr, Np;
	Nt = N_poloidal;
	Nr = N_radial;
	Np = 25;

	std::vector<std::vector<double>> TIx, TIy;
	TIx.assign(Nt, std::vector<double>(Nr, 0.0));
	TIy.assign(Nt, std::vector<double>(Nr, 0.0));
	dTip.assign(Nt, std::vector<double>(Nr, 0.0));
	dTir.assign(Nt, std::vector<double>(Nr, 0.0));

	dTix.assign(Nt, std::vector<double>(Nr, 0.0));
	dTiy.assign(Nt, std::vector<double>(Nr, 0.0));
	dTiz.assign(Nt, std::vector<double>(Nr, 0.0));

	int i, j;

	for (i = 0; i <= Nt - 1; i = i + 1)
	{
		for (j = 0; j <= Nr - 1; j = j + 1)
		{
			dTiz[i][j] = 0;
		}
	}

	double x1, x2, x3, x4, y1, y2, y3, y4, d1, d2, d3, d4;
	double dt1, dt2, dt3, dt4;

	for (i = 1; i <= Nt - 2; i = i + 1)
	{
		for (j = 1; j <= Nr - 2; j = j + 1)
		{
			x1 = crx[i + 1][j] - crx[i][j];
			x2 = crx[i][j + 1] - crx[i][j];
			x3 = crx[i - 1][j] - crx[i][j];
			x4 = crx[i][j - 1] - crx[i][j];

			y1 = cry[i + 1][j] - cry[i][j];
			y2 = cry[i][j + 1] - cry[i][j];
			y3 = cry[i - 1][j] - cry[i][j];
			y4 = cry[i][j - 1] - cry[i][j];

			d1 = sqrt(x1 * x1 + y1 * y1);
			d2 = sqrt(x2 * x2 + y2 * y2);
			d3 = sqrt(x3 * x3 + y3 * y3);
			d4 = sqrt(x4 * x4 + y4 * y4);

			dt1 = Ti[i + 1][j] - Ti[i][j];
			dt2 = Ti[i][j + 1] - Ti[i][j];
			dt3 = -Ti[i - 1][j] + Ti[i][j];
			dt4 = -Ti[i][j - 1] + Ti[i][j];

			TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
			TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
		}
	}

	i = 0;
	for (j = 1; j <= Nr - 2; j = j + 1)
	{
		x1 = crx[i + 1][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		// x3 = crx[i - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[i + 1][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		// y3 = cry[i - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		// d3 = sqrt(x3*x3 + y3*y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		dt1 = Ti[i + 1][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		// dt3 = -Ti[i-1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = dt1 / d1;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	i = Nt - 1;
	for (j = 1; j <= Nr - 2; j = j + 1)
	{
		// x1 = crx[i + 1][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[i - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		// y1 = cry[i + 1][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[i - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		// d1 = sqrt(x1*x1 + y1*y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		// dt1 = Ti[i+1][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[i - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = dt3 / d3;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	j = 0;
	for (i = 1; i <= Nt - 2; i = i + 1)
	{
		x1 = crx[i + 1][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[i - 1][j] - crx[i][j];
		// x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[i + 1][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[i - 1][j] - cry[i][j];
		// y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		// d4 = sqrt(x4*x4 + y4*y4);

		dt1 = Ti[i + 1][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[i - 1][j] + Ti[i][j];
		// dt4 = -Ti[i][j-1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = dt2 / d2;
	}

	j = Nr - 1;
	for (i = 1; i <= Nt - 2; i = i + 1)
	{
		x1 = crx[i + 1][j] - crx[i][j];
		// x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[i - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[i + 1][j] - cry[i][j];
		// y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[i - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		// d4 = sqrt(x4*x4 + y4*y4);

		dt1 = Ti[i + 1][j] - Ti[i][j];
		// dt2 = Ti[i][j+1] - Ti[i][j];
		dt3 = -Ti[i - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = dt4 / d4;
	}

	i = 0;
	j = 0;
	x1 = crx[i + 1][j] - crx[i][j];
	x2 = crx[i][j + 1] - crx[i][j];

	y1 = cry[i + 1][j] - cry[i][j];
	y2 = cry[i][j + 1] - cry[i][j];

	d1 = sqrt(x1 * x1 + y1 * y1);
	d2 = sqrt(x2 * x2 + y2 * y2);

	dt1 = Ti[i + 1][j] - Ti[i][j];
	dt2 = Ti[i][j + 1] - Ti[i][j];

	TIx[i][j] = dt1 / d1;
	TIy[i][j] = dt2 / d2;

	i = 0;
	j = Nr - 1;
	x1 = crx[i + 1][j] - crx[i][j];
	x4 = crx[i][j - 1] - crx[i][j];

	y1 = cry[i + 1][j] - cry[i][j];
	y4 = cry[i][j - 1] - cry[i][j];

	d1 = sqrt(x1 * x1 + y1 * y1);
	d4 = sqrt(x4 * x4 + y4 * y4);

	dt1 = Ti[i + 1][j] - Ti[i][j];
	dt4 = -Ti[i][j - 1] + Ti[i][j];

	TIx[i][j] = dt1 / d1;
	TIy[i][j] = dt4 / d4;

	i = Nt - 1;
	j = 0;
	x2 = crx[i][j + 1] - crx[i][j];
	x3 = crx[i - 1][j] - crx[i][j];

	y2 = cry[i][j + 1] - cry[i][j];
	y3 = cry[i - 1][j] - cry[i][j];

	d2 = sqrt(x2 * x2 + y2 * y2);
	d3 = sqrt(x3 * x3 + y3 * y3);

	dt2 = Ti[i][j + 1] - Ti[i][j];
	dt3 = -Ti[i - 1][j] + Ti[i][j];

	TIx[i][j] = dt3 / d3;
	TIy[i][j] = dt2 / d2;

	i = Nt - 1;
	j = Nr - 1;

	x3 = crx[i - 1][j] - crx[i][j];
	x4 = crx[i][j - 1] - crx[i][j];

	y3 = cry[i - 1][j] - cry[i][j];
	y4 = cry[i][j - 1] - cry[i][j];

	d3 = sqrt(x3 * x3 + y3 * y3);
	d4 = sqrt(x4 * x4 + y4 * y4);

	dt3 = -Ti[i - 1][j] + Ti[i][j];
	dt4 = -Ti[i][j - 1] + Ti[i][j];

	TIx[i][j] = dt3 / d3;
	TIy[i][j] = dt4 / d4;

	// i=25,74(c里的 Np-1  Nt-Np)处 私有区接到了芯部
	i = Np - 1; // inner private
	int N_sp = 19;
	for (j = 1; j <= N_sp; j = j + 1)
	{
		x1 = crx[Nt - Np][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[i - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[Nt - Np][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[i - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		dt1 = Ti[Nt - Np][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[i - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	i = Np; // inner core
	for (j = 1; j <= N_sp; j = j + 1)
	{
		x1 = crx[i + 1][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[Nt - Np - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[i + 1][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[Nt - Np - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		dt1 = Ti[i + 1][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[Nt - Np - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	i = Nt - Np - 1; // outer core
	for (j = 1; j <= N_sp; j = j + 1)
	{
		x1 = crx[Np][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[i - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[Np][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[i - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		dt1 = Ti[Np][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[i - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	i = Nt - Np; // outer private
	for (j = 1; j <= N_sp; j = j + 1)
	{
		x1 = crx[i + 1][j] - crx[i][j];
		x2 = crx[i][j + 1] - crx[i][j];
		x3 = crx[Np - 1][j] - crx[i][j];
		x4 = crx[i][j - 1] - crx[i][j];

		y1 = cry[i + 1][j] - cry[i][j];
		y2 = cry[i][j + 1] - cry[i][j];
		y3 = cry[Np - 1][j] - cry[i][j];
		y4 = cry[i][j - 1] - cry[i][j];

		d1 = sqrt(x1 * x1 + y1 * y1);
		d2 = sqrt(x2 * x2 + y2 * y2);
		d3 = sqrt(x3 * x3 + y3 * y3);
		d4 = sqrt(x4 * x4 + y4 * y4);

		dt1 = Ti[i + 1][j] - Ti[i][j];
		dt2 = Ti[i][j + 1] - Ti[i][j];
		dt3 = -Ti[Np - 1][j] + Ti[i][j];
		dt4 = -Ti[i][j - 1] + Ti[i][j];

		TIx[i][j] = (dt1 / d1 + dt3 / d3) / 2.0;
		TIy[i][j] = (dt2 / d2 + dt4 / d4) / 2.0;
	}

	for (i = 0; i <= Nt - 1; i = i + 1)
	{
		for (j = 0; j <= Nr - 1; j = j + 1)
		{
			dTip[i][j] = TIx[i][j];
			dTir[i][j] = TIy[i][j];
			// TIy[i][j] = 0;
			// TIx[i][j] = 0;
			dTix[i][j] = TIx[i][j] * epx[i][j] + TIy[i][j] * erx[i][j];
			dTiy[i][j] = TIx[i][j] * epy[i][j] + TIy[i][j] * ery[i][j];
		}
	}

	// for revised 2

	// for (j = 0; j <= Nr - 1; j = j + 1)
	//{
	//	dTip[1][j] = dTip[2][j];
	// }

	/*ofstream f1, f2;
	f1.open(Outputpath + "dTix.txt");
	f2.open(Outputpath + "dTiy.txt");

	for (i = 0; i <= Nt - 1; i = i + 1)
	{
			for (j = 0; j <= Nr - 1; j = j + 1)
			{
					f1 << dTix[i][j] << '\t';
					f2 << dTiy[i][j] << '\t';
			}
			f1 << endl;
			f2 << endl;
	}
	f1.close();
	f2.close();

	ofstream f3, f4;
	f3.open(Outputpath + "dTip.txt");
	f4.open(Outputpath + "dTir.txt");

	for (i = 0; i <= Nt - 1; i = i + 1)
	{
			for (j = 0; j <= Nr - 1; j = j + 1)
			{
					f3 << dTip[i][j] << '\t';
					f4 << dTir[i][j] << '\t';
			}
			f3 << endl;
			f4 << endl;
	}
	f3.close();
	f4.close();*/
}

void CX_DT_Fix()
{
	/// Output the collision cross-section
	double k1[301], k2[301], ave1, ave2, ave3, ave4, sum1, sum2, Count, variance1, variance2, variance3, variance4, K1 = 3, K2 = 3;
	variance1 = 0;
	variance2 = 0;
	ave1 = 0;
	ave2 = 0;
	for (double j = 3; j < 10; j += 0.01)
	{
		Ratio_D_Coll = pow(Dmass / Hmass, -1. / j);
		Ratio_T_Coll = pow(Tmass / Hmass, -1. / j);
		Count = 0;
		sum1 = 0;
		sum2 = 0;
		for (double i = 0.1; i < 100; i *= 1.1)
		{
			k1[(int)(i * 10)] = Tools::cds(log10(CCD96_H.cal(i, 1e20) * Ratio_D_Coll) / log10(CCD96_D.cal(i, 1e20)));
			ave1 += k1[(int)(i * 10)];
			k2[(int)(i * 10)] = Tools::cds(log10(CCD96_H.cal(i, 1e20) * Ratio_T_Coll) / log10(CCD96_T.cal(i, 1e20)));
			ave2 += k2[(int)(i * 10)];
			Count += 1;
			// std::cout << i << '\t' << CCD96_H.cal(i, 1e20) << '\t' << CCD96_D.cal(i, 1e20) << '\t'
			//		  << CCD96_T.cal(i, 1e20) << '\t' << CCD96_H.cal(i, 1e20) * Ratio_D_Coll << '\t'
			//		  << CCD96_H.cal(i, 1e20) * Ratio_T_Coll << '\t' << CCD96_H.cal(i, 1e20) * Ratio_DT_Coll << endl;
		}
		ave1 /= Count;
		ave2 /= Count;
		for (int i = 0; i < (int)Count; i++)
		{
			sum1 += Tools::sqr(k1[i] - ave1);
			sum2 += Tools::sqr(k2[i] - ave2);
		}
		variance3 = sqrt(sum1);
		variance4 = sqrt(sum2);

		/*if (j != 3)
		{
			if (variance1 > variance3)
				K1 = j;

			if (variance2 > variance4)
				K2 = j;
		}*/
		if (j != 3)
		{
			if (ave1 < ave3)
				K1 = j;

			if (ave2 < ave4)
				K2 = j;
		}
		variance1 = variance3;
		variance2 = variance4;
		ave3 = ave1;
		ave4 = ave2;
		// std::cout << j << '\t' << variance1 << '\t' << variance2 << endl
		//		  << K1 << '\t' << K2 << endl;
	}
	Ratio_D_Coll = pow(Dmass / Hmass, -1. / K1);
	Ratio_T_Coll = pow(Tmass / Hmass, -1. / K2);
	Ratio_DT_Coll = pow(2 * Dmass * Tmass / (Dmass + Tmass) / Hmass, -1. / ((K1 + K2) / 2));
	// std::cout << K1 << '\t' << K2 << '\t' << Ratio_DT_Coll << endl;
}
