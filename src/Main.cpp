#include "Global.h"

#include <chrono>

namespace
{
using Clock = std::chrono::steady_clock;

double elapsedSeconds(const Clock::time_point &start, const Clock::time_point &end)
{
	return std::chrono::duration<double>(end - start).count();
}

void printRuntimeReport(const std::vector<std::pair<std::string, double>> &stages, double total_seconds)
{
	std::cout << std::endl;
	std::cout << "===== Runtime report =====" << std::endl;
	for (const auto &stage : stages)
	{
		std::cout << std::left << std::setw(12) << stage.first << std::right << std::fixed << std::setprecision(3)
		          << stage.second << " s" << std::endl;
	}
	std::cout << std::left << std::setw(12) << "Total" << std::right << std::fixed << std::setprecision(3)
	          << total_seconds << " s" << std::endl;
	std::cout << "==========================" << std::endl;

	if (!Outputpath.empty())
	{
		std::ofstream out(Outputpath + "runtime_report.txt");
		if (out.is_open())
		{
			out << "Stage Seconds" << std::endl;
			for (const auto &stage : stages)
			{
				out << stage.first << " " << std::fixed << std::setprecision(6) << stage.second << std::endl;
			}
			out << "Total " << std::fixed << std::setprecision(6) << total_seconds << std::endl;
		}
	}
}

void printVarianceReductionReport()
{
	std::cout << std::endl;
	std::cout << "===== Variance reduction report =====" << std::endl;
	std::cout << "Russian roulette: " << (K_Roulette ? "on" : "off")
			  << ", trials=" << RouletteTrials << ", survived=" << RouletteSurvived
			  << ", killed=" << RouletteKilled << std::endl;
	std::cout << "Regional importance sampling: " << (K_Splitting ? "on" : "off")
			  << ", events=" << SplitEvents << ", children=" << SplitChildrenCreated
			  << ", max_pending=" << SplitMaxPendingStates << std::endl;
	std::cout << "Importance regions: outside=" << RegionImportance[0]
			  << ", divertor=" << RegionImportance[1]
			  << ", main_chamber=" << RegionImportance[2]
			  << ", main_poloidal=[" << ImportanceMainPoloidalBegin
			  << ',' << ImportanceMainPoloidalEnd << ']' << std::endl;
	std::cout << "Importance roulette: trials=" << ImportanceRouletteTrials
			  << ", survived=" << ImportanceRouletteSurvived
			  << ", killed=" << ImportanceRouletteKilled
			  << ", split_suppressed_min=" << SplitSuppressedByMinWeight << std::endl;
	std::cout << "=====================================" << std::endl;

	if (!Outputpath.empty())
	{
		std::ofstream out(Outputpath + "variance_reduction_report.txt");
		if (out.is_open())
		{
			out << "Parameter Value\n"
				<< "K_Roulette " << K_Roulette << '\n'
				<< "K_Splitting " << K_Splitting << '\n'
				<< "W_RouletteMin " << W_RouletteMin << '\n'
				<< "W_RouletteTarget " << W_RouletteTarget << '\n'
				<< "W_SplitMax " << W_SplitMax << '\n'
				<< "W_SplitTarget " << W_SplitTarget << '\n'
				<< "W_SplitMin " << W_SplitMin << '\n'
				<< "MaxSplit " << MaxSplit << '\n'
				<< "MaxSplitDepth " << MaxSplitDepth << '\n'
				<< "ImportanceOutside " << RegionImportance[0] << '\n'
				<< "ImportanceDivertor " << RegionImportance[1] << '\n'
				<< "ImportanceMainChamber " << RegionImportance[2] << '\n'
				<< "ImportanceMainPoloidalBegin " << ImportanceMainPoloidalBegin << '\n'
				<< "ImportanceMainPoloidalEnd " << ImportanceMainPoloidalEnd << '\n'
				<< "RouletteTrials " << RouletteTrials << '\n'
				<< "RouletteSurvived " << RouletteSurvived << '\n'
				<< "RouletteKilled " << RouletteKilled << '\n'
				<< "SplitEvents " << SplitEvents << '\n'
				<< "SplitChildrenCreated " << SplitChildrenCreated << '\n'
				<< "SplitMaxPendingStates " << SplitMaxPendingStates << '\n'
				<< "SplitSuppressedByMinWeight " << SplitSuppressedByMinWeight << '\n'
				<< "ImportanceRouletteTrials " << ImportanceRouletteTrials << '\n'
				<< "ImportanceRouletteSurvived " << ImportanceRouletteSurvived << '\n'
				<< "ImportanceRouletteKilled " << ImportanceRouletteKilled << '\n';
		}
	}
}
}

int main(int argc, char *argv[])

{
	std::vector<std::pair<std::string, double>> runtime_stages;
	auto total_start = Clock::now();
	auto stage_start = total_start;

	Initialize(argc, argv);
	auto stage_end = Clock::now();
	runtime_stages.push_back({"Initialize", elapsedSeconds(stage_start, stage_end)});
	stage_start = stage_end;

	Read();
	stage_end = Clock::now();
	runtime_stages.push_back({"Read", elapsedSeconds(stage_start, stage_end)});
	stage_start = stage_end;

	Prepare();
	stage_end = Clock::now();
	runtime_stages.push_back({"Prepare", elapsedSeconds(stage_start, stage_end)});
	stage_start = stage_end;

	Moncar();
	stage_end = Clock::now();
	runtime_stages.push_back({"Moncar", elapsedSeconds(stage_start, stage_end)});
	stage_start = stage_end;

	Output();
	stage_end = Clock::now();
	runtime_stages.push_back({"Output", elapsedSeconds(stage_start, stage_end)});

	printRuntimeReport(runtime_stages, elapsedSeconds(total_start, stage_end));
	printVarianceReductionReport();
	return 0;
}
// CalPec();
/*ofstream oup("doc/cross_CX_H_1e18.txt");
	double n___ = 1e18;
	for (double i = 0.1; i < 100; i += 0.1)
	{
		oup << i << "\t" << CCD96_H.cal(i, n___) << endl;
		// oup << i << "\t" << R3_2_3r_H4.cal(n___, i) << "\t" << R3_2_3d_H4.cal(n___, i) << "\t" << R3_2_3i_H4.cal(n___, i) << "\t" << R2_2_17r_H4.cal(n___, i) << "\t" << R2_2_17d_H4.cal(n___, i) << endl;
	}
	oup.close();

	ofstream oup("doc/cross-section_CX_n-n.txt");
	double n___ = 1e18;
	for (double i = 0.1; i < 100; i += 0.1)
	{
		oup << i << "\t" << R_H_H.cal(n___, i) << "\t" << R_H_H2.cal(n___, i) << "\t" << R_H2_H.cal(n___, i) << "\t" << R_H2_H2.cal(n___, i) << endl;
		// oup << i << "\t" << R3_2_3r_H4.cal(n___, i) << "\t" << R3_2_3d_H4.cal(n___, i) << "\t" << R3_2_3i_H4.cal(n___, i) << "\t" << R2_2_17r_H4.cal(n___, i) << "\t" << R2_2_17d_H4.cal(n___, i) << endl;
	}
	oup.close();
	ofstream oup("doc/cross_Mar_H2.txt");
	for (double i = 0.01; i < 10; i += 0.01)
	{
		oup << i << "\t" << R3_2_3r_H4.cal(1e16, i) << "\t" << R3_2_3r_H4.cal(5e16, i) << "\t" << R3_2_3r_H4.cal(1e17, i) << "\t" << R3_2_3r_H4.cal(5e17, i) << "\t" << R3_2_3r_H4.cal(1e18, i) << "\t" << R3_2_3r_H4.cal(5e18, i) << endl;
		// oup << i << "\t" << R3_2_3r_H4.cal(n___, i) << "\t" << R3_2_3d_H4.cal(n___, i) << "\t" << R3_2_3i_H4.cal(n___, i) << "\t" << R2_2_17r_H4.cal(n___, i) << "\t" << R2_2_17d_H4.cal(n___, i) << endl;
	}
	oup.close();*/
/*
	std::vector<double> V_temp(3);
	ofstream oup("doc/cos_test.txt");
	for (int i = 0; i < 100; i++)
	{
		V_temp[0] = 0;
		V_temp[1] = 11;
		V_temp[2] = 1;
		Tools::calculateReflectionVelocity(V_temp, 0, -1, 0);
		std::cout << i << "\t" << V_temp[0] << "\t" << V_temp[1] << "\t" << V_temp[2] << "\t" << V_temp[0] * V_temp[0] + V_temp[1] * V_temp[1] + V_temp[2] * V_temp[2] << endl;
	}
	oup.close();
 */
