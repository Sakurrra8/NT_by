#include "Particle.h"

#include <cmath>
#include <iostream>
#include <queue>

namespace
{
bool nearlyEqual(double a, double b, double tolerance)
{
	return std::abs(a - b) <= tolerance;
}
}

int main()
{
	N_poloidal = 98;
	N_radial = 38;
	K_Splitting = 1;
	MaxSplit = 4;
	W_SplitMin = 0.05;
	ImportanceMainPoloidalBegin = 25;
	ImportanceMainPoloidalEnd = 72;
	RegionImportance = {{1.0, 2.0, 1.0}};

	Particle split_particle;
	split_particle.SetWeight(10.0);
	split_particle.SetXY(0, 40);
	split_particle.SetXY(1, 10);
	std::queue<Particle::State> pending;
	split_particle.ApplyRegionalImportance(pending);
	if (split_particle.PhysicalImportanceRegion() != 2 || !pending.empty())
	{
		std::cerr << "main-chamber region classification failed\n";
		return 1;
	}

	split_particle.SetXY(0, 10);
	split_particle.ApplyRegionalImportance(pending);
	if (split_particle.PhysicalImportanceRegion() != 1 ||
		pending.size() != 1 || !nearlyEqual(split_particle.Weight(), 5.0, 1e-12))
	{
		std::cerr << "divertor importance split failed\n";
		return 1;
	}
	if (!nearlyEqual(
			split_particle.Weight() * static_cast<double>(pending.size() + 1),
			10.0, 1e-12))
	{
		std::cerr << "regional split weight conservation failed\n";
		return 1;
	}

	Particle outside_particle;
	outside_particle.SetXY(0, -1);
	outside_particle.SetXY(1, -1);
	if (outside_particle.PhysicalImportanceRegion() != 0)
	{
		std::cerr << "outside region classification failed\n";
		return 1;
	}

	Particle min_weight_particle;
	min_weight_particle.SetWeight(0.05);
	min_weight_particle.SetXY(0, 40);
	min_weight_particle.SetXY(1, 10);
	std::queue<Particle::State> min_pending;
	min_weight_particle.ApplyRegionalImportance(min_pending);
	min_weight_particle.SetXY(0, 10);
	min_weight_particle.ApplyRegionalImportance(min_pending);
	if (!min_pending.empty() ||
		!nearlyEqual(min_weight_particle.Weight(), 0.05, 1e-12) ||
		SplitSuppressedByMinWeight == 0)
	{
		std::cerr << "regional minimum child weight guard failed\n";
		return 1;
	}

	constexpr int samples = 100000;
	double final_weight_sum = 0.0;
	for (int i = 0; i < samples; ++i)
	{
		Particle roulette_particle;
		roulette_particle.SetWeight(1.0);
		roulette_particle.SetXY(0, 10);
		roulette_particle.SetXY(1, 10);
		std::queue<Particle::State> roulette_pending;
		roulette_particle.ApplyRegionalImportance(roulette_pending);
		roulette_particle.SetXY(0, 40);
		roulette_particle.ApplyRegionalImportance(roulette_pending);
		final_weight_sum += roulette_particle.Weight();
	}
	const double average_weight = final_weight_sum / samples;
	if (!nearlyEqual(average_weight, 1.0, 0.02))
	{
		std::cerr << "regional roulette expectation check failed: "
				  << average_weight << '\n';
		return 1;
	}

	std::cout << "regional importance tests passed\n";
	return 0;
}
