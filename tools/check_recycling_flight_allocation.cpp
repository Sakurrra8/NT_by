#include "RecyclingFlightAllocation.h"

#include <cstdlib>
#include <iostream>

namespace
{
void Check(double atomic_source, double molecular_source, int total_flights,
	int expected_atomic, int expected_molecular)
{
	const auto allocation = AllocateRecyclingFlights(
		atomic_source, molecular_source, total_flights);
	if (allocation.first != expected_atomic ||
		allocation.second != expected_molecular)
	{
		std::cerr << "allocation mismatch: got " << allocation.first << ", "
				  << allocation.second << "; expected " << expected_atomic << ", "
				  << expected_molecular << '\n';
		std::exit(EXIT_FAILURE);
	}
}
}

int main()
{
	Check(0.0, 0.0, 50000, 0, 0);
	Check(100.0, 0.0, 100, 100, 0);
	Check(0.0, 50.0, 100, 0, 100);
	Check(100.0, 50.0, 100, 50, 50);

	// Latest 2 MW, 5e19 target source: D atoms/s and D2 molecules/s.
	Check(1.347521728e23, 1.594935830e22, 50000, 40429, 9571);

	std::cout << "recycling flight allocation checks passed\n";
	return EXIT_SUCCESS;
}
