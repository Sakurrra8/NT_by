#ifndef RECYCLING_FLIGHT_ALLOCATION_H
#define RECYCLING_FLIGHT_ALLOCATION_H

#include <algorithm>
#include <cmath>
#include <utility>

inline std::pair<int, int> AllocateRecyclingFlights(
	double atomic_source, double molecular_source, int total_flights)
{
	atomic_source = std::max(0.0, atomic_source);
	molecular_source = std::max(0.0, molecular_source);
	if (total_flights <= 0 || atomic_source + molecular_source <= 0.0)
		return {0, 0};
	if (atomic_source <= 0.0)
		return {0, total_flights};
	if (molecular_source <= 0.0)
		return {total_flights, 0};

	int atomic_flights = static_cast<int>(std::llround(
		total_flights * atomic_source / (atomic_source + molecular_source)));
	atomic_flights = std::max(1, std::min(total_flights - 1, atomic_flights));
	return {atomic_flights, total_flights - atomic_flights};
}

#endif
