#include "iostream"
#include <ctime>
#include <cstdlib>
#include "sys/timeb.h"
#include <random>

#define qe 1.602e-19
using namespace std;

double Random()
{
/*double p;
// srand((unsigned)time(NULL));
do
{
	p = (double)(rand() / (double)RAND_MAX);
} while (p <= 1.0E-30);
// p = (double)(rand() / (double)RAND_MAX);
return p;*/
#define N 99999999
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_int_distribution<int> distr(1, N);
	// srand((unsigned)time(NULL));
	return (double)distr(eng) / (double)(N + 1);

	/*int min = 1, max = 9999;
	random_device seed;
	ranlux48 engine(seed());
	uniform_int_distribution<> distrib(min, max);
	int random = distrib(engine);
	return (double)random / (double)(max + 1);*/
}

double Maxwell(double E, double Mass)
{
	double RR, FuncA, FuncB, Vel;
	do
	{
		RR = Random();
		FuncA = RR * RR;
		RR = Random();
		FuncB = -2.71828 * RR * log(RR);
	} while (FuncA > FuncB);
	Vel = sqrt(-3. * qe * E * log(RR) / Mass);
	return Vel;
}

double intersect(double *A, double *B, int i)
{
	if (i == 1)
		return A[1] * B[2] - A[2] * B[1];
	else if (i == 2)
		return A[2] * B[0] - A[0] * B[2];
	else if (i == 3)
		return A[0] * B[1] - A[1] * B[0];
	else
		throw std::invalid_argument("intersect(): i must be 1, 2, or 3");
}