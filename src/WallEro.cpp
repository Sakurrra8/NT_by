#include "WallEro.h"

int WallEro::CollSource()
{
    return CollSource_;
}

int WallEro::num_wall()
{
    return num_wall_;
}

int WallEro::N_p()
{
    return N_p_;
}
int WallEro::N_r()
{
    return N_r_;
}

void WallEro::setCollSource(int CollSource)
{
    CollSource_ = CollSource;
}

/// @brief Initialization function of class WallEro
/// @param n: 1 for Wall; 2 for Core boundry; 3 for charged particle on target
void WallEro::WallEroInit(int n, int N_p, int N_r)
{
    CollSource_ = 0;
    num_wall_ = n;
    N_p_ = N_p;
    N_r_ = N_r;

    All_PartoWall.assign(num_wall_, 0.0);
    All_E_PartoWall.assign(num_wall_, 0.0);
    const size_t cells = static_cast<size_t>(num_wall_) * N_p_ * N_r_;
    particleCount_.assign(cells, 0.0);
    energySum_.assign(cells, 0.0);
}

size_t WallEro::index(int wall, int i, int j) const
{
    return (static_cast<size_t>(wall) * N_p_ + i) * N_r_ + j;
}

WallEro::WallEro()
{
}

WallEro::~WallEro()
{
}

void WallEro::AddPar(int num_Ero_wall, double n, double T, int source, int xy[2])
{
    if (num_Ero_wall < 0 || num_Ero_wall >= num_wall_)
        return;

    All_PartoWall[num_Ero_wall] += n;
    All_E_PartoWall[num_Ero_wall] += n * T;

    int i = xy[0], j = xy[1];
    if (i < 0 || i >= N_p_ || j < 0 || j >= N_r_)
        return;

    const size_t cell = index(num_Ero_wall, i, j);
    particleCount_[cell] += n;
    energySum_[cell] += n * T;

}

double WallEro::N(int wall, int i, int j) const
{
    return particleCount_[index(wall, i, j)];
}

double WallEro::T(int wall, int i, int j) const
{
    const size_t cell = index(wall, i, j);
    const double count = particleCount_[cell];
    return count == 0.0 ? 0.0 : energySum_[cell] / count;
}

double WallEro::T_all(int wall) const
{
    const double count = All_PartoWall[wall];
    return count == 0.0 ? 0.0 : All_E_PartoWall[wall] / count;
}

void WallEro::stat()
{
    // Averages are computed on demand by T() and T_all().
}
