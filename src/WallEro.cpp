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
    // std::cout << "jianli: " << n << endl;
    CollSource_ = 0;
    num_wall_ = n;
    N_p_ = N_p;
    N_r_ = N_r;
    // 1D 汇总
    All_PartoWall.assign(num_wall_, 0.0);
    All_E_PartoWall.assign(num_wall_, 0.0);

    // 工具：生成 3D [wall][Np][Nr] 并置零
    auto make3 = [&](std::vector<std::vector<std::vector<double>>> &a)
    {
        a.assign(num_wall_, std::vector<std::vector<double>>(
                                N_p_, std::vector<double>(N_r_, 0.0)));
    };
    make3(sum_PartoWall);
    make3(sum_E_PartoWall);
    make3(PartoWall);
    make3(E_PartoWall);
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
    int i = xy[0], j = xy[1];
    if (i < 0 || i >= N_p_ || j < 0 || j >= N_r_)
        return;

    sum_PartoWall[num_Ero_wall][i][j] += n;
    sum_E_PartoWall[num_Ero_wall][i][j] += n * T;
    PartoWall[num_Ero_wall][i][j] += n;
    E_PartoWall[num_Ero_wall][i][j] += n * T;

    All_PartoWall[num_Ero_wall] += n;
    All_E_PartoWall[num_Ero_wall] += n * T;
}

void WallEro::stat()
{
    for (int i = 0; i < num_wall_; i++)
    {
        for (int j = 0; j < N_p_; j++)
        {
            for (int k = 0; k < N_r_; k++)
            {
                PartoWall[i][j][k] = sum_PartoWall[i][j][k];
                if (sum_PartoWall[i][j][k] != 0)
                    E_PartoWall[i][j][k] = sum_E_PartoWall[i][j][k] / sum_PartoWall[i][j][k];
            }
        }
        if (All_PartoWall[i] != 0)
            All_E_PartoWall[i] /= All_PartoWall[i];
    }
}