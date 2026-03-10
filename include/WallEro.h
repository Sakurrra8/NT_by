#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
// #include <string.h>
#include <string>
#include <cassert>
#include "Reflect.h"
#include "EIRENE.h"
#include <malloc.h>
#include <unistd.h>
#include <iomanip>
#include "ADAS.h"
#include <math.h>

class WallEro
{
private:
    int num_wall_;
    int N_p_;
    int N_r_;
    int CollSource_;
    std::vector<std::vector<std::vector<double>>> sum_PartoWall;
    std::vector<std::vector<std::vector<double>>> sum_E_PartoWall;

public:
    WallEro();
    ~WallEro();
    vector<double> All_PartoWall;
    vector<double> All_E_PartoWall;
    std::vector<std::vector<std::vector<double>>> PartoWall;
    std::vector<std::vector<std::vector<double>>> E_PartoWall;

    int num_wall();
    int N_p();
    int N_r();
    int CollSource();
    void setCollSource(int CollSource);
    void WallEroInit(int n, int N_p, int N_r);
    void AddPar(int charge, double n, double T, int source, int xy[2]);
    void stat();
};