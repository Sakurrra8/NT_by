#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include "math.h"

using namespace std;

#define Dmass 3.34e-27
#define D2mass 6.695e-27
#define CD4mass 2.665e-26
#define CD3mass 2.498e-26
#define CD2mass 2.330e-26
#define CD1mass 2.163e-26
#define Cmass 1.99518e-26
#define pi 3.14159265
#define qe 1.602e-19
#define e 2.71828

ifstream &seekline(ifstream &in, int line);

class EIRENE
{
private:
    string name_;
    string database_;
    int Fit_;
    int Num_;
    double data_[90];

public:
    EIRENE(int Fit, int Num, const string &database = "amjuel");
    double cal(double n, double T, double *H2dataforH8 = NULL);
    int fit() const;
    double data(int i);
    void out();
};
