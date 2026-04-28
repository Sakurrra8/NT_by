#include "utils.h"
#include <iostream>

void Point::PointIndex(int x, int y)
{
    Index_[0] = x;
    Index_[1] = y;
}

void Point::SetPoint(double Point_x, double Point_y)
{
    x_ = Point_x;
    y_ = Point_y;
}

int *Point::Index()
{
    return Index_;
}

double Point::x()
{
    return x_;
}
double Point::y()
{
    return y_;
}

void Point::SetX(double x)
{
    x_ = x;
}
void Point::SetY(double y)
{
    y_ = y;
}

namespace Tools
{
    // 判断点 A, B, C 的方向
    int orientation(Point &A, Point &B, Point &C)
    {
        double cross = (B.x() - A.x()) * (C.y() - A.y()) - (B.y() - A.y()) * (C.x() - A.x());

        if (cross > 1e-10)
            return 1; // 逆时针 (CCW)
        else if (cross < -1e-10)
            return -1; // 顺时针 (CW)
        else
            return 0; // 共线
    }

    /////////////////////////////////////
    /// Function of Point class
    /////////////////////////////////////

    std::mt19937 &GetEngine()
    {
        // static 保证了 rd 和 gen 只会在第一次进入函数时创建
        // 之后每次调用都会返回同一个 gen 实例
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return gen;
    }

    double Random()
    {
        std::uniform_real_distribution<double> dist(0, 1);
        return dist(GetEngine());
    }

    double Maxwell(double E, double Mass)
    {
        double RR, FuncA, FuncB, Vel;
        do
        {
            RR = Tools::Random();
            FuncA = RR * RR;
            RR = Tools::Random();
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

    int PointandLine(double line_x1, double line_y1, double line_x2, double line_y2, double point_x, double point_y)
    {
        double x1, x2, y1, y2;
        x1 = line_x2 - line_x1;
        y1 = line_y2 - line_y1;
        x2 = point_x - line_x1;
        y2 = point_y - line_y1;
        if (x1 * y2 - x2 * y1 >= 0)
            return 1;
        else
            return 0;
    }

    double randomSign()
    {
        if (Tools::Random() < 0.5)
        {
            return 1.;
        }
        else
        {
            return -1.;
        }
    }

    void Vinit_cosine(std::vector<double> &Vt, double Rcos, double Rsin)
    {
        do
        {
            double ztheta = 6.283185307 * Tools::Random();
            double zcost = cos(ztheta);
            double zsint = sin(ztheta);
            double zpsi = Tools::Random() * 3.1415927;
            double zcosp = abs(cos(zpsi));
            double zsinp = sin(zpsi);
            double Vx = zsint * zsinp;
            double Vy = zcosp;
            double Vz = zsinp * zcost;
            // std::cout << pow(Vx, 2) << ", " << pow(Vy, 2) << ", " << pow(Vz, 2) << "  ";

            Vt[0] = Rcos * Vx - Rsin * Vy;
            Vt[1] = Rsin * Vx + Rcos * Vy;
            Vt[2] = Vz;
            // std::cout << pow(Vt[0], 2) << "," << pow(Vt[1], 2) << ", " << pow(Vt[2], 2) << endl;
        } while (Vt[2] == 1 || Vt[0] < -3e8);
    }

    int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
                              double p2_x, double p2_y, double p3_x, double p3_y,
                              double *i_x, double *i_y)
    {
        double s02_x, s02_y, s10_x, s10_y, s32_x, s32_y, s_numer, t_numer, denom, t;
        s10_x = p1_x - p0_x;
        s10_y = p1_y - p0_y;
        s32_x = p3_x - p2_x;
        s32_y = p3_y - p2_y;

        denom = s10_x * s32_y - s32_x * s10_y;
        if (denom == 0)
            return 0; // Collinear
        int denomPositive;
        if (denom > 0)
        {
            denomPositive = 1;
        }
        else
        {
            denomPositive = 0;
        }

        s02_x = p0_x - p2_x;
        s02_y = p0_y - p2_y;
        s_numer = s10_x * s02_y - s10_y * s02_x;
        t_numer = s32_x * s02_y - s32_y * s02_x;
        t = t_numer / denom;
        if (i_x != NULL)
            *i_x = p0_x + (t * s10_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s10_y);

        if ((s_numer < 0) == denomPositive)
            return 0; // No collision

        if ((t_numer < 0) == denomPositive)
            return 0; // No collision

        if (fabs(s_numer) > fabs(denom) || fabs(t_numer) > fabs(denom))
            return 0; // No collision
                      // Collision detected

        return 1;
    }

    double polygonarea(std::vector<Point> &polygon, int N)
    {
        int i, j;
        double area = 0;
        for (i = 0; i < N; i++)
        {
            j = (i + 1) % N;
            area += polygon[i].x() * polygon[j].y();
            area -= polygon[i].y() * polygon[j].x();
        }
        area /= 2;
        return (area < 0 ? -area : area);
    }

    double cds(double a)
    {
        if (a < 1)
            return 1. / a;
        return a;
    }
}
