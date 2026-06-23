#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>

#define qe 1.602e-19

void get_point_to_y_line(double p1x, double p1y, double p2x, double p2y, double *y0);

class Point
{
private:
    double x_;
    double y_;
    int Index_[2];

public:
    Point() {};
    int *Index();
    double x();
    double y();
    void SetX(double x);
    void SetY(double y);
    void PointIndex(int x, int y);
    void SetPoint(double Point_x, double Point_y);
};

namespace Tools
{
    inline double sqr(double x) { return x * x; }
    int orientation(Point &A, Point &B, Point &C);

    // 获取全局唯一的引擎引用（不直接暴露变量，更安全）
    std::mt19937 &GetEngine();

    // 快捷函数：生成 [min, max] 范围内的浮点数
    double Random();
    double Maxwell(double E, double Mass);
    double MaxwellianFluxSpeed(double temperature_eV, double mass);
    double intersect(double *A, double *B, int i);
    int PointandLine(double line_x1, double line_y1, double line_x2, double line_y2, double point_x, double point_y);
    double randomSign();
    // void Vinit_cosine(std::vector<double> &Vt, double Rcos, double Rsin);
    int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
                              double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y);
    double polygonarea(std::vector<Point> &polygon, int N);
    double cds(double a);
    void calculateReflectionVelocity(std::vector<double> &Vt, double Rcos, double Rsin, int type);
    int getLineSegmentIntersection(double lx1, double ly1, double lx2, double ly2,
                                   double sx1, double sy1, double sx2, double sy2,
                                   double &ix, double &iy);
    double CalBFieldToWallNormalAngle(double bx, double by, double bz, double cos_wall, double sin_wall);
    void AdjustIncidentVelocity(std::vector<double> &V, double bx, double by, double bz, double cos_wall, double sin_wall);
}

#endif
