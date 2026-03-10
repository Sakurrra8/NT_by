#include "utils.h"
#include <iostream>

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
