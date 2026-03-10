#include <vector>

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

int orientation(Point &A, Point &B, Point &C);