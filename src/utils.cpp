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
        thread_local std::mt19937 gen(std::random_device{}());
        return gen;
    }

    double Random()
    {
        thread_local std::uniform_real_distribution<double> dist(1e-5, 1.0 - 1e-5);
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
        // 确保 Vt 有足够的空间 (防御性编程)
        if (Vt.size() < 3)
        {
            Vt.resize(3);
        }

        // 1. 生成 [0, 1) 的均匀随机数
        double U1 = Tools::Random();
        double U2 = Tools::Random();

        // 2. 余弦分布逆变换采样 (假设法线方向为 Y 轴)
        double phi = 6.283185307179586 * U1; // 方位角 2*pi*U1

        // Y轴方向 (法向)
        double Vy = std::sqrt(U2);

        // 水平面的投影半径 (sin_theta)
        double sin_theta = std::sqrt(1.0 - U2);

        // X轴和Z轴方向
        double Vx = sin_theta * std::cos(phi);
        double Vz = sin_theta * std::sin(phi);

        // 3. 坐标旋转 (X-Y平面)
        Vt[0] = Rcos * Vx - Rsin * Vy;
        Vt[1] = Rsin * Vx + Rcos * Vy;
        Vt[2] = Vz;
    }

    int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
                              double p2_x, double p2_y, double p3_x, double p3_y,
                              double *i_x, double *i_y)
    {
        if (std::max(p0_x, p1_x) < std::min(p2_x, p3_x) ||
            std::max(p2_x, p3_x) < std::min(p0_x, p1_x) ||
            std::max(p0_y, p1_y) < std::min(p2_y, p3_y) ||
            std::max(p2_y, p3_y) < std::min(p0_y, p1_y))
            return 0;

        const double s10_x = p1_x - p0_x;
        const double s10_y = p1_y - p0_y;
        const double s32_x = p3_x - p2_x;
        const double s32_y = p3_y - p2_y;
        const double denom = s10_x * s32_y - s32_x * s10_y;
        if (denom == 0.)
            return 0; // Collinear

        const bool denomPositive = denom > 0.;
        const double s02_x = p0_x - p2_x;
        const double s02_y = p0_y - p2_y;
        const double s_numer = s10_x * s02_y - s10_y * s02_x;
        if ((s_numer < 0.) == denomPositive)
            return 0;

        const double t_numer = s32_x * s02_y - s32_y * s02_x;
        if ((t_numer < 0.) == denomPositive)
            return 0;

        const double abs_denom = std::abs(denom);
        if (std::abs(s_numer) > abs_denom || std::abs(t_numer) > abs_denom)
            return 0;

        const double t = t_numer / denom;
        if (i_x != NULL)
            *i_x = p0_x + t * s10_x;
        if (i_y != NULL)
            *i_y = p0_y + t * s10_y;
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

    /**
     * @brief 计算粒子撞击二维对称壁面后的反射速度向量 (原地修改)
     * * @param Vt    [输入/输出] 粒子的速度向量 (要求 std::vector 大小至少为3)。
     * - 漫反射时：原有速度被舍弃，直接覆盖为新的随机反弹方向。
     * - 镜面反射时：作为入射速度被读取，随后覆盖为出射速度。
     * @param Rcos  壁面倾斜角的余弦值 (cos(alpha))
     * @param Rsin  壁面倾斜角的正弦值 (sin(alpha))
     * @param type  反射类型：
     * - 0: 余弦漫反射 (Cosine Diffuse)
     * - 1: 镜面反射 (Specular)
     */
    void calculateReflectionVelocity(std::vector<double> &Vt, double Rcos, double Rsin, int type)
    {
        // 防御性编程：确保 Vt 至少有 3 个元素，防止内存越界
        if (Vt.size() < 3)
        {
            Vt.resize(3, 0.0);
        }

        if (type == 0)
        {
            // --- 0: 余弦漫反射 ---
            // 注意：多线程环境 (如 OpenMP) 下，thread_local 必不可少
            thread_local std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<double> dist(0.0, 1.0);

            double U1 = dist(rng);
            double U2 = dist(rng);

            // 1. 局部坐标系采样 (假设 Y 轴为正法向)
            double phi = 2.0 * M_PI * U1;
            double Vy = std::sqrt(U2);
            double sin_theta = std::sqrt(1.0 - U2);

            double Vx = sin_theta * std::cos(phi);
            double Vz = sin_theta * std::sin(phi);

            // 2. 绕 Z 轴进行 2D 旋转，对齐到真实的壁面法线，并直接写入 Vt
            Vt[0] = Rcos * Vx - Rsin * Vy;
            Vt[1] = Rsin * Vx + Rcos * Vy;
            Vt[2] = Vz;

            // 注: 局部采样本身已是单位向量，旋转不改变模长，因此无需再次归一化
        }
        else if (type == 1)
        {
            // --- 1: 镜面反射 ---
            // 1. 通过直接传入的 Rcos 和 Rsin，隐式推导壁面法向量 n = (-sin(a), cos(a), 0)
            double nx = -Rsin;
            double ny = Rcos;
            double nz = 0.0;

            // 2. 计算点乘 (Vt · n) (即计算入射速度在法向上的投影)
            double dot_prod = Vt[0] * nx + Vt[1] * ny + Vt[2] * nz;

            // 3. 物理有效性检查：仅当粒子朝着壁面运动 (点乘 <= 0) 时才发生反弹
            if (dot_prod <= 0.0)
            {
                // 镜面反射矢量公式： v_out = v_in - 2 * (v_in · n) * n
                Vt[0] -= 2.0 * dot_prod * nx;
                Vt[1] -= 2.0 * dot_prod * ny;
                Vt[2] -= 2.0 * dot_prod * nz;
            }

            // 4. 数值漂移修正：重新归一化以保证速度依然是严格的单位向量
            double len = std::sqrt(Vt[0] * Vt[0] + Vt[1] * Vt[1] + Vt[2] * Vt[2]);
            if (len > 1e-9)
            {
                Vt[0] /= len;
                Vt[1] /= len;
                Vt[2] /= len;
            }
        }
    }
    /**
     * @brief 判断直线与线段是否相交，并计算交点坐标
     * * @param lx1 直线上第一点的 x 坐标
     * @param ly1 直线上第一点的 y 坐标
     * @param lx2 直线上第二点的 x 坐标
     * @param ly2 直线上第二点的 y 坐标
     * @param sx1 线段第一个端点的 x 坐标
     * @param sy1 线段第一个端点的 y 坐标
     * @param sx2 线段第二个端点的 x 坐标
     * @param sy2 线段第二个端点的 y 坐标
     * @param ix  [输出] 交点的 x 坐标 (仅在返回 1 时有效)
     * @param iy  [输出] 交点的 y 坐标 (仅在返回 1 时有效)
     * @return int 1 表示相交，0 表示不相交
     */
    int getLineSegmentIntersection(double lx1, double ly1, double lx2, double ly2,
                                   double sx1, double sy1, double sx2, double sy2,
                                   double &ix, double &iy)
    {
        // 防御性编程：如果定义直线的两点重合，无法构成直线
        if (lx1 == lx2 && ly1 == ly2)
        {
            return 0;
        }

        // 1. 计算叉积 cp3 = (P2 - P1) x (S1 - P1)
        double cp3 = (lx2 - lx1) * (sy1 - ly1) - (ly2 - ly1) * (sx1 - lx1);

        // 2. 计算叉积 cp4 = (P2 - P1) x (S2 - P1)
        double cp4 = (lx2 - lx1) * (sy2 - ly1) - (ly2 - ly1) * (sx2 - lx1);

        // 3. 判断两个叉积是否异号（或其中一个为 0），确定是否相交
        if ((cp3 >= 0 && cp4 <= 0) || (cp3 <= 0 && cp4 >= 0))
        {

            double denominator = cp3 - cp4;

            // 特殊情况：如果分母接近 0，说明 cp3 == cp4 == 0，即线段完全共线于直线
            if (std::abs(denominator) < 1e-9)
            {
                // 共线时有无数个交点，这里默认返回线段的起点
                ix = sx1;
                iy = sy1;
            }
            else
            {
                // 核心算法：利用叉积的比例进行参数化插值 (t 在 0 到 1 之间)
                double t = cp3 / denominator;
                ix = sx1 + t * (sx2 - sx1);
                iy = sy1 + t * (sy2 - sy1);
            }

            return 1; // 确定相交并成功写入交点
        }

        return 0; // 不相交
    }
    /**
     * @brief 计算磁场方向与器壁法线之间的夹角（返回 0~90 度）
     * * @param bx         磁场 X 分量
     * @param by         磁场 Y 分量
     * @param bz         磁场 Z 分量
     * @param cos_wall   器壁倾斜角的 cos 值
     * @param sin_wall   器壁倾斜角的 sin 值
     * @return double    返回夹角（角度制：Degree）
     */
    double CalBFieldToWallNormalAngle(double bx, double by, double bz, double cos_wall, double sin_wall)
    {
        // 1. 计算三维磁场向量的模长 |B|
        double b_mod_sq = bx * bx + by * by + bz * bz;

        // 安全校验：防止强磁场为 0 区域导致除零崩溃
        if (b_mod_sq <= 1e-12)
        {
            return 0.0;
        }
        double b_mod = std::sqrt(b_mod_sq);

        // 2. 计算点积 (Dot Product)
        // 根据你第一道题的物理模型：器壁方向为 (cos, sin, 0)，则法线方向为 (sin, -cos, 0)
        // 磁场 B(bx, by, bz) 与法线 n(sin, -cos, 0) 的点积为： bx * sin - by * cos
        double dot_product_abs = std::abs(bx * sin_wall - by * cos_wall);

        /* 💡 备用情况提醒：
           如果你传入的 (cos_wall, sin_wall, 0) 本身就已经直接是【法线向量】了，
           请将上面的代码改成下面这一行：
           double dot_product_abs = std::abs(bx * cos_wall + by * sin_wall);
        */

        // 3. 计算 cos(theta) 值，器壁法线的模长默认已归一化（cos^2 + sin^2 = 1）
        double cos_theta = dot_product_abs / b_mod;

        // 4. 防御性截断：防止浮点数精度误差（如变成 1.00000001）导致 acos 吐出 NaN
        cos_theta = std::clamp(cos_theta, 0.0, 1.0);

        // 5. 将弧度转换为角度（适合 TRIM 数据库查表）
        const double PI = 3.14159265358979323846;
        return (180.0 / PI) * std::acos(cos_theta);
    }
    /**
     * @brief 校验并修正入射离子速度，防止其远离靶板
     * @param V          离子速度数组，长度为 3 (V[0]=Vx, V[1]=Vy, V[2]=Vz) -> 将在函数内直接修改
     * @param bx, by, bz 磁场三维分量
     * @param cos_wall   器壁倾斜角 cos 值
     * @param sin_wall   器壁倾斜角 sin 值
     */
    void AdjustIncidentVelocity(std::vector<double> &V, double bx, double by, double bz, double cos_wall, double sin_wall)
    {
        // 1. 定义器壁的法向量 n = (sin_wall, -cos_wall, 0)
        // 注意：该法向需要与你网格定义的“里外”一致。这里沿用你之前的点积公式格式：
        double nx = sin_wall;
        double ny = -cos_wall;

        // 2. 计算速度 V 和磁场 B 在器壁法线方向的投影（点积）
        double Vn = V[0] * nx + V[1] * ny;
        double Bn = bx * nx + by * ny;

        // 3. 核心物理逻辑判断：
        // 如果 Vn * Bn < 0，说明磁场引导的方向是指向靶板的，但粒子的法向速度却背离了靶板（两者反向了）
        if (Vn * Bn < 0.0)
        {
            // 4. 修正速度：将法向速度分量直接反转（相当于让它朝靶板方向推进）
            // 物理公式：V_new = V - 2 * Vn * n
            V[0] = V[0] - 2.0 * Vn * nx;
            V[1] = V[1] - 2.0 * Vn * ny;

            // Z 轴速度平行于二维靶板表面，通常不需要反转
        }
    }
}

void get_point_to_y_line(double p1x, double p1y, double p2x, double p2y, double *y0)
{
    if (p1y == p2y)
    {
        *y0 = p1y;
        return;
    }

    if (p1x == p2x)
    {
        *y0 = 1e32;
        printf("In function : get_point_to_y_line, p1x == p2x, it is not good \n");
        return;
    }

    double a, b, c;
    a = p2y - p1y;
    b = p1x - p2x;
    c = p2x * p1y - p1x * p2y;

    *y0 = -1.0 * c / b;
}