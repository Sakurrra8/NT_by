#include "WallImpactTracker.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <H5Cpp.h>

namespace
{
    constexpr double PI = 3.14159265358979323846;
}

// ② 一步到位构造（可选）
WallImpactTracker::WallImpactTracker(int nWall,
                                     int nEbin, double E_min, double E_max,
                                     int nAbin, double mu_min, double mu_max)
{
    init(nWall, nEbin, E_min, E_max, nAbin, mu_min, mu_max);
}

// ③ 延迟初始化
bool WallImpactTracker::init(int nWall,
                             int nEbin, double E_min, double E_max,
                             int nAbin, double mu_min, double mu_max)
{
    if (nWall <= 0 || nEbin <= 0 || nAbin <= 0 || !(E_max > E_min) || !(mu_max > mu_min))
    {
        std::cerr << "[WallImpactTracker] bad init args\n";
        return false;
    }
    NWALL_ = nWall;
    NBIN_E_ = nEbin;
    NBIN_MU_ = nAbin;
    E_MIN_ = E_min;
    E_MAX_ = E_max;
    MU_MIN_ = mu_min;
    MU_MAX_ = mu_max;

    energyHist_.assign(NWALL_, std::vector<double>(NBIN_E_, 0.0));
    angleHist_.assign(NWALL_, std::vector<double>(NBIN_MU_, 0.0));
    totalFlux_.assign(NWALL_, 0.0);
    sputterEnergyHist_.assign(NWALL_, std::vector<double>(NBIN_E_, 0.0));
    sputterAngleHist_.assign(NWALL_, std::vector<double>(NBIN_MU_, 0.0));
    totalSputter_.assign(NWALL_, 0.0);

    inited_ = true;
    return true;
}

// ④ 变更墙段数（可保留已有统计）
bool WallImpactTracker::resizeWalls(int nWall, bool preserve)
{
    if (!inited_)
        return init(nWall, NBIN_E_, E_MIN_, E_MAX_, NBIN_MU_, MU_MIN_, MU_MAX_);
    if (nWall <= 0)
    {
        std::cerr << "[WallImpactTracker] resizeWalls: nWall<=0\n";
        return false;
    }
    if (!preserve)
    {
        return init(nWall, NBIN_E_, E_MIN_, E_MAX_, NBIN_MU_, MU_MIN_, MU_MAX_);
    }
    // preserve=true: 扩/缩容并拷贝交集部分
    std::vector<std::vector<double>> eh(nWall, std::vector<double>(NBIN_E_, 0.0));
    std::vector<std::vector<double>> ah(nWall, std::vector<double>(NBIN_MU_, 0.0));
    std::vector<double> tf(nWall, 0.0);
    std::vector<std::vector<double>> seh(nWall, std::vector<double>(NBIN_E_, 0.0));
    std::vector<std::vector<double>> sah(nWall, std::vector<double>(NBIN_MU_, 0.0));
    std::vector<double> ts(nWall, 0.0);

    const int copyW = std::min(NWALL_, nWall);
    for (int w = 0; w < copyW; ++w)
    {
        std::copy(energyHist_[w].begin(), energyHist_[w].end(), eh[w].begin());
        std::copy(angleHist_[w].begin(), angleHist_[w].end(), ah[w].begin());
        std::copy(sputterEnergyHist_[w].begin(), sputterEnergyHist_[w].end(), seh[w].begin());
        std::copy(sputterAngleHist_[w].begin(), sputterAngleHist_[w].end(), sah[w].begin());
        tf[w] = totalFlux_[w];
        ts[w] = totalSputter_[w];
    }
    energyHist_.swap(eh);
    angleHist_.swap(ah);
    sputterEnergyHist_.swap(seh);
    sputterAngleHist_.swap(sah);
    totalFlux_.swap(tf);
    totalSputter_.swap(ts);
    NWALL_ = nWall;
    return true;
}

// ⑤ 清零
void WallImpactTracker::zero()
{
    if (!inited_)
        return;
    for (int w = 0; w < NWALL_; ++w)
    {
        std::fill(energyHist_[w].begin(), energyHist_[w].end(), 0.0);
        std::fill(angleHist_[w].begin(), angleHist_[w].end(), 0.0);
        std::fill(sputterEnergyHist_[w].begin(), sputterEnergyHist_[w].end(), 0.0);
        std::fill(sputterAngleHist_[w].begin(), sputterAngleHist_[w].end(), 0.0);
    }
    std::fill(totalFlux_.begin(), totalFlux_.end(), 0.0);
    std::fill(totalSputter_.begin(), totalSputter_.end(), 0.0);
}

// ====== 统计 ======
void WallImpactTracker::accumulate(int wall_id, double E, double mu, double weight)
{
    if (!inited_)
        return;
    if (wall_id < 0 || wall_id >= NWALL_)
        return;
    if (E < E_MIN_ || E > E_MAX_)
        return;
    if (mu < MU_MIN_ || mu > MU_MAX_)
        return;

    const int ebin = getEnergyBin(E);
    const int abin = getAngleBin(mu);

    energyHist_[wall_id][ebin] += weight;
    angleHist_[wall_id][abin] += weight;
    totalFlux_[wall_id] += weight;

    if (hasYield_)
    {
        const double mu_c = std::max(-1.0, std::min(1.0, mu));
        const double theta_deg = std::acos(mu_c) * 180.0 / PI;
        double Y = interpYield(E, theta_deg);
        if (Y < 0.0)
            Y = 0.0;
        const double sput = weight * Y;
        /*if (E > 230)
            std::cout << E << ", " << theta_deg << ", " << Y << ", " << sput << std::endl;*/
        sputterEnergyHist_[wall_id][ebin] += sput;
        sputterAngleHist_[wall_id][abin] += sput;
        totalSputter_[wall_id] += sput;
    }
}

int WallImpactTracker::getEnergyBin(double E) const
{
    int bin = static_cast<int>((E - E_MIN_) / (E_MAX_ - E_MIN_) * NBIN_E_);
    return std::min(std::max(bin, 0), NBIN_E_ - 1);
}
int WallImpactTracker::getAngleBin(double mu) const
{
    int bin = static_cast<int>((mu - MU_MIN_) / (MU_MAX_ - MU_MIN_) * NBIN_MU_);
    return std::min(std::max(bin, 0), NBIN_MU_ - 1);
}
int WallImpactTracker::findBracket(const std::vector<double> &v, double x)
{
    if (v.size() < 2)
        return 0;
    if (x <= v.front())
        return 0;
    if (x >= v.back())
        return static_cast<int>(v.size()) - 2;
    auto it = std::upper_bound(v.begin(), v.end(), x);
    int i = static_cast<int>(it - v.begin()) - 1;
    return std::max(0, std::min(i, static_cast<int>(v.size()) - 2));
}
double WallImpactTracker::interpYield(double E_eV, double theta_deg) const
{
    if (!hasYield_ || Y_E_.size() < 2 || Y_A_.size() < 2)
    {
        std::cout << "sputter calculation error!" << std::endl;
        return 0.0;
    }
    int i = findBracket(Y_E_, E_eV);
    int j = findBracket(Y_A_, theta_deg);
    const double x0 = Y_E_[i], x1 = Y_E_[i + 1];
    const double y0 = Y_A_[j], y1 = Y_A_[j + 1];
    const double t = (x1 == x0) ? 0.0 : (E_eV - x0) / (x1 - x0);
    const double u = (y1 == y0) ? 0.0 : (theta_deg - y0) / (y1 - y0);
    auto at = [&](int ii, int jj) -> double
    {
        return Y_table_[ii * static_cast<int>(Y_A_.size()) + jj];
    };
    const double Y00 = at(i, j);
    const double Y10 = at(i + 1, j);
    const double Y01 = at(i, j + 1);
    const double Y11 = at(i + 1, j + 1);
    const double val = (1 - t) * (1 - u) * Y00 + t * (1 - u) * Y10 + (1 - t) * u * Y01 + t * u * Y11;
    return (val < 0.0) ? 0.0 : val;
}

bool WallImpactTracker::loadYieldTable(const std::string &path)
{
    Y_E_.clear();
    Y_A_.clear();
    Y_table_.clear();
    hasYield_ = false;

    std::ifstream fin(path.c_str());
    if (!fin)
    {
        std::cerr << "[WallImpactTracker] cannot open yield table: " << path << "\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string s;
    while (std::getline(fin, s))
    {
        if (!s.empty())
            lines.push_back(s);
    }
    if (lines.size() < 2)
    {
        std::cerr << "[WallImpactTracker] yield table content too short\n";
        return false;
    }

    auto isNumber = [](const std::string &z) -> bool
    {
        char *end = nullptr;
        std::strtod(z.c_str(), &end);
        return end && *end == '\0';
    };
    auto split = [](const std::string &ln)
    {
        std::istringstream iss(ln);
        std::vector<std::string> tok;
        std::string w;
        while (iss >> w)
            tok.push_back(w);
        return tok;
    };

    int data_start = 1;
    bool angles_ok = false;

    // 先找 "# energy(eV) 0 15 ..." 形式
    for (auto &ln : lines)
    {
        std::string t = ln;
        t.erase(t.begin(), std::find_if(t.begin(), t.end(), [](int ch)
                                        { return !std::isspace(ch); }));
        if (!t.empty() && t[0] == '#')
        {
            auto tok = split(t.substr(1));
            std::vector<std::string> cleaned;
            for (auto &x : tok)
            {
                std::string xl = x;
                std::transform(xl.begin(), xl.end(), xl.begin(), ::tolower);
                if (xl.find("energy") != std::string::npos || xl.find("(ev)") != std::string::npos)
                    continue;
                cleaned.push_back(x);
            }
            if (!cleaned.empty() && std::all_of(cleaned.begin(), cleaned.end(), isNumber))
            {
                for (auto &x : cleaned)
                    Y_A_.push_back(std::stod(x));
                angles_ok = true;
                data_start = 0;
                break;
            }
        }
    }
    // 否则视首行为角度
    if (!angles_ok)
    {
        auto t0 = split(lines[0]);
        auto t1 = split(lines[1]);
        bool ok0 = !t0.empty() && std::all_of(t0.begin(), t0.end(), isNumber);
        bool ok1 = !t1.empty() && std::all_of(t1.begin(), t1.end(), isNumber);
        if (!(ok0 && ok1 && (t1.size() == t0.size() + 1)))
        {
            std::cerr << "[WallImpactTracker] cannot parse angles header\n";
            return false;
        }
        for (auto &x : t0)
            Y_A_.push_back(std::stod(x));
        data_start = 1;
    }
    if (Y_A_.size() < 2)
    {
        std::cerr << "[WallImpactTracker] need >=2 angles\n";
        return false;
    }

    // 读数据行：E + Y(angles...)
    for (size_t li = data_start; li < lines.size(); ++li)
    {
        std::string ln = lines[li];
        std::string t = ln;
        t.erase(t.begin(), std::find_if(t.begin(), t.end(), [](int ch)
                                        { return !std::isspace(ch); }));
        if (!t.empty() && t[0] == '#')
            continue;

        auto tok = split(ln);
        if (tok.empty() || !isNumber(tok[0]))
            continue;

        if (tok.size() < Y_A_.size() + 1)
            tok.resize(Y_A_.size() + 1, "0");
        if (tok.size() > Y_A_.size() + 1)
            tok.resize(Y_A_.size() + 1);

        Y_E_.push_back(std::stod(tok[0]));
        for (size_t j = 0; j < Y_A_.size(); ++j)
        {
            const std::string &val = tok[1 + j];
            double y = 0.0;
            if (val == "NA" || val == "NaN" || val == "nan")
                y = 0.0;
            else
                y = std::stod(val);
            if (y < 0.0)
                y = 0.0;
            Y_table_.push_back(y);
        }
    }
    if (Y_E_.size() < 2)
    {
        std::cerr << "[WallImpactTracker] need >=2 energies\n";
        return false;
    }

    // 按能量升序重排
    std::vector<int> idx(Y_E_.size());
    for (int i = 0; i < (int)idx.size(); ++i)
        idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b)
              { return Y_E_[a] < Y_E_[b]; });

    std::vector<double> E_sorted(Y_E_.size());
    std::vector<double> Y_sorted(Y_table_.size());
    const int nA = static_cast<int>(Y_A_.size());
    for (int i = 0; i < (int)Y_E_.size(); ++i)
    {
        E_sorted[i] = Y_E_[idx[i]];
        for (int j = 0; j < nA; ++j)
        {
            Y_sorted[i * nA + j] = Y_table_[idx[i] * nA + j];
        }
    }
    Y_E_.swap(E_sorted);
    Y_table_.swap(Y_sorted);
    hasYield_ = true;
    return true;
}

double WallImpactTracker::Y_E(int i)
{
    return Y_E_[i];
}

double WallImpactTracker::Y_A(int i)
{
    return Y_A_[i];
}

double WallImpactTracker::Y_table(int i, int j)
{
    return Y_table_[i * Y_A_.size() + j];
}

// cos(theta) for a wall segment in XY-plane and ion velocity (vx, vy)
double cosIncidence(double x1, double y1, double x2, double y2, double vx, double vy)
{
    const double dx = x2 - x1, dy = y2 - y1;
    const double L = std::hypot(dx, dy);
    const double vL = std::hypot(vx, vy);
    if (L == 0.0 || vL == 0.0)
        return 0.0; // 退化情况

    // 单位切向 & 一个单位法向（切向逆时针旋转 90°）
    const double tx = dx / L, ty = dy / L;
    double nx = -ty, ny = tx;

    // 单位速度
    const double vxn = vx / vL, vyn = vy / vL;

    // 选择让粒子“朝向墙”的法向，使 -v·n >= 0
    if (vxn * nx + vyn * ny > 0.0)
    {
        nx = -nx;
        ny = -ny;
    }

    double mu = -(vxn * nx + vyn * ny);      // mu = cos(theta)
    return std::max(0.0, std::min(1.0, mu)); // clamp 到 [0,1]
}

double cosIncidence3D(double x1, double y1, double x2, double y2, double vx, double vy, double vz)
{
    const double vL = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (vL == 0.0)
        return 0.0;
    const double mu2D = cosIncidence(x1, y1, x2, y2, vx, vy);
    const double vL2D = std::hypot(vx, vy);
    // n 在 XY 平面：mu = - (v · n_hat)/|v| = (|v_xy|/|v|) * mu2D
    return (vL2D > 0.0) ? (vL2D / vL) * mu2D : 0.0;
}
