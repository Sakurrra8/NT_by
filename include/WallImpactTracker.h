#pragma once
#include <vector>
#include <string>

double cosIncidence(double x1, double y1, double x2, double y2, double vx, double vy);
double cosIncidence3D(double x1, double y1, double x2, double y2, double vx, double vy, double vz);

class WallImpactTracker
{
public:
    // ① 默认构造：尚未初始化
    WallImpactTracker() = default;

    // ② 一步到位构造（可选）
    WallImpactTracker(int nWall,
                      int nEbin, double E_min, double E_max,
                      int nAbin, double mu_min, double mu_max);

    // ③ 延迟初始化（推荐先默认构造，再在几何就绪时调用）
    bool init(int nWall,
              int nEbin, double E_min, double E_max,
              int nAbin, double mu_min, double mu_max);

    // ④ 若墙段数后续改变：可扩/缩容；preserve=true 时保留已有统计数据，否则清零
    bool resizeWalls(int nWall, bool preserve = false);

    // ⑤ 清零当前所有统计（维度不变）
    void zero();

    // ⑥ 统计与 IO
    bool loadYieldTable(const std::string &path);
    void accumulate(int wall_id, double E, double mu, double weight);
    void accumulate(double E, double mu, double weight) { accumulate(0, E, mu, weight); }
    double Y_E(int i);
    double Y_A(int i);
    double Y_table(int i, int j);

private:
    // 维度与范围
    int NWALL_{0}, NBIN_E_{0}, NBIN_MU_{0};
    double E_MIN_{0}, E_MAX_{0}, MU_MIN_{0}, MU_MAX_{0};
    bool inited_{false};

    // 每段墙的直方图/总量：[NWALL][NBIN]
    std::vector<std::vector<double>> energyHist_;        // [NWALL][NBIN_E]
    std::vector<std::vector<double>> angleHist_;         // [NWALL][NBIN_MU]
    std::vector<double> totalFlux_;                      // [NWALL]
    std::vector<std::vector<double>> sputterEnergyHist_; // [NWALL][NBIN_E]
    std::vector<std::vector<double>> sputterAngleHist_;  // [NWALL][NBIN_MU]
    std::vector<double> totalSputter_;                   // [NWALL]

    // 溅射产额表
    bool hasYield_{false};
    std::vector<double> Y_E_, Y_A_, Y_table_; // 行主序: iE*Y_A_.size()+jA

    // 工具
    int getEnergyBin(double E) const;
    int getAngleBin(double mu) const;
    static int findBracket(const std::vector<double> &v, double x);
    double interpYield(double E_eV, double theta_deg) const;
};
