#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <vector>
// #include <string.h>
#include <string>
#include <cassert>
#include <malloc.h>
#include <unistd.h>
#include <iomanip>
#include <math.h>
#include "float.h"
#include <unordered_map>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "utils.h"

using namespace std;

void copyFile(std::string source_filename, std::string target_filename);
int isIntersect(double x1, double y1, double x2, double y2,
                double x3, double y3, double x4, double y4);
int pointInPolygon(int polySides, vector<double> &polyX, vector<double> &polyY, double x, double y);
int isPolygonOverlap(vector<double> &polygon1X, vector<double> &polygon1Y, int n1,
                     vector<double> &polygon2X, vector<double> &polygon2Y, int n2);
void read_2D_array(vector<vector<double>> &p, int Nt, int Nr, string a);
void fprint_2d(int Nt, int Nr, vector<vector<double>> &value, string path);

class grid
{
private:
    double num_Point_;
    vector<vector<double>> Grid_Point_;

public:
    grid(int i = 4);
    void Setgrid(vector<double> &griddata);
    double Grid_Point(int i, int xy);
    int Ifingrid(double x, double y);
    int Ifingrid1(double x, double y);
};

class GRID
{
private:
    int N_poloidal_;
    int N_radial_;
    string Casepath_;

    vector<vector<grid>> Rectangular_Grid_;
    double Rectangular_Start_[2];
    double Rectangular_Step_[2];
    double Rectangular_End_[2];
    vector<vector<vector<int>>> Correlation_Grid_Index_;
    vector<vector<grid>> Plasma_Grid_;

    vector<vector<double>> Wall_;
    vector<vector<double>> PLasma_Grid_Boundry_;
    int Target1_Index_[2];
    int Target2_Index_[2];
    vector<vector<double>> Core_Boundry_;

    vector<int> Rectangular_num_;
    vector<int> Plasma_num_;
    int Wall_num_;
    int PLasma_Grid_Boundry_num_;
    int Core_Boundry_num_;
    vector<vector<vector<double>>> griddata_;

    vector<double> Sin_target_;
    vector<double> Cos_target_;

public:
    GRID(int num_x, int num_y);
    grid Plasma_Grid(int i, int j);
    void GRID_Read(int N_poloidal, int N_radial, string Casepath);
    void PointRead(string Path1, string Path2, int n);
    vector<double> &griddata(int i, int j);
    int IfinIndex(int i, int j, double x, double y);
    int IfinIndex1(int i, int j, double x, double y);
    int IfinWall(double x, double y);
    int IfinPlasmaGrid(double x, double y);
    int IfinCore(double x, double y);
    int PLasma_Grid_Boundry_num();
    double PLasma_Grid_Boundry(int i, int xy);
    int Wall_num();
    int Core_Boundry_num();
    double Core_Boundry(int i, int xy);
    void RectangularGridGeneration();
    void Find(double X[], int *Zone, int XY[]);
    double Rectangular_Start(int i);
    double Rectangular_Step(int i);
    double Rectangular_End(int i);
    double Wall(int num, int xy);
    int Target1_Index(int i);
    int Target2_Index(int i);
    void CalAngleTarget();
    double Sin_target(int i);
    double Cos_target(int i);
};

class Grid_extern
{
private:
    int nx_; // matrix mesh
    int ny_; // matrix mesh
    int nCv_;
    int nVx_;
    int nFc_;
    int nCmxVx_;
    double dnx_;
    double dny_;
    int Nbc_;
    int sp_Ntotal_;

    int Nwall_;

    string inputpath_mesh_;
    string Casepath_;
    string Outputpath_;

    std::vector<std::vector<double>> fcLbl_;
    std::vector<std::vector<double>> cvX_;
    std::vector<std::vector<double>> cvY_;
    std::vector<std::vector<double>> vxX_;
    std::vector<std::vector<double>> vxY_;
    std::vector<std::vector<double>> cvVx_;
    std::vector<std::vector<double>> cvVxP_;
    std::vector<std::vector<double>> fcCv_;
    std::vector<std::vector<double>> fcVx_;

    std::vector<std::vector<double>> cvBb_;
    std::vector<std::vector<double>> cvEb_;
    std::vector<std::vector<double>> cvVol_;
    std::vector<std::vector<double>> cvReg_;

    std::vector<double> boundary_core_x_e_;
    std::vector<double> boundary_core_y_e_;

    std::vector<double> boundary_first_wall_x_e_;
    std::vector<double> boundary_first_wall_y_e_;

    std::vector<double> erx_e_;
    std::vector<double> ery_e_;

    std::vector<int> core_cv_;

    std::vector<std::vector<int>> i2j2_cv_num_;
    std::vector<std::vector<std::vector<int>>> i2j2_cv_;

    vector<double> nx_point_;
    vector<double> ny_point_;
    vector<vector<double>> cv_VxX_;
    vector<vector<double>> cv_VxY_;

public:
    Grid_extern();
    ~Grid_extern();
    Grid_extern(string Casepath, string Outputpath);
    void Grid_externInit(string Casepath, string Outputpath);

    int nx(); // matrix mesh
    int ny(); // matrix mesh
    int nCv();
    int nVx();
    int nFc();
    int nCmxVx();
    int sp_Ntotal();

    string inputpath_mesh();
    string Casepath();
    string Outputpath();

    void Set_nx(int nx);
    void Set_ny(int ny);
    void Set_nCv(int nCv);
    void Set_nVx(int nVx);
    void Set_nFc(int nFc);
    void Set_nCmxVx(int nCmxVx);

    double cvBb(int i, int j);
    double cvEb(int i, int j);

    void Set_cvBb(int i, int j, double Var);
    void Set_cvEb(int i, int j, double Var);

    void read_extend_mesh();
    void get_nx_ny_point();
    void get_cv_VxX_VxY();
    void get_core_line();
    void get_first_wall_line();
    void get_core_cv();
    void get_i2j2_cv();
    void get_max_min_of_Vx(double *max_vxX, double *max_vxY, double *min_vxX, double *min_vxY);
};
// 兼容性：如果没有 C++17，可以取消下面的注释提供 clamp 的后备实现
// #ifndef __cpp_lib_clamp
// namespace std { template <class T> constexpr const T& clamp(const T& v, const T& lo, const T& hi){ return (v<lo)?lo:((hi<v)?hi:v); } }
// #endif

namespace eirene
{

    struct Node
    {
        double r{}, z{};
    };

    struct Tri
    {
        int v[3] = {-1, -1, -1};
        int neigh[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    };

    struct B2Geom
    {
        int nx = 0, ny = 0;
        std::vector<double> Rcc, Zcc; // size = nx*ny
    };

    struct Fort44
    {
        int nx = 0, ny = 0;
        int natm = 0, nmol = 0;
        std::unordered_map<std::string, std::vector<double>> fields; // key -> ny*nx
        static size_t idx3(int iy, int ix, int is, int ny, int nx)
        {
            return ((size_t)is * ny + (size_t)iy) * nx + (size_t)ix;
        }
    };

    class WALL
    {
    private:
        string Wallpath_;
        vector<double> x_;
        vector<double> y_;
        int Wall_num_;
        vector<double> Sin_wall_;
        vector<double> Cos_wall_;

    public:
        WALL();
        WALL(string Casepath);
        ~WALL();
        void Init(string casepath);
        void Wallread();
        int IfinWall(double x, double y);
        int Wall_num();
        double P(int num, int xy);
        void CalAngleWall();
        double Sin_Wall(int i);
        double Cos_Wall(int i);
    };

    class TriMesh
    {
    public:
        // 网格
        int num_Nodes_;
        int num_tris_;
        std::vector<Node> nodes_;                                    // fort.33
        std::vector<Tri> tris_;                                      // fort.34
        WALL Wall_;                                                  // wall of device
        std::vector<double> triArea_;                                // 每三角面积
        std::vector<double> triVolume_;                              // 每三角体积
        std::vector<double> triRc_, triZc_;                          // 每三角质心
        std::vector<std::vector<double>> Sin_trimesh_, Cos_trimesh_; // 边界的三角函数值

        int N_poloidal_; // number of poloidal index in plasma mesh
        int N_radial_;   // number of radial index in plasma mesh

        // 外部数据
        B2Geom b2g_; // b2fgmtry
        Fort44 f44_; // fort.44

        // === 对外接口 ===
        void load_eirene_mesh(const std::string &fort33, const std::string &fort34, const std::string &fort35, const std::string &wall);
        double TritoB2(std::vector<std::vector<double>> &Var, int i);
        double TritoB2(double **Var, int i);
        void judge_orientation();

        void read_b2fgmtry(const std::string &path); // 解析 b2fgmtry 并调用 set_b2_geometry
        void set_b2_geometry(int nx, int ny, std::vector<double> Rcc, std::vector<double> Zcc);

        void read_fort44(const std::string &path); // 启发式读取常见字段

        // 将 fort.44 的 2D 场映射到三角网（k=1 最近邻；否则反距离加权）
        std::vector<double> map_fort44_to_triangles(const std::string &field_key, int k_neighbors = 1, double p = 2.0) const;

        void mesh_find(int N_poloidal, int N_radial);
        int b2_index(int i, int j);
        int lines_info(int i, int j);
        int if_in_plasmagrid(int i);
        int Num_Target();
        double Target(int i, int j);
        int targetIndex(int i, int j); // indx of target
        double Sin_Target(int i);
        double Cos_Target(int i);
        int num_Nodes();
        int num_tris();
        double Sin_trimesh(int i, int j);
        double Cos_trimesh(int i, int j);
        bool Ifingrid(int Tri_Index, double px, double py);
        double triVolume(int i);
        void SetB2Index_find(int i, int m, int n);

    private:
        // 网格处理
        int Num_Target_;
        std::vector<std::vector<int>> TargetIndex_; // indx of target
        std::vector<std::vector<double>> Target_;   // indx of target
        std::vector<double> Sin_Target_;
        std::vector<double> Cos_Target_;
        std::vector<std::vector<int>> b2_index_; // index of b2 grid
        /// lines information: 1st parameter represent if the line is special;
        // 2nd: 0-nothing; 1-outer parget; 2-inner target; 3-other wall; 4-plasma boundrary; 5-core region. 6-PFR wall
        std::vector<std::vector<int>> lines_info_;

        void b2_find();
        void target_find(int N_poloidal, int N_radial);
        void lines_find();

        // 小工具
        static void require(bool cond, const std::string &msg);

        // fort.33/34/35
        void read_fort33(const std::string &path);
        void read_fort34(const std::string &path);
        void read_fort35(const std::string &path);
        void compute_centroids_areas();

        // fort.44
        static bool find_token_value(const std::string &s, const std::string &key, int &out);
        void readFort44_Heuristic(const std::string &path);

        // B2 近邻桶
        struct NN
        {
            int ix, iy;
            double dist2;
        };
        int nbx = 0, nbz = 0;
        double Rmin = 0, Rmax = 0, Zmin = 0, Zmax = 0, dRb = 1, dZb = 1;
        std::vector<std::vector<int>> bins; // 每桶存储 (iy*nx+ix)
        void build_b2_bins();
        std::pair<int, int> bin_of(double R, double Z) const;
        std::vector<NN> query_b2_neighbors(double R, double Z, int k) const;
    };

} // namespace eirene
