#include "Grid.h"

grid::grid(int num_Point)
{
    num_Point_ = num_Point;

    for (int i = 0; i < num_Point_; i++)
    {
        Grid_Point_.push_back(vector<double>());
        for (int j = 0; j < 2; j++)
        {
            Grid_Point_[i].push_back(double());
        }
    }
}

void grid::Setgrid(vector<double> &griddata)
{
    Grid_Point_[0][0] = griddata[0];
    Grid_Point_[0][1] = griddata[4];
    Grid_Point_[1][0] = griddata[1];
    Grid_Point_[1][1] = griddata[5];
    Grid_Point_[2][0] = griddata[2];
    Grid_Point_[2][1] = griddata[6];
    Grid_Point_[3][0] = griddata[3];
    Grid_Point_[3][1] = griddata[7];
}

double grid::Grid_Point(int i, int xy)
{
    return Grid_Point_[i][xy];
}

int grid::Ifingrid(double x, double y)
{
    /// judge the relationship between the points and the grid by the cross product
    /*double x1, x2, y1, y2, a[4];
    for (int i = 0; i < num_Point_; i++)
    {
        int j = i + 1;
        if (i == num_Point_ - 1)
        {
            j = 0;
        }
        // std::cout << i << '\t' << j << endl;
        x1 = Grid_Point_[j].x() - Grid_Point_[i].x();
        y1 = Grid_Point_[j].y() - Grid_Point_[i].y();
        x2 = x - Grid_Point_[i].x();
        y2 = y - Grid_Point_[i].y();
        a[i] = x1 * y2 - x2 * y1;
        if (a[i] < 0)
        {
            // std::cout << "不在" << endl;
            return 0;
        }
    }
    return 1;*/

    ///
    int count = 0;
    double secx, secy;

    for (int i = 0; i < num_Point_; i++)
    {
        double m;
        if (i < num_Point_ - 1)
        {
            m = i + 1;
        }
        else
        {
            m = 0;
        }
        if (min(Grid_Point_[i][1], Grid_Point_[m][1]) < y && max(Grid_Point_[i][1], Grid_Point_[m][1]) > y)
        {
            if (max(Grid_Point_[i][0], Grid_Point_[m][0]) < x)
            {
                count += 1;
            }
            else if (Grid_Point_[i][0] < x || Grid_Point_[m][0] < x)
            {
                if (Tools::get_line_intersection(x, y, -1, y, Grid_Point_[i][0], Grid_Point_[i][1], Grid_Point_[m][0], Grid_Point_[m][1], &secx, &secy))
                    count += 1;
            }
        }
    }
    if (count % 2 == 1)
        return true;
    else
        return false;
}

int grid::Ifingrid1(double x, double y)
{
    /// judge the relationship between the points and the grid by the cross product
    double x1, x2, y1, y2, a[4];
    for (int i = 0; i < num_Point_; i++)
    {
        int j = i + 1;
        if (i == num_Point_ - 1)
        {
            j = 0;
        }
        // std::cout << i << '\t' << j << endl;
        x1 = Grid_Point_[j][0] - Grid_Point_[i][0];
        y1 = Grid_Point_[j][1] - Grid_Point_[i][1];
        x2 = x - Grid_Point_[i][0];
        y2 = y - Grid_Point_[i][1];
        a[i] = x1 * y2 - x2 * y1;
        if (a[i] < 0)
        {
            // std::cout << "不在" << endl;
            return 0;
        }
    }
    return 1;
}

GRID::GRID(int num_x, int num_y)
{
    Plasma_num_.push_back(num_x);
    Plasma_num_.push_back(num_y);

    Rectangular_num_.push_back(100);
    Rectangular_num_.push_back(100);
    for (int i = 0; i < 100; i++)
    {
        Rectangular_Grid_.push_back(vector<grid>());
        for (int j = 0; j < 100; j++)
        {
            Rectangular_Grid_[i].push_back(grid(4));
        }
    }
    for (int i = 0; i < Plasma_num_[0]; i++)
    {
        Plasma_Grid_.push_back(vector<grid>());
        for (int j = 0; j < Plasma_num_[1]; j++)
        {
            Plasma_Grid_[i].push_back(grid(4));
        }
    }
}

grid GRID::Plasma_Grid(int i, int j)
{
    return Plasma_Grid_[i][j];
}

void GRID::GRID_Read(int N_poloidal, int N_radial, string Casepath)
{
    N_poloidal_ = N_poloidal;
    N_radial_ = N_radial;
    for (int i = 0; i < 2 * N_radial_; i++)
    {
        Sin_target_.push_back(0);
        Cos_target_.push_back(0);
    }
    Casepath_ = Casepath;
    string Path_Grid_1, Path_Grid_2;
    for (int i = 0; i < N_poloidal_; i++)
    {
        griddata_.push_back(vector<vector<double>>());
        for (int j = 0; j < N_radial_; j++)
        {
            griddata_[i].push_back(vector<double>());
            for (int k = 0; k < 8; k++)
            {
                griddata_[i][j].push_back(double());
            }
        }
    }
    Path_Grid_1 = Casepath_ + "shapedata/com_nx1";
    Path_Grid_2 = Casepath_ + "shapedata/com_ny1";
    PointRead(Path_Grid_1, Path_Grid_2, 0);
    Path_Grid_1 = Casepath_ + "shapedata/com_nx2";
    Path_Grid_2 = Casepath_ + "shapedata/com_ny2";
    PointRead(Path_Grid_1, Path_Grid_2, 1);
    Path_Grid_1 = Casepath_ + "shapedata/com_nx3";
    Path_Grid_2 = Casepath_ + "shapedata/com_ny3";
    PointRead(Path_Grid_1, Path_Grid_2, 2);
    Path_Grid_1 = Casepath_ + "shapedata/com_nx4";
    Path_Grid_2 = Casepath_ + "shapedata/com_ny4";
    PointRead(Path_Grid_1, Path_Grid_2, 3);
    for (int i = 0; i < Plasma_num_[0]; i++)
    {
        for (int j = 0; j < Plasma_num_[1]; j++)
        {
            Plasma_Grid_[i][j].Setgrid(griddata_[i][j]);
        }
    }

    /// Tokamak wall position reading(Clockwise)
    string Path_Wall = Casepath_ + "shapedata/wall_segment.txt";
    std::ifstream fp;
    fp.open(Path_Wall, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        std::cerr << "This file READING for wall have some problem!!!\n";
    }
    fp >> Wall_num_;
    Wall_num_ -= 1;
    // std::cout << Wall_num_ << endl;
    double a, b;
    Sin_Wall_.resize(Wall_num_);
    Cos_Wall_.resize(Wall_num_);
    for (int i = 0; i < Wall_num_; i++)
    {
        Wall_.push_back(vector<double>());
        Wall_[i].push_back(double());
        Wall_[i].push_back(double());
        fp >> Wall_[i][0] >> Wall_[i][1];
        // std::cout << Wall_[i][0] << Wall_[i][1] << endl;
    }
    Wall_.push_back(vector<double>());
    Wall_[Wall_num_].push_back(Wall_[0][0]);
    Wall_[Wall_num_].push_back(Wall_[0][1]);
    fp.close();
    Wall_segments_.resize(Wall_num_);
    Wall_bounds_.resize(Wall_num_);
    for (int i = 0; i < Wall_num_; i++)
    {
        Wall_segments_[i] = {Wall_[i][0], Wall_[i][1], Wall_[i + 1][0], Wall_[i + 1][1]};
        Wall_bounds_[i] = {std::min(Wall_segments_[i][0], Wall_segments_[i][2]), std::max(Wall_segments_[i][0], Wall_segments_[i][2]), std::min(Wall_segments_[i][1], Wall_segments_[i][3]), std::max(Wall_segments_[i][1], Wall_segments_[i][3])};
        double Var_temp = sqrt((Wall_[i + 1][0] - Wall_[i][0]) * (Wall_[i + 1][0] - Wall_[i][0]) + (Wall_[i + 1][1] - Wall_[i][1]) * (Wall_[i + 1][1] - Wall_[i][1]));
        Cos_Wall_[i] = (Wall_[i + 1][0] - Wall_[i][0]) / Var_temp;
        Sin_Wall_[i] = (Wall_[i + 1][1] - Wall_[i][1]) / Var_temp;
    }
    BuildWallCandidateBins();
    // std::cout << Wall_.size() - 1 << endl;
    ofstream out("doc/Wall_.txt");
    for (int i = 0; i < Wall_.size(); i++)
    {
        out << Wall_[i][0] << '\t' << Wall_[i][1] << endl;
    }
    out.close();

    RectangularGridGeneration();

    int count = 0;
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    PLasma_Grid_Boundry_.push_back(vector<double>());
    PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[1][1].Grid_Point(0, 0));
    PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[1][1].Grid_Point(0, 1));
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    Target1_Index_[0] = count;
    for (int i = 1; i < N_radial_ - 1; i++)
    {
        Target1_Index_[1] = count;
        PLasma_Grid_Boundry_.push_back(vector<double>());
        PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[1][i].Grid_Point(3, 0));
        PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[1][i].Grid_Point(3, 1));
    }
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    for (int i = 1; i < N_poloidal_ - 1; i++)
    {
        PLasma_Grid_Boundry_.push_back(vector<double>());
        PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[i][N_radial_ - 2].Grid_Point(2, 0));
        PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[i][N_radial_ - 2].Grid_Point(2, 1));
    }
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    Target2_Index_[0] = count;
    for (int i = N_radial_ - 2; i > 0; i--)
    {
        Target2_Index_[1] = count;
        PLasma_Grid_Boundry_.push_back(vector<double>());
        PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(1, 0));
        PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(1, 1));
    }
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    for (int i = N_poloidal_ - 2; i > 72; i--)
    {
        PLasma_Grid_Boundry_.push_back(vector<double>());
        PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[i][1].Grid_Point(0, 0));
        PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[i][1].Grid_Point(0, 1));
    }
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    for (int i = 24; i > 0; i--)
    {
        PLasma_Grid_Boundry_.push_back(vector<double>());
        PLasma_Grid_Boundry_[count].push_back(Plasma_Grid_[i][1].Grid_Point(0, 0));
        PLasma_Grid_Boundry_[count++].push_back(Plasma_Grid_[i][1].Grid_Point(0, 1));
    }
    PLasma_Grid_Boundry_num_ = PLasma_Grid_Boundry_.size() - 1;
    Plasma_boundary_segments_.resize(PLasma_Grid_Boundry_num_);
    Plasma_boundary_bounds_.resize(PLasma_Grid_Boundry_num_);
    for (int i = 0; i < PLasma_Grid_Boundry_num_; i++)
    {
        Plasma_boundary_segments_[i] = {PLasma_Grid_Boundry_[i][0], PLasma_Grid_Boundry_[i][1], PLasma_Grid_Boundry_[i + 1][0], PLasma_Grid_Boundry_[i + 1][1]};
        Plasma_boundary_bounds_[i] = {std::min(Plasma_boundary_segments_[i][0], Plasma_boundary_segments_[i][2]), std::max(Plasma_boundary_segments_[i][0], Plasma_boundary_segments_[i][2]), std::min(Plasma_boundary_segments_[i][1], Plasma_boundary_segments_[i][3]), std::max(Plasma_boundary_segments_[i][1], Plasma_boundary_segments_[i][3])};
    }
    // std::cout << PLasma_Grid_Boundry_.size() << endl;
    /*out.open("doc/plasmaGridBoundry.txt");
    for (int i = 0; i < PLasma_Grid_Boundry_.size(); i++)
    {
        out << PLasma_Grid_Boundry_[i][0] << '\t' << PLasma_Grid_Boundry_[i][1] << endl;
    }
    out.close();*/

    count = 0;
    Core_Boundry_.push_back(vector<double>());
    Core_Boundry_[count].push_back(Plasma_Grid_[25][1].Grid_Point(0, 0));
    Core_Boundry_[count++].push_back(Plasma_Grid_[25][1].Grid_Point(0, 1));
    for (int i = 25; i < 73; i++)
    {
        Core_Boundry_.push_back(vector<double>());
        Core_Boundry_[count].push_back(Plasma_Grid_[i][1].Grid_Point(1, 0));
        Core_Boundry_[count++].push_back(Plasma_Grid_[i][1].Grid_Point(1, 1));
    }
    Core_Boundry_num_ = Core_Boundry_.size() - 1;
    /*out.open("doc/EAST_CoreBoundry.txt");
    for (int i = 0; i < Core_Boundry_.size(); i++)
    {
        out << Core_Boundry_[i][0] << '\t' << Core_Boundry_[i][1] << endl;
    }
    out.close();*/
}

void GRID::PointRead(string Path1, string Path2, int n)
{
    std::ifstream fp;             // file stream
    fp.open(Path1, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path1 + " have some problem!!!\n";
    }
    for (int i = 0; i < N_poloidal_; i++)
        for (int j = 0; j < N_radial_; j++)
        {
            fp >> griddata_[i][j][n];
        }
    fp.close();
    fp.open(Path2, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        cerr << "This file READING for " + Path2 + " have some problem!!!\n";
    }
    for (int i = 0; i < N_poloidal_; i++)
        for (int j = 0; j < N_radial_; j++)
        {
            fp >> griddata_[i][j][n + 4];
        }
    fp.close();
}

vector<double> &GRID::griddata(int i, int j)
{
    return griddata_[i][j];
}

int GRID::IfinIndex(int i, int j, double x, double y)
{
    if (i < 0 && j < 0)
    {
        if (IfinPlasmaGrid(x, y))
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else if (Plasma_Grid_[i][j].Ifingrid(x, y))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int GRID::IfinIndex1(int i, int j, double x, double y)
{
    if (i < 0 && j < 0)
    {
        if (IfinPlasmaGrid(x, y))
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else if (Plasma_Grid_[i][j].Ifingrid1(x, y))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int GRID::IfinWall(double x, double y)
{
    int count = 0;
    double secx, secy;
    for (int i = 0; i < Wall_.size() - 1; i++)
    {
        if (min(Wall_[i][1], Wall_[i + 1][1]) < y && max(Wall_[i][1], Wall_[i + 1][1]) > y)
        {
            if (max(Wall_[i][0], Wall_[i + 1][0]) < x)
            {
                count += 1;
            }
            else if (Wall_[i][0] < x || Wall_[i + 1][0] < x)
            {
                if (Tools::get_line_intersection(x, y, -1, y, Wall_[i][0], Wall_[i][1], Wall_[i + 1][0], Wall_[i + 1][1], &secx, &secy))
                    count += 1;
            }
        }
    }
    if (count % 2 == 1)
        return true;
    else
        return false;
}

int GRID::IfinPlasmaGrid(double x, double y)
{
    int count = 0;
    double secx, secy;
    for (int i = 0; i < PLasma_Grid_Boundry_.size() - 1; i++)
    {
        if (min(PLasma_Grid_Boundry_[i][1], PLasma_Grid_Boundry_[i + 1][1]) < y && max(PLasma_Grid_Boundry_[i][1], PLasma_Grid_Boundry_[i + 1][1]) > y)
        {
            if (max(PLasma_Grid_Boundry_[i][0], PLasma_Grid_Boundry_[i + 1][0]) < x)
            {
                count += 1;
            }
            else if (PLasma_Grid_Boundry_[i][0] < x || PLasma_Grid_Boundry_[i + 1][0] < x)
            {
                if (Tools::get_line_intersection(x, y, -1, y, PLasma_Grid_Boundry_[i][0], PLasma_Grid_Boundry_[i][1], PLasma_Grid_Boundry_[i + 1][0], PLasma_Grid_Boundry_[i + 1][1], &secx, &secy))
                    count += 1;
            }
        }
    }
    if (count % 2 == 1)
        return true;
    else
        return false;
}

int GRID::IfinCore(double x, double y)
{
    int count = 0;
    double secx, secy;
    for (int i = 0; i < Core_Boundry_.size() - 1; i++)
    {
        if (min(Core_Boundry_[i][1], Core_Boundry_[i + 1][1]) < y && max(Core_Boundry_[i][1], Core_Boundry_[i + 1][1]) > y)
        {
            if (max(Core_Boundry_[i][0], Core_Boundry_[i + 1][0]) < x)
            {
                count += 1;
            }
            else if (Core_Boundry_[i][0] < x || Core_Boundry_[i + 1][0] < x)
            {
                if (Tools::get_line_intersection(x, y, -1, y, Core_Boundry_[i][0], Core_Boundry_[i][1], Core_Boundry_[i + 1][0], Core_Boundry_[i + 1][1], &secx, &secy))
                    count += 1;
            }
        }
    }
    if (count % 2 == 1)
        return true;
    else
        return false;
}

void GRID::RectangularGridGeneration()
{
    Rectangular_Start_[0] = Wall_[0][0];
    Rectangular_End_[0] = Wall_[0][0];
    Rectangular_Start_[1] = Wall_[0][1];
    Rectangular_End_[1] = Wall_[0][1];
    for (int i = 1; i < Wall_.size(); i++)
    {
        Rectangular_Start_[0] = min(Rectangular_Start_[0], Wall_[i][0]);
        Rectangular_End_[0] = max(Rectangular_End_[0], Wall_[i][0]);
        Rectangular_Start_[1] = min(Rectangular_Start_[1], Wall_[i][1]);
        Rectangular_End_[1] = max(Rectangular_End_[1], Wall_[i][1]);
    }
    Rectangular_Step_[0] = (Rectangular_End_[0] - Rectangular_Start_[0]) / (double)Rectangular_num_[0];
    Rectangular_Step_[1] = (Rectangular_End_[1] - Rectangular_Start_[1]) / (double)Rectangular_num_[1];
    for (int i = 0; i < Rectangular_num_[0]; i++)
    {
        Correlation_Grid_Index_.push_back(vector<vector<int>>());
        for (int j = 0; j < Rectangular_num_[1]; j++)
        {
            Correlation_Grid_Index_[i].push_back(vector<int>());
        }
    }
    double index_x1, index_y1, index_x2, index_y2;
    for (int i = 1; i < Plasma_num_[0] - 1; i++)
    {
        for (int j = 1; j < Plasma_num_[1] - 1; j++)
        {
            index_x1 = (Plasma_Grid_[i][j].Grid_Point(0, 0) - Rectangular_Start_[0]) / (double)Rectangular_Step_[0];
            index_x2 = (Plasma_Grid_[i][j].Grid_Point(0, 0) - Rectangular_Start_[0]) / (double)Rectangular_Step_[0];
            index_y1 = (Plasma_Grid_[i][j].Grid_Point(0, 1) - Rectangular_Start_[1]) / (double)Rectangular_Step_[1];
            index_y2 = (Plasma_Grid_[i][j].Grid_Point(0, 1) - Rectangular_Start_[1]) / (double)Rectangular_Step_[1];
            for (int k = 1; k < 4; k++)
            {
                index_x1 = min(index_x1, (Plasma_Grid_[i][j].Grid_Point(k, 0) - Rectangular_Start_[0]) / (double)Rectangular_Step_[0]);
                index_x2 = max(index_x2, (Plasma_Grid_[i][j].Grid_Point(k, 0) - Rectangular_Start_[0]) / (double)Rectangular_Step_[0]);
                index_y1 = min(index_y1, (Plasma_Grid_[i][j].Grid_Point(k, 1) - Rectangular_Start_[1]) / (double)Rectangular_Step_[1]);
                index_y2 = max(index_y2, (Plasma_Grid_[i][j].Grid_Point(k, 1) - Rectangular_Start_[1]) / (double)Rectangular_Step_[1]);
            }
            for (int m = index_x1; m <= index_x2; m++)
            {
                for (int n = index_y1; n <= index_y2; n++)
                {
                    Correlation_Grid_Index_[m][n].push_back(i * Plasma_num_[1] + j);
                }
            }
        }
    }
    /*for (int i = 0; i < Rectangular_num_[0]; i++)
    {
        for (int j = 0; j < Rectangular_num_[1]; j++)
        {
            for (int k = 0; k < Correlation_Grid_Index_[i][j].size(); k++)
            {
                std::cout << Correlation_Grid_Index_[i][j][k] << '\t';
            }
            std::cout << endl;
        }
    }*/

    // std::cout << a << '\t' << b << '\t' << c << '\t' << d << endl;
}

int GRID::PLasma_Grid_Boundry_num()
{
    return PLasma_Grid_Boundry_num_;
}

double GRID::PLasma_Grid_Boundry(int i, int xy)
{
    return PLasma_Grid_Boundry_[i][xy];
}

const std::array<double, 4> &GRID::PlasmaBoundarySegment(int num) const
{
    return Plasma_boundary_segments_[num];
}

const std::array<double, 4> &GRID::PlasmaBoundaryBounds(int num) const
{
    return Plasma_boundary_bounds_[num];
}

int GRID::Wall_num()
{
    return Wall_num_;
}

int GRID::Core_Boundry_num()
{
    return Core_Boundry_num_;
}

double GRID::Core_Boundry(int i, int xy)
{
    return Core_Boundry_[i][xy];
}

const std::array<double, 4> &GRID::CoreBoundarySegment(int num) const
{
    return Core_boundary_segments_[num];
}

const std::array<double, 4> &GRID::CoreBoundaryBounds(int num) const
{
    return Core_boundary_bounds_[num];
}

void GRID::Find(double X[], int *Zone, int XY[])
{
    int IfFind = 0;
    int index_x, index_y;

    // test if the particle is in the near grid

    if (XY[0] > 0 && XY[0] < 97 && XY[1] > 0 && XY[1] < 37)
    {
        if (IfinIndex(XY[0], XY[1], X[0], X[1]))
        {
            return;
        }
    }
    if (XY[0] > 1 && XY[0] < N_poloidal_ - 2 && XY[1] > 1 && XY[1] < N_radial_ - 2)
    {
        if (IfinIndex(XY[0] + 1, XY[1], X[0], X[1]))
        {
            XY[0] += 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0] - 1, XY[1], X[0], X[1]))
        {
            XY[0] -= 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0], XY[1] + 1, X[0], X[1]))
        {
            XY[1] += 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0], XY[1] - 1, X[0], X[1]))
        {
            XY[1] -= 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0] + 1, XY[1] + 1, X[0], X[1]))
        {
            XY[0] += 1;
            XY[1] += 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0] - 1, XY[1] - 1, X[0], X[1]))
        {
            XY[0] -= 1;
            XY[1] -= 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0] - 1, XY[1] + 1, X[0], X[1]))
        {
            XY[0] -= 1;
            XY[1] += 1;
            IfFind = 1;
        }
        else if (IfinIndex(XY[0] + 1, XY[1] - 1, X[0], X[1]))
        {
            XY[0] += 1;
            XY[1] -= 1;
            IfFind = 1;
        }
        if (IfFind)
        {
            if (XY[1] >= N_radial_ / 2 && XY[1] <= N_radial_ - 2)
                *Zone = 3;
            else if (XY[0] <= 24)
                *Zone = 4;
            else if (XY[0] <= 72)
                *Zone = 2;
            else
                *Zone = 5;
            return;
        }
    }

    /// Calculate the index of the particle
    if (X[0] <= Rectangular_Start_[0] || X[0] >= Rectangular_End_[0] || X[1] <= Rectangular_Start_[1] || X[1] >= Rectangular_End_[1])
    {
        *Zone = 7;
        XY[0] = -1;
        XY[1] = -1;
        IfFind += 1;
        return;
    }
    else
    {
        index_x = (X[0] - Rectangular_Start_[0]) / Rectangular_Step_[0];
        index_y = (X[1] - Rectangular_Start_[1]) / Rectangular_Step_[1];

        /*if (Correlation_Grid_Index_[index_x][index_y].empty())
        {
            *Zone = 7;
            XY[0] = -1;
            XY[1] = -1;
            IfFind += 1;
            return;
        }*/
        for (int i = 0; i < Correlation_Grid_Index_[index_x][index_y].size(); i++)
        {
            if (IfinIndex(Correlation_Grid_Index_[index_x][index_y][i] / Plasma_num_[1], Correlation_Grid_Index_[index_x][index_y][i] % Plasma_num_[1], X[0], X[1]))
            {
                XY[0] = Correlation_Grid_Index_[index_x][index_y][i] / Plasma_num_[1];
                XY[1] = Correlation_Grid_Index_[index_x][index_y][i] % Plasma_num_[1];

                if (XY[1] >= N_radial_ / 2 && XY[1] <= N_radial_ - 2)
                    *Zone = 3;
                else if (XY[0] <= 24)
                    *Zone = 4;
                else if (XY[0] <= 72)
                    *Zone = 2;
                else
                    *Zone = 5;
                IfFind += 1;
                return;

                // std::cout << IfFind << ',' << *Zone << ',' << XY[0] << ',' << XY[1] << '\t';
            }
        }
    }
    // std::cout << endl;

    if (IfinCore(X[0], X[1]))
    {
        *Zone = 1;
        XY[0] = -1;
        XY[1] = -1;
        IfFind += 1;
        return;
    }
    if (IfinWall(X[0], X[1]))
    {
        *Zone = 6;
        XY[0] = -1;
        XY[1] = -1;
        IfFind += 1;
        return;
    }
    /*else if (!IfinPlasmaGrid(X[0], X[1]))
    {
        *Zone = 6;
        XY[0] = -1;
        XY[1] = -1;
        IfFind += 1;
        return;
    }*/
    else
    {
        *Zone = 7;
        XY[0] = -1;
        XY[1] = -1;
        IfFind += 1;
        return;
    }
    /*if (IfFind)
        std::cout << *Zone << ',' << XY[0] << ',' << XY[1] << ',' << IfFind << endl;
    else
        std::cout << *Zone << ',' << "Particle was not Found" << endl;*/
}

double GRID::Rectangular_Start(int i)
{
    return Rectangular_Start_[i];
}

double GRID::Rectangular_Step(int i)
{
    return Rectangular_Step_[i];
}

double GRID::Rectangular_End(int i)
{
    return Rectangular_End_[i];
}

/// @brief Wall
/// @param num element of wall
/// @param xy 0 for x, 1 for y
double GRID::Wall(int num, int xy)
{
    if (num > Wall_num_ || num < 0)
    {
        std::cout << "Wall index error" << endl;
        exit(0);
    }
    return Wall_[num][xy];
}

const std::array<double, 4> &GRID::WallSegment(int num) const
{
    return Wall_segments_[num];
}

const std::array<double, 4> &GRID::WallBounds(int num) const
{
    return Wall_bounds_[num];
}

void GRID::BuildWallCandidateBins()
{
    Wall_bins_.clear();
    if (Wall_bounds_.empty())
        return;

    double min_x = Wall_bounds_[0][0];
    double max_x = Wall_bounds_[0][1];
    double min_y = Wall_bounds_[0][2];
    double max_y = Wall_bounds_[0][3];
    for (const auto &bounds : Wall_bounds_)
    {
        min_x = std::min(min_x, bounds[0]);
        max_x = std::max(max_x, bounds[1]);
        min_y = std::min(min_y, bounds[2]);
        max_y = std::max(max_y, bounds[3]);
    }

    const int bin_count = std::max(8, static_cast<int>(std::sqrt(static_cast<double>(Wall_bounds_.size()))));
    Wall_bin_nx_ = bin_count;
    Wall_bin_ny_ = bin_count;
    Wall_bin_min_x_ = min_x;
    Wall_bin_min_y_ = min_y;
    Wall_bin_dx_ = (max_x - min_x) / Wall_bin_nx_;
    Wall_bin_dy_ = (max_y - min_y) / Wall_bin_ny_;
    if (Wall_bin_dx_ <= 0.0 || Wall_bin_dy_ <= 0.0)
    {
        Wall_bins_.clear();
        return;
    }

    Wall_bins_.assign(Wall_bin_nx_ * Wall_bin_ny_, std::vector<int>());
    for (int i = 0; i < Wall_num_; i++)
    {
        const auto &bounds = Wall_bounds_[i];
        int ix0 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((bounds[0] - Wall_bin_min_x_) / Wall_bin_dx_)));
        int ix1 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((bounds[1] - Wall_bin_min_x_) / Wall_bin_dx_)));
        int iy0 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((bounds[2] - Wall_bin_min_y_) / Wall_bin_dy_)));
        int iy1 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((bounds[3] - Wall_bin_min_y_) / Wall_bin_dy_)));
        for (int ix = ix0; ix <= ix1; ix++)
            for (int iy = iy0; iy <= iy1; iy++)
                Wall_bins_[ix * Wall_bin_ny_ + iy].push_back(i);
    }
}

void GRID::WallCandidateSegments(double min_x, double max_x, double min_y, double max_y, std::vector<int> &candidates) const
{
    candidates.clear();
    if (Wall_bins_.empty() || Wall_bin_dx_ <= 0.0 || Wall_bin_dy_ <= 0.0)
    {
        candidates.reserve(Wall_num_);
        for (int i = 0; i < Wall_num_; i++)
            candidates.push_back(i);
        return;
    }

    int ix0 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((min_x - Wall_bin_min_x_) / Wall_bin_dx_)));
    int ix1 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((max_x - Wall_bin_min_x_) / Wall_bin_dx_)));
    int iy0 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((min_y - Wall_bin_min_y_) / Wall_bin_dy_)));
    int iy1 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((max_y - Wall_bin_min_y_) / Wall_bin_dy_)));

    static thread_local std::vector<int> seen;
    static thread_local int stamp = 1;
    if (static_cast<int>(seen.size()) < Wall_num_)
        seen.assign(Wall_num_, 0);
    if (stamp == std::numeric_limits<int>::max())
    {
        std::fill(seen.begin(), seen.end(), 0);
        stamp = 1;
    }
    const int current_stamp = stamp++;

    for (int ix = ix0; ix <= ix1; ix++)
    {
        for (int iy = iy0; iy <= iy1; iy++)
        {
            const auto &bin = Wall_bins_[ix * Wall_bin_ny_ + iy];
            for (int index : bin)
            {
                if (seen[index] != current_stamp)
                {
                    seen[index] = current_stamp;
                    candidates.push_back(index);
                }
            }
        }
    }
}

int GRID::Target1_Index(int i)
{
    return Target1_Index_[i];
}

int GRID::Target2_Index(int i)
{
    return Target2_Index_[i];
}
void GRID::CalAngleTarget()
{
    for (int i = 0; i < N_radial_; i++)
    {
        double deltaX, deltaY, deltaL;
        deltaX = Plasma_Grid_[1][i].Grid_Point(3, 0) - Plasma_Grid_[1][i].Grid_Point(0, 0);
        deltaY = Plasma_Grid_[1][i].Grid_Point(3, 1) - Plasma_Grid_[1][i].Grid_Point(0, 1);
        deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
        Cos_target_[i] = -deltaX / deltaL;
        Sin_target_[i] = -deltaY / deltaL;
        deltaX = Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(1, 0) - Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(2, 0);
        deltaY = Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(1, 1) - Plasma_Grid_[N_poloidal_ - 2][i].Grid_Point(2, 1);
        deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
        Cos_target_[i + N_radial_] = -deltaX / deltaL;
        Sin_target_[i + N_radial_] = -deltaY / deltaL;
    }
    /*ofstream out_temp3("doc/angle_target2.txt");
    for (int i = 0; i < N_radial_ * 2; i++)
    {
        out_temp3 << i << "\t" << Cos_target_[i] << "\t" << Sin_target_[i] << endl;
    }
    out_temp3.close();*/
}
double GRID::Sin_target(int i)
{
    return Sin_target_[i];
}
double GRID::Cos_target(int i)
{
    return Cos_target_[i];
}
double GRID::Sin_Wall(int i)
{
    return Sin_Wall_[i];
}
double GRID::Cos_Wall(int i)
{
    return Cos_Wall_[i];
}
/*int Tools::get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
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
}*/

Grid_extern::Grid_extern() {}
Grid_extern::~Grid_extern() {}

Grid_extern::Grid_extern(string Casepath, string Outputpath)
{
    Grid_externInit(Casepath, Outputpath);
}

void Grid_extern::Grid_externInit(string Casepath, string Outputpath)
{
    Casepath_ = Casepath;
    inputpath_mesh_ = Casepath_ + "configuration/";
    Outputpath_ = Outputpath;
}

void Grid_extern::read_extend_mesh()
{
    read_2D_array(fcLbl_, nFc_, 1, inputpath_mesh_ + "fcLbl.txt");
    read_2D_array(cvX_, nCv_, 1, inputpath_mesh_ + "cvX.txt");
    read_2D_array(cvY_, nCv_, 1, inputpath_mesh_ + "cvY.txt");
    read_2D_array(vxX_, nVx_, 1, inputpath_mesh_ + "vxX.txt");
    read_2D_array(vxY_, nVx_, 1, inputpath_mesh_ + "vxY.txt");
    read_2D_array(cvVx_, nCmxVx_, 1, inputpath_mesh_ + "cvVx.txt");
    read_2D_array(cvVxP_, nCmxVx_, 2, inputpath_mesh_ + "cvVxP.txt");
    read_2D_array(fcCv_, nFc_, 2, inputpath_mesh_ + "fcCv.txt");
    read_2D_array(fcVx_, nFc_, 2, inputpath_mesh_ + "fcVx.txt");
    read_2D_array(cvBb_, nCv_, 4, inputpath_mesh_ + "cvBb.txt");
    read_2D_array(cvEb_, nCv_, 3, inputpath_mesh_ + "cvEb.txt");
    read_2D_array(cvVol_, nCv_, 1, inputpath_mesh_ + "cvVol.txt");
    read_2D_array(cvReg_, nCv_, 1, inputpath_mesh_ + "cvReg.txt");

    for (int i = 0; i < nCv_; i++)
    {
        erx_e_.push_back(-cvEb_[i][1] / sqrt(cvEb_[i][0] * cvEb_[i][0] + cvEb_[i][1] * cvEb_[i][1]));
        ery_e_.push_back(cvEb_[i][0] / sqrt(cvEb_[i][0] * cvEb_[i][0] + cvEb_[i][1] * cvEb_[i][1]));
    }

    string str_in, str_out;

    str_in = inputpath_mesh_ + "fcLbl.txt";
    str_out = Outputpath_ + "mesh/fcLbl.txt";
    std::cout << str_in << "\t" << str_out << endl;
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvX.txt";
    str_out = Outputpath_ + "mesh/cvX.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvY.txt";
    str_out = Outputpath_ + "mesh/cvY.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "vxX.txt";
    str_out = Outputpath_ + "mesh/vxX.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "vxY.txt";
    str_out = Outputpath_ + "mesh/vxY.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvVx.txt";
    str_out = Outputpath_ + "mesh/cvVx.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvVxP.txt";
    str_out = Outputpath_ + "mesh/cvVxP.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "fcCv.txt";
    str_out = Outputpath_ + "mesh/fcCv.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "fcVx.txt";
    str_out = Outputpath_ + "mesh/fcVx.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvBb.txt";
    str_out = Outputpath_ + "mesh/cvBb.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvEb.txt";
    str_out = Outputpath_ + "mesh/cvEb.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvVol.txt";
    str_out = Outputpath_ + "mesh/cvVol.txt";
    copyFile(str_in, str_out);

    str_in = inputpath_mesh_ + "cvReg.txt";
    str_out = Outputpath_ + "mesh/cvReg.txt";
    copyFile(str_in, str_out);
}

void Grid_extern::get_nx_ny_point()
{
    int nx, ny;
    nx = nx_ + 1;
    ny = ny_ + 1;
    for (int i = 0; i < nx; i++)
    {
        nx_point_.push_back(0);
    }
    for (int i = 0; i < ny; i++)
    {
        ny_point_.push_back(0);
    }

    double max_vxX, max_vxY, min_vxX, min_vxY;
    get_max_min_of_Vx(&max_vxX, &max_vxY, &min_vxX, &min_vxY);

    for (int i = 0; i <= nx - 1; i = i + 1)
    {
        nx_point_[i] = min_vxX + (double)i / (nx - 1) * (max_vxX - min_vxX);
    }

    for (int i = 0; i <= ny - 1; i = i + 1)
    {
        ny_point_[i] = min_vxY + (double)i / (ny - 1) * (max_vxY - min_vxY);
    }

    dnx_ = 1.0 / (double)(nx - 1) * (max_vxX - min_vxX);
    dny_ = 1.0 / (double)(ny - 1) * (max_vxY - min_vxY);

    string path1, path2;
    path1 = Outputpath_ + "//mesh//nx_point.txt";
    path2 = Outputpath_ + "//mesh//ny_point.txt";
    std::ofstream fout1(path1, std::ios::out);
    std::ofstream fout2(path2, std::ios::out);
    for (int i = 0; i <= nx - 1; i = i + 1)
    {
        fout1 << nx_point_[i] << endl;
        ;
    }
    for (int i = 0; i <= ny - 1; i = i + 1)
    {
        fout2 << nx_point_[i] << endl;
    }
    fout1.close();
    fout2.close();
}

void Grid_extern::get_cv_VxX_VxY()
{

    for (int i = 0; i < nCv_; i++)
    {
        cv_VxX_.push_back(vector<double>());
        cv_VxY_.push_back(vector<double>());
    }

    string path1 = Outputpath_ + "/mesh//cv_VxX.txt";
    string path2 = Outputpath_ + "/mesh//cv_VxY.txt";
    std::ofstream fout1(path1, std::ios::out);
    std::ofstream fout2(path2, std::ios::out);
    for (int i = 0; i <= nCv_ - 1; i = i + 1)
    {
        for (int j = 0; j <= (int)cvVxP_[i][1] - 1; j = j + 1)
        {
            cv_VxX_[i].push_back(vxX_[(int)cvVx_[(int)cvVxP_[i][0] + j - 1][0] - 1][0]);
            cv_VxY_[i].push_back(vxY_[(int)cvVx_[(int)cvVxP_[i][0] + j - 1][0] - 1][0]);
            fout1 << cv_VxX_[i][j] << "\t";
            fout2 << cv_VxY_[i][j] << "\t";
        }
        fout1 << endl;
        fout2 << endl;
    }
    fout1.close();
    fout2.close();
}

void Grid_extern::get_core_line()
{
    int nEdges = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] == -21 || fcLbl_[i][0] == -25)
        {
            nEdges = nEdges + 1;
        }
    }

    Nbc_ = nEdges;
    for (int i = 0; i < nEdges; i++)
    {
        boundary_core_x_e_.push_back(0.);
        boundary_core_y_e_.push_back(0.);
    }

    vector<int> Rows(nEdges);
    vector<int> Visited(nEdges);
    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        Visited[i] = 0;
    }
    Visited[0] = 1;

    vector<vector<int>> Edges;
    for (int i = 0; i < nEdges; i++)
    {
        Edges.push_back(vector<int>());
        for (int j = 0; j < 2; j++)
        {
            Edges[i].push_back(0);
        }
    }

    nEdges = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] == -21 || fcLbl_[i][0] == -25)
        {
            Rows[nEdges] = i;
            Edges[nEdges][0] = (int)fcVx_[i][0];
            Edges[nEdges][1] = (int)fcVx_[i][1];

            nEdges = nEdges + 1;
        }
    }

    int currentVertex;
    currentVertex = Edges[0][0];

    int nextVertex;

    vector<int> sequence((nEdges + 1));
    sequence[0] = Edges[0][0];

    for (int i = 1; i <= nEdges; i = i + 1)
    {
        if (i == nEdges)
        {
            Visited[0] = 0;
        }

        int idx = 0;
        for (int j = 0; j <= nEdges - 1; j = j + 1)
        {
            if ((Edges[j][0] == currentVertex || Edges[j][1] == currentVertex) && Visited[j] == 0)
            {
                idx = j;
                continue;
            }
        }

        nextVertex = Edges[idx][0] + Edges[idx][1] - currentVertex;
        currentVertex = nextVertex;
        sequence[i] = nextVertex;
        Visited[idx] = 1;
    }

    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        boundary_core_x_e_[i] = vxX_[sequence[i] - 1][0];
        boundary_core_y_e_[i] = vxY_[sequence[i] - 1][0];
    }

    string path1 = Outputpath_ + "/mesh//boundary_core_x_e.txt";
    string path2 = Outputpath_ + "/mesh//boundary_core_y_e.txt";
    std::ofstream fout1(path1, std::ios::out);
    std::ofstream fout2(path2, std::ios::out);
    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        fout1 << boundary_core_x_e_[i] << endl;
        fout2 << boundary_core_y_e_[i] << endl;
    }
    fout1.close();
    fout2.close();
}

void Grid_extern::get_first_wall_line()
{
    int nEdges = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] >= 1 && fcLbl_[i][0] <= 14)
        {
            nEdges = nEdges + 1;
        }
    }

    Nwall_ = nEdges;
    for (int i = 0; i < nEdges; i++)
    {
        boundary_first_wall_x_e_.push_back(0);
        boundary_first_wall_y_e_.push_back(0);
    }

    vector<int> Rows(nEdges);
    vector<int> Visited(nEdges);
    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        Visited[i] = 0;
    }
    Visited[0] = 1;

    vector<vector<int>> Edges;
    for (int i = 0; i < nEdges; i++)
    {
        Edges.push_back(vector<int>());
        for (int j = 0; j < 2; j++)
        {
            Edges[i].push_back(0);
        }
    }

    nEdges = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] >= 1 && fcLbl_[i][0] <= 14)
        {
            Rows[nEdges] = i;
            Edges[nEdges][0] = (int)fcVx_[i][0];
            Edges[nEdges][1] = (int)fcVx_[i][1];

            nEdges = nEdges + 1;
        }
    }

    int currentVertex;
    currentVertex = Edges[0][0];

    int nextVertex;

    vector<int> sequence(nEdges + 1);
    sequence[0] = Edges[0][0];

    for (int i = 1; i <= nEdges; i = i + 1)
    {
        if (i == nEdges)
        {
            Visited[0] = 0;
        }

        int idx = 0;
        for (int j = 0; j <= nEdges - 1; j = j + 1)
        {
            if ((Edges[j][0] == currentVertex || Edges[j][1] == currentVertex) && Visited[j] == 0)
            {
                idx = j;
                continue;
            }
        }

        nextVertex = Edges[idx][0] + Edges[idx][1] - currentVertex;
        currentVertex = nextVertex;
        sequence[i] = nextVertex;
        Visited[idx] = 1;
    }

    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        boundary_first_wall_x_e_[i] = vxX_[sequence[i] - 1][0];
        boundary_first_wall_y_e_[i] = vxY_[sequence[i] - 1][0];
    }

    string path1 = Outputpath_ + "/mesh//boundary_first_wall_x_e.txt";
    string path2 = Outputpath_ + "/mesh//boundary_first_wall_y_e.txt";
    std::ofstream fout1(path1, std::ios::out);
    std::ofstream fout2(path2, std::ios::out);
    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        fout1 << boundary_first_wall_x_e_[i] << endl;
        fout2 << boundary_first_wall_y_e_[i] << endl;
    }
    fout1.close();
    fout2.close();
}

void Grid_extern::get_core_cv()
{
    int nEdges = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] == -21 || fcLbl_[i][0] == -25)
        {
            nEdges = nEdges + 1;
        }
    }
    for (int i = 0; i < nEdges; i++)
    {
        core_cv_.push_back(0);
    }
    vector<vector<int>> core_cv_tmp;

    for (int i = 0; i < nEdges; i++)
    {
        core_cv_tmp.push_back(vector<int>());
        for (int j = 2; j < 2; j++)
        {
            core_cv_tmp[i].push_back(0);
        }
    }

    int num = 0;
    for (int i = 0; i <= nFc_ - 1; i = i + 1)
    {
        if (fcLbl_[i][0] == -21 || fcLbl_[i][0] == -25)
        {
            core_cv_tmp[num][0] = (int)fcCv_[i][0];
            core_cv_tmp[num][1] = (int)fcCv_[i][1];

            if (cvVxP_[core_cv_tmp[num][0]][1] > 2)
            {
                core_cv_[num] = core_cv_tmp[num][0];
            }
            else
            {
                core_cv_[num] = core_cv_tmp[num][1];
            }
            num = num + 1;
        }
    }

    string path1 = Outputpath_ + "/mesh//core_cv_tmp.txt";
    string path2 = Outputpath_ + "/mesh//core_cv.txt";
    std::ofstream fout1(path1, std::ios::out);
    std::ofstream fout2(path2, std::ios::out);
    for (int i = 0; i <= nEdges - 1; i = i + 1)
    {
        fout1 << core_cv_tmp[i][0] << "\t" << core_cv_tmp[i][1] << endl;
        fout2 << core_cv_[i] << endl;
    }
    fout1.close();
    fout2.close();
}

void Grid_extern::get_i2j2_cv()
{
    for (int i = 0; i < nx_; i++)
    {
        i2j2_cv_.push_back(vector<vector<int>>());
        for (int j = 0; j < ny_; j++)
        {
            i2j2_cv_[i].push_back(vector<int>());
        }
    }

    i2j2_cv_num_.resize(nx_);
    for (int i = 0; i <= nx_ - 1; i = i + 1)
    {
        i2j2_cv_num_[i].resize(ny_);
    }

    int n1 = 4; // nx_ponit
    for (int i = 0; i <= nx_ - 1; i = i + 1)
    {
        for (int j = 0; j <= ny_ - 1; j = j + 1)
        {

            int num = 0;
            vector<double> polygon1X(4);
            vector<double> polygon1Y(4);
            polygon1X[0] = nx_point_[i];
            polygon1X[1] = nx_point_[i];
            polygon1X[3] = nx_point_[i + 1];
            polygon1X[2] = nx_point_[i + 1];
            polygon1Y[0] = ny_point_[j];
            polygon1Y[1] = ny_point_[j + 1];
            polygon1Y[3] = ny_point_[j];
            polygon1Y[2] = ny_point_[j + 1];
            for (int iCv = 0; iCv <= nCv_ - 1; iCv = iCv + 1)
            {

                if (i == 2 && j == 48 && iCv == 4140)
                {
                    int a = 1;
                }
                int n2 = cvVxP_[iCv][1];
                if (n2 > 2)
                {
                    if (isPolygonOverlap(polygon1X, polygon1Y, n1, cv_VxX_[iCv], cv_VxY_[iCv], n2))
                    {
                        num = num + 1;
                    }
                }
            }

            i2j2_cv_num_[i][j] = num;
            i2j2_cv_[i][j].resize(num);

            num = 0;
            for (int iCv = 0; iCv <= nCv_ - 1; iCv = iCv + 1)
            {
                int n2 = cvVxP_[iCv][1];
                if (n2 > 2)
                {
                    if (isPolygonOverlap(polygon1X, polygon1Y, n1, cv_VxX_[iCv], cv_VxY_[iCv], n2))
                    {
                        i2j2_cv_[i][j][num] = iCv;
                        num = num + 1;
                    }
                }
            }
        }
    }

    string path1 = Outputpath_ + "/mesh//i2j2_cv.txt";
    ofstream fout1(path1, std::ios::out);

    for (int i = 0; i <= nx_ - 1; i = i + 1)
    {
        for (int j = 0; j <= ny_ - 1; j = j + 1)
        {
            for (int n = 0; n <= i2j2_cv_num_[i][j] - 1; n = n + 1)
            {
                fout1 << i << "\t" << j << "\t" << n << "\t" << i2j2_cv_[i][j][n] << endl;
            }
        }
    }
    fout1.close();
}

void Grid_extern::get_max_min_of_Vx(double *max_vxX, double *max_vxY, double *min_vxX, double *min_vxY)
{
    int nVx = nVx_; //?
    double max = -DBL_MAX;
    double min = DBL_MAX;

    for (int i = 0; i <= nVx - 1; i = i + 1)
    {
        double value = vxX_[i][0];

        if (value > max)
        {
            max = value;
        }
        if (value < min)
        {
            min = value;
        }
    }
    *max_vxX = max;
    *min_vxX = min;

    max = -DBL_MAX;
    min = DBL_MAX;

    for (int i = 0; i <= nVx - 1; i = i + 1)
    {
        double value = vxY_[i][0];

        if (value > max)
        {
            max = value;
        }
        if (value < min)
        {
            min = value;
        }
    }
    *max_vxY = max;
    *min_vxY = min;
}

void read_2D_array(vector<vector<double>> &p, int Nt, int Nr, string a)
{
    std::ifstream fp;         // file stream
    fp.open(a, std::ios::in); // ios::in means read
    if (!fp.is_open())
    {
        std::cerr << a + " is not exist, it will be a all zero array" << endl;
        for (int i = 0; i < Nt; i++)
        {
            p.push_back(vector<double>());
            for (int j = 0; j < Nr; j++)
            {
                p[i].push_back(0);
            }
        }
        return;
    }
    for (int i = 0; i < Nt; i++)
    {
        p.push_back(vector<double>());
        for (int j = 0; j < Nr; j++)
        {
            p[i].push_back(0);
        }
    }

    for (int i = 0; i <= Nt - 1; i = i + 1)
    {
        for (int j = 0; j <= Nr - 1; j = j + 1)
        {
            fp >> p[i][j];
        }
    }

    /*for (int i = 0; i <= Nt - 1; i = i + 1)
    {
        for (int j = 0; j <= Nr - 1; j = j + 1)
        {
            printf("%.3lf	", *(*(*p + i) + j));
        }
        printf("\n");
    }*/

    fp.close();
}
void copyFile(std::string source_filename, std::string target_filename)
{
    std::ifstream fin(source_filename, std::ios::in);
    std::ofstream fout(target_filename, std::ios::out);

    if (fin && fout)
    {
        std::string line;
        // 逐行读取源文件的内容并写入目标文件
        while (getline(fin, line))
        {
            fout << line << std::endl;
        }
    }
    else
    {
        if (!fin)
        {
            std::cerr << "in error;" << endl;
        }
        else if (!fout)
        {
            std::cerr << "out error;" << endl;
        }
    }

    fin.close();
    fout.close();
}

void fprint_2d(int Nt, int Nr, vector<vector<double>> &value, string path)
{
    std::ofstream fout(path, std::ios::out);
    for (int i = 0; i <= Nt - 1; i = i + 1)
    {
        for (int j = 0; j <= Nr - 1; j = j + 1)
        {
            fout << value[i][j] << "\t";
        }
        fout << endl;
    }
    fout.close();
}

// extend_mesh
int isPolygonOverlap(vector<double> &polygon1X, vector<double> &polygon1Y, int n1,
                     vector<double> &polygon2X, vector<double> &polygon2Y, int n2)
{
    // �������1�Ķ����Ƿ��ڶ����2��
    for (int i = 0; i < n1; i++)
    {
        if (pointInPolygon(n2, polygon2X, polygon2Y, polygon1X[i], polygon1Y[i]))
        {
            return 1;
        }
    }

    // �������2�Ķ����Ƿ��ڶ����1��
    for (int i = 0; i < n2; i++)
    {
        if (pointInPolygon(n1, polygon1X, polygon1Y, polygon2X[i], polygon2Y[i]))
        {
            return 1;
        }
    }

    // ������εı��Ƿ��ཻ
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            if (isIntersect(polygon1X[i], polygon1Y[i], polygon1X[(i + 1) % n1], polygon1Y[(i + 1) % n1],
                            polygon2X[j], polygon2Y[j], polygon2X[(j + 1) % n2], polygon2Y[(j + 1) % n2]))
            {
                return 1;
            }
        }
    }

    return 0;
}
int pointInPolygon(int polySides, vector<double> &polyX, vector<double> &polyY, double x, double y)
{
    int i, j = polySides - 1;
    int oddNodes = 0;

    for (i = 0; i < polySides; i++)
    {
        if ((polyY[i] < y && polyY[j] >= y || polyY[j] < y && polyY[i] >= y) && (polyX[i] <= x || polyX[j] <= x))
        {
            if (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x)
            {
                oddNodes = !oddNodes;
            }
        }
        j = i;
    }

    return oddNodes;
}
int isIntersect(double x1, double y1, double x2, double y2,
                double x3, double y3, double x4, double y4)
{
    double cross1 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    double cross2 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
    double cross3 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
    double cross4 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);

    return (cross1 * cross2 < 0) && (cross3 * cross4 < 0);
}

int Grid_extern::nx() { return nx_; }
int Grid_extern::ny() { return ny_; }
int Grid_extern::nCv() { return nCv_; }
int Grid_extern::nVx() { return nVx_; }
int Grid_extern::nFc() { return nFc_; }
int Grid_extern::nCmxVx() { return nCmxVx_; }
int Grid_extern::sp_Ntotal() { return sp_Ntotal_; }

string Grid_extern::inputpath_mesh()
{
    return inputpath_mesh_;
}

string Grid_extern::Casepath()
{
    return Casepath_;
}

string Grid_extern::Outputpath()
{
    return Outputpath_;
}

void Grid_extern::Set_nx(int nx) { nx_ = nx; }
void Grid_extern::Set_ny(int ny) { ny_ = ny; }
void Grid_extern::Set_nCv(int nCv) { nCv_ = nCv; }
void Grid_extern::Set_nVx(int nVx) { nVx_ = nVx_; }
void Grid_extern::Set_nFc(int nFc) { nFc_ = nFc; }
void Grid_extern::Set_nCmxVx(int nCmxVx) { nCmxVx_ = nCmxVx; }

double Grid_extern::cvBb(int i, int j)
{
    return cvBb_[i][j];
}
double Grid_extern::cvEb(int i, int j)
{
    return cvEb_[i][j];
}

void Grid_extern::Set_cvBb(int i, int j, double Var)
{
    cvBb_[i][j] = Var;
}
void Grid_extern::Set_cvEb(int i, int j, double Var)
{
    cvEb_[i][j] = Var;
}

namespace eirene
{

    void TriMesh::require(bool cond, const std::string &msg)
    {
        if (!cond)
            throw std::runtime_error(msg);
    }

    // ========== 公有接口 ==========
    void TriMesh::load_eirene_mesh(const std::string &f33, const std::string &f34, const std::string &f35, const std::string &wall)
    {
        read_fort33(f33);
        read_fort34(f34);
        read_fort35(f35);
        compute_centroids_areas();
        Wall_.Init(wall);
    }

    void TriMesh::read_b2fgmtry(const std::string &path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        int nx = -1, ny = -1;
        std::string line;
        std::streampos after_dims;
        while (std::getline(in, line))
        {
            std::istringstream iss(line);
            int a, b;
            if (iss >> a >> b)
            {
                nx = a;
                ny = b;
                after_dims = in.tellg();
                break;
            }
        }
        if (nx <= 0 || ny <= 0)
            throw std::runtime_error("b2fgmtry: failed to read nx, ny");
        const int N = nx * ny;
        std::vector<double> Rcc;
        Rcc.reserve(N);
        std::vector<double> Zcc;
        Zcc.reserve(N);
        in.clear();
        if (after_dims != std::streampos(-1))
            in.seekg(after_dims);
        auto read_N = [&](int Nvals, std::vector<double> &out)
        { double v; while((int)out.size()<Nvals && (in>>v)) out.push_back(v); if((int)out.size()!=Nvals) throw std::runtime_error("b2fgmtry: not enough doubles"); };
        bool foundR = false, foundZ = false;
        while (std::getline(in, line))
        {
            std::string U = line;
            std::transform(U.begin(), U.end(), U.begin(), ::toupper);
            if (!foundR && U.find("R") != std::string::npos && U.find("COORDINATES") != std::string::npos)
            {
                read_N(N, Rcc);
                foundR = true;
            }
            if (!foundZ && U.find("Z") != std::string::npos && U.find("COORDINATES") != std::string::npos)
            {
                read_N(N, Zcc);
                foundZ = true;
            }
            if (foundR && foundZ)
                break;
        }
        if (!foundR || !foundZ)
            throw std::runtime_error("b2fgmtry: failed to find R/Z COORDINATES blocks");
        set_b2_geometry(nx, ny, std::move(Rcc), std::move(Zcc));
    }

    void TriMesh::set_b2_geometry(int nx, int ny, std::vector<double> Rcc, std::vector<double> Zcc)
    {
        if ((int)Rcc.size() != nx * ny || (int)Zcc.size() != nx * ny)
            throw std::runtime_error("B2 geometry size mismatch");
        b2g_.nx = nx;
        b2g_.ny = ny;
        b2g_.Rcc = std::move(Rcc);
        b2g_.Zcc = std::move(Zcc);
        build_b2_bins();
    }
    void TriMesh::read_fort44(const std::string &path) { readFort44_Heuristic(path); }

    std::vector<double> TriMesh::map_fort44_to_triangles(const std::string &field_key, int k_neighbors, double p) const
    {
        auto it = f44_.fields.find(field_key);
        if (it == f44_.fields.end())
            throw std::runtime_error("Fort.44 field not found: " + field_key);
        if (b2g_.nx == 0 || b2g_.ny == 0)
            throw std::runtime_error("B2 geometry not set. Call read_b2fgmtry/set_b2_geometry first.");

        const auto &fld = it->second;
        if ((int)fld.size() != f44_.nx * f44_.ny)
            throw std::runtime_error("Field size != nx*ny.");

        std::vector<double> out(tris_.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t t = 0; t < tris_.size(); ++t)
        {
            const double R = triRc_[t], Z = triZc_[t];
            auto neigh = query_b2_neighbors(R, Z, k_neighbors);
            if (neigh.empty())
                continue;
            if (k_neighbors <= 1 || neigh.size() == 1)
            {
                auto nn = neigh[0];
                int ix = nn.ix, iy = nn.iy;
                out[t] = fld[(size_t)iy * f44_.nx + (size_t)ix];
            }
            else
            {
                double wsum = 0.0, vsum = 0.0;
                for (const auto &q : neigh)
                {
                    double d = std::sqrt(std::max(1e-18, q.dist2));
                    double w = 1.0 / std::pow(d, p);
                    double v = fld[(size_t)q.iy * f44_.nx + (size_t)q.ix];
                    wsum += w;
                    vsum += w * v;
                }
                out[t] = vsum / std::max(1e-300, wsum);
            }
        }
        return out;
    }
    // ========== 私有：fort.33/34/35 ==========
    void TriMesh::read_fort33(const std::string &path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        in >> num_Nodes_;
        // std::cout << num_nodes;
        for (int i = 0; i < num_Nodes_; i++)
        {
            nodes_.push_back(Node());
            in >> nodes_[i].r;
        }
        for (int i = 0; i < num_Nodes_; i++)
        {
            in >> nodes_[i].z;
        }
        for (int i = 0; i < num_Nodes_; i++)
        {
            nodes_[i].r /= 100.;
            nodes_[i].z /= 100.;
        }
        /*std::ofstream out_temp("doc/fort.33");
        for (int i = 0; i < nodes_.size(); i++)
        {
            out_temp << nodes_[i].r << "\t" << nodes_[i].z << endl;
        }
        out_temp.close();*/
    }

    void TriMesh::read_fort34(const std::string &path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        int val_temp;

        in >> num_tris_;
        int V_temp[3];
        // std::cout << num_tris;
        for (int i = 0; i < num_tris_; i++)
        {
            tris_.push_back(Tri());
            // in >> val_temp >> tris_[i].v[0] >> tris_[i].v[1] >> tris_[i].v[2];
            in >> val_temp >> V_temp[0] >> V_temp[1] >> V_temp[2];
            for (int j = 0; j < 3; j++)
            {
                tris_[i].v[j] = V_temp[j] - 1;
            }
        }
        /*std::ofstream out_temp("doc/fort.34");
        for (int i = 0; i < tris_.size(); i++)
        {
            out_temp << tris_[i].v[0] << "\t" << tris_[i].v[1] << "\t" << tris_[i].v[2] << "\t" << endl;
        }
        out_temp.close();*/
    }
    void TriMesh::read_fort35(const std::string &path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        int num_temp0, num_temp[3];
        in >> num_temp0;

        for (int i = 0; i < num_tris_; i++)
        {
            in >> num_temp0;

            for (int m = 0; m < 3; m++)
            {
                in >> num_temp[0] >> num_temp[1] >> num_temp[2];
                for (int n = 0; n < 3; n++)
                    tris_[i].neigh[m * 3 + n] = num_temp[n] - 1;
                tris_[i].neigh[m * 3 + 2] += 1;
            }

            in >> num_temp[0] >> num_temp[1];
            /// Something is wierd.
            if (num_temp[0] > 74)
                tris_[i].neigh[9] = num_temp[0] - 2; // poloidal index
            else if (num_temp[0] > 25)
                tris_[i].neigh[9] = num_temp[0] - 1; // poloidal index
            else
                tris_[i].neigh[9] = num_temp[0];
            tris_[i].neigh[10] = num_temp[1]; // radial index
        }
        /*std::ofstream out_temp("doc/fort.35");
        for (int i = 0; i < tris_.size(); i++)
        {
            for (int m = 0; m < 11; m++)
            {
                out_temp << tris_[i].neigh[m] << "\t";
            }
            out_temp << endl;
        }
        out_temp.close();*/
    }

    void TriMesh::compute_centroids_areas()
    {
        triArea_.resize(tris_.size());
        triVolume_.resize(tris_.size());
        triRc_.resize(tris_.size());
        triZc_.resize(tris_.size());
        for (size_t t = 0; t < tris_.size(); ++t)
        {
            const Node &A = nodes_[tris_[t].v[0]], &B = nodes_[tris_[t].v[1]], &C = nodes_[tris_[t].v[2]];
            double area = 0.5 * std::abs((B.r - A.r) * (C.z - A.z) - (B.z - A.z) * (C.r - A.r));
            triArea_[t] = area;
            triRc_[t] = (A.r + B.r + C.r) / 3.0;
            triZc_[t] = (A.z + B.z + C.z) / 3.0;
            triVolume_[t] = 2.0 * M_PI * triRc_[t] * triArea_[t];
        }
        /*std::ofstream out_temp("doc/Triz_Info.txt");
        for (int i = 0; i < tris_.size(); i++)
        {
            out_temp << i << "\t" << triArea_[i] << "\t" << triRc_[i] << "\t" << triZc_[i] << endl;
        }
        out_temp.close();*/
    }

    // ========== 私有：fort.44 ==========
    bool TriMesh::find_token_value(const std::string &s, const std::string &key, int &out)
    {
        auto p = s.find(key);
        if (p == std::string::npos)
            return false;
        p = s.find_first_of("=:", p);
        if (p == std::string::npos)
            return false;
        std::istringstream iss(s.substr(p + 1));
        int v;
        if (iss >> v)
        {
            out = v;
            return true;
        }
        return false;
    }

    void TriMesh::readFort44_Heuristic(const std::string &path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        {
            std::string line;
            int hits = 0;
            std::streampos pos = in.tellg();
            for (int i = 0; i < 100 && std::getline(in, line); ++i)
            {
                hits += (find_token_value(line, "NDX", f44_.nx) ? 1 : 0);
                hits += (find_token_value(line, "NDY", f44_.ny) ? 1 : 0);
                hits += (find_token_value(line, "NATM", f44_.natm) ? 1 : 0);
                hits += (find_token_value(line, "NMOL", f44_.nmol) ? 1 : 0);
                if (hits >= 2 && f44_.nx > 0 && f44_.ny > 0)
                    break;
            }
            in.clear();
            in.seekg(pos);
        }
        if (f44_.nx <= 0 || f44_.ny <= 0)
            throw std::runtime_error("Fort.44: failed to detect NDX/NDY from header");
        auto read2D = [&](const std::string &key)
        { const size_t n=(size_t)f44_.nx*(size_t)f44_.ny; std::vector<double> v; v.reserve(n); double x; while(v.size()<n && (in>>x)) v.push_back(x); if(v.size()<n) throw std::runtime_error("Fort.44: not enough data for "+key); f44_.fields[key]=std::move(v); };
        read2D("dab2");
        read2D("tab2");
        if (f44_.nmol > 0)
        {
            read2D("dmb2");
            read2D("tmb2");
        }
        read2D("rfluxa");
        read2D("pfluxa");
        if (f44_.nmol > 0)
        {
            read2D("rfluxm");
            read2D("pfluxm");
        }
    }

    // ========== 私有：B2 近邻桶 ==========
    void TriMesh::build_b2_bins()
    {
        if (b2g_.nx == 0 || b2g_.ny == 0)
            return;
        Rmin = Zmin = std::numeric_limits<double>::infinity();
        Rmax = Zmax = -std::numeric_limits<double>::infinity();
        const size_t N = (size_t)b2g_.nx * b2g_.ny;
        for (size_t i = 0; i < N; ++i)
        {
            Rmin = std::min(Rmin, b2g_.Rcc[i]);
            Rmax = std::max(Rmax, b2g_.Rcc[i]);
            Zmin = std::min(Zmin, b2g_.Zcc[i]);
            Zmax = std::max(Zmax, b2g_.Zcc[i]);
        }
        int nb = std::max(1, (int)std::sqrt((double)N));
        nbx = std::max(1, nb);
        nbz = std::max(1, nb / 2);
        dRb = (Rmax - Rmin) / std::max(1, nbx);
        dZb = (Zmax - Zmin) / std::max(1, nbz);
        bins.assign((size_t)nbx * nbz, {});
        for (int iy = 0; iy < b2g_.ny; ++iy)
        {
            for (int ix = 0; ix < b2g_.nx; ++ix)
            {
                size_t id = (size_t)iy * b2g_.nx + (size_t)ix;
                auto pr = bin_of(b2g_.Rcc[id], b2g_.Zcc[id]);
                int bx = pr.first, bz = pr.second;
                bins[(size_t)bz * nbx + (size_t)bx].push_back((int)id);
            }
        }
    }

    std::pair<int, int> TriMesh::bin_of(double R, double Z) const
    {
        int bx = (dRb > 0) ? (int)std::floor((R - Rmin) / dRb) : 0;
        int bz = (dZb > 0) ? (int)std::floor((Z - Zmin) / dZb) : 0;
        bx = std::clamp(bx, 0, nbx - 1);
        bz = std::clamp(bz, 0, nbz - 1);
        return {bx, bz};
    }

    std::vector<TriMesh::NN> TriMesh::query_b2_neighbors(double R, double Z, int k) const
    {
        if (k <= 0)
            k = 1;
        k = std::min(k, b2g_.nx * b2g_.ny);
        auto center = bin_of(R, Z);
        int bx = center.first, bz = center.second;
        std::vector<NN> cand;
        cand.reserve((size_t)k * 4);
        int rad = 0, maxrad = std::max(nbx, nbz);
        while ((int)cand.size() < k && rad < maxrad)
        {
            for (int dz = -rad; dz <= rad; ++dz)
            {
                for (int dx = -rad; dx <= rad; ++dx)
                {
                    if (std::max(std::abs(dx), std::abs(dz)) != rad)
                        continue;
                    int cx = bx + dx, cz = bz + dz;
                    if (cx < 0 || cz < 0 || cx >= nbx || cz >= nbz)
                        continue;
                    const auto &vec = bins[(size_t)cz * nbx + (size_t)cx];
                    for (int id : vec)
                    {
                        int ix = id % b2g_.nx;
                        int iy = id / b2g_.nx;
                        double dR = b2g_.Rcc[(size_t)id] - R;
                        double dZ = b2g_.Zcc[(size_t)id] - Z;
                        cand.push_back({ix, iy, dR * dR + dZ * dZ});
                    }
                }
            }
            ++rad;
        }
        if (cand.empty())
            return cand;
        std::nth_element(cand.begin(), cand.begin() + std::min((int)cand.size(), k) - 1, cand.end(), [](const NN &a, const NN &b)
                         { return a.dist2 < b.dist2; });
        cand.resize(std::min((int)cand.size(), k));
        std::sort(cand.begin(), cand.end(), [](const NN &a, const NN &b)
                  { return a.dist2 < b.dist2; });
        return cand;
    }

    int TriMesh::num_Nodes()
    {
        return num_Nodes_;
    }

    int TriMesh::num_tris()
    {
        return num_tris_;
    }

    /// @brief // find the target in the triangle
    /// @param N_poloidal number of the poloidal grid
    /// @param N_radial number of the radial grid
    void TriMesh::target_find(int N_poloidal, int N_radial)
    {
        N_poloidal_ = N_poloidal;
        N_radial_ = N_radial;
        Num_Target_ = N_radial * 2;

        TargetIndex_.resize(Num_Target_);
        Cos_Target_.resize(Num_Target_);
        Sin_Target_.resize(Num_Target_);
        Vol_Target_.resize(Num_Target_);

        Target_.resize(Num_Target_);
        Mid_Target_.resize(Num_Target_);

        for (int i = 0; i < TargetIndex_.size(); i++)
        {
            TargetIndex_[i].push_back(-1);
            TargetIndex_[i].push_back(-1);
            Mid_Target_[i].resize(2);
            Target_[i].resize(4);
        }
        for (int i = 0; i < num_tris_; i++)
        {
            if (tris_[i].neigh[2] != 0)
            {
                if (tris_[i].neigh[9] == 1)
                {
                    TargetIndex_[tris_[i].neigh[10]][0] = i;
                    TargetIndex_[tris_[i].neigh[10]][1] = 0;
                }
                else if (tris_[i].neigh[9] == N_poloidal_ - 2)
                {
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][0] = i;
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][1] = 0;
                }
            }
            else if (tris_[i].neigh[5] != 0)
            {
                if (tris_[i].neigh[9] == 1)
                {
                    TargetIndex_[tris_[i].neigh[10]][0] = i;
                    TargetIndex_[tris_[i].neigh[10]][1] = 1;
                }
                else if (tris_[i].neigh[9] == N_poloidal_ - 2)
                {
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][0] = i;
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][1] = 1;
                }
            }
            else if (tris_[i].neigh[8] != 0)
            {
                if (tris_[i].neigh[9] == 1)
                {
                    TargetIndex_[tris_[i].neigh[10]][0] = i;
                    TargetIndex_[tris_[i].neigh[10]][1] = 2;
                }
                else if (tris_[i].neigh[9] == N_poloidal_ - 2)
                {
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][0] = i;
                    TargetIndex_[tris_[i].neigh[10] + N_radial_][1] = 2;
                }
            }
        }
        for (int i = 0; i < Num_Target_; i++)
        {
            if (TargetIndex_[i][0] != -1)
            {
                Target_[i][0] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1]]].r;
                Target_[i][1] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1]]].z;
                if (TargetIndex_[i][1] == 0 || TargetIndex_[i][1] == 1)
                {
                    Target_[i][2] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1] + 1]].r;
                    Target_[i][3] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1] + 1]].z;
                }
                else if (TargetIndex_[i][1] == 2)
                {
                    Target_[i][2] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1] - 2]].r;
                    Target_[i][3] = nodes_[tris_[TargetIndex_[i][0]].v[TargetIndex_[i][1] - 2]].z;
                }
                else
                {
                    std::cerr << "Error in targer_find()!" << endl;
                }
            }
        }
        double deltaX, deltaY, deltaL;
        for (int i = 0; i < Num_Target_; i++)
        {
            Mid_Target_[i][0] = 0.5 * (Target_[i][0] + Target_[i][2]);
            Mid_Target_[i][1] = 0.5 * (Target_[i][1] + Target_[i][3]);

            double deltaX = Target_[i][2] - Target_[i][0];
            double deltaY = Target_[i][3] - Target_[i][1];
            double deltaL = std::sqrt(deltaX * deltaX + deltaY * deltaY);
            if (deltaL > 0) // 防御性检查，防止除以 0
            {
                Cos_Target_[i] = deltaX / deltaL;
                Sin_Target_[i] = deltaY / deltaL;
            }
            else
            {
                // 异常网格处理：如果两个点完全重合，给一个默认的切向
                Cos_Target_[i] = 0.0;
                Sin_Target_[i] = 0.0;
            }
        }
        /*ofstream out_temp3("doc/Grid4_target_cos.txt");
        for (int i = 0; i < N_radial_ * 2; i++)
        {
            if (TargetIndex_[i][0] != -1)
                out_temp3 << i << "\t" << Cos_trimesh_[TargetIndex_[i][0]][TargetIndex_[i][1]] << "\t" << Cos_Target_[i] << "\t" << Sin_trimesh_[TargetIndex_[i][0]][TargetIndex_[i][1]] << "\t" << Sin_Target_[i] << endl;
        }
        out_temp3.close();*/
        /*ofstream out_temp3("doc/Grid4_target_1.txt");
        for (int i = 0; i < N_radial_ * 2; i++)
        {
            out_temp3 << i << "\t" << Target_[i][0] << "\t" << Target_[i][1] << "\t" << Target_[i][2] << "\t" << Target_[i][3] << "\t" << Cos_Target_[i] << "\t" << Sin_Target_[i] << endl;
        }
        out_temp3.close();*/
        /*ofstream out_temp3("doc/angle_target1.txt");
        for (int i = 0; i < N_radial_ * 2; i++)
        {
            out_temp3 << i << "\t" << Cos_Target_[i] << "\t" << Sin_Target_[i] << endl;
        }
        out_temp3.close();*/
        /*ofstream out_temp1("doc/tri_target_index.txt");
        ofstream out_temp2("doc/tri_target.txt");
        for (int i = 0; i < N_radial_ * 2; i++)
        {
            out_temp1 << i << "\t" << TargetIndex_[i][0] << "\t" << TargetIndex_[i][1] << endl;
            out_temp2 << i << "\t" << Target_[i][0] << "\t" << Target_[i][1] << endl;
        }
        out_temp1.close();
        out_temp2.close();*/
    }

    ///@brief find the b2 index;
    void TriMesh::b2_find()
    {
        b2_index_.resize(num_tris_);
        for (int i = 0; i < num_tris_; i++)
        {
            if (if_in_plasmagrid(i))
            {
                b2_index_[i].push_back(tris_[i].neigh[9]);
                b2_index_[i].push_back(tris_[i].neigh[10]);
            }
            else
            {
                for (int j = 0; j < 2; j++)
                    b2_index_[i].push_back(-1);
            }
        }
        /*ofstream out_temp("doc/b2_index.txt");
        for (int i = 0; i < num_tris_; i++)
        {
            out_temp << i << "\t" << b2_index_[i][0] << "\t" << b2_index_[i][1] << endl;
        }
        out_temp.close();*/
    }

    int TriMesh::Num_Target() { return Num_Target_; }

    double TriMesh::Target(int i, int j)
    {
        return Target_[i][j];
    }

    int TriMesh::targetIndex(int i, int j)
    {
        return TargetIndex_[i][j];
    }

    /// @brief // process the information in the triangle
    /// @param N_poloidal number of the poloidal grid
    /// @param N_radial number of the radial grid
    void TriMesh::mesh_find(int N_poloidal, int N_radial)
    {
        b2_find();
        lines_find();
        target_find(N_poloidal, N_radial);
        build_wall_candidates_by_tri();
    }

    int TriMesh::b2_index(int i, int j)
    {
        return b2_index_[i][j];
    }

    /// @brief find the line information
    void TriMesh::lines_find()
    {
        lines_info_.resize(num_tris_);
        for (int i = 0; i < num_tris_; i++)
        {
            lines_info_[i].resize(3);
            for (int j = 0; j < 3; j++)
            {
                lines_info_[i][j] = 0;
            }
        }

        for (int i = 0; i < num_tris_; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (tris_[i].neigh[j * 3 + 2] != 0)
                {
                    if (tris_[i].neigh[j * 3 + 2] == -1)
                    {
                        lines_info_[i][j] = -1;
                    }
                    else if (tris_[i].neigh[j * 3 + 2] == 3)
                    {
                        lines_info_[i][j] = 3;
                    }
                    else if (tris_[i].neigh[j * 3 + 2] == 4)
                    {
                        lines_info_[i][j] = 4;
                    }
                    else if (tris_[i].neigh[j * 3 + 2] == -2)
                    {
                        lines_info_[i][j] = -2;
                    }
                    else if (tris_[i].neigh[j * 3 + 2] == 1)
                    {
                        lines_info_[i][j] = 1;
                    }
                    else
                    {
                        std::cerr << "lines information find function error!" << endl;
                        std::cerr << "lines information " << i << ", " << j << " is: " << tris_[i].neigh[j * 3 + 2] << endl;
                        std::cerr << nodes_[tris_[i].v[0]].r << ", " << nodes_[tris_[i].v[0]].z << endl;
                        pause();
                    }
                }
                else
                {
                    if (if_in_plasmagrid(i) == 1)
                    {
                        if (if_in_plasmagrid(tris_[i].neigh[3 * j]) != 1)
                            lines_info_[i][j] = 2;
                    }
                    if (if_in_plasmagrid(i) != 1)
                    {
                        if (if_in_plasmagrid(tris_[i].neigh[3 * j]) == 1)
                            lines_info_[i][j] = 2;
                    }
                }
            }
        }
        Cos_trimesh_.resize(num_tris_);
        Sin_trimesh_.resize(num_tris_);
        double deltaX, deltaY, deltaL, Cos_temp, Sin_temp;
        for (int i = 0; i < num_tris_; i++)
        {
            deltaX = nodes_[tris_[i].v[1]].r - nodes_[tris_[i].v[0]].r;
            deltaY = nodes_[tris_[i].v[1]].z - nodes_[tris_[i].v[0]].z;
            deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
            Cos_trimesh_[i].push_back(deltaX / deltaL);
            Sin_trimesh_[i].push_back(deltaY / deltaL);

            deltaX = nodes_[tris_[i].v[2]].r - nodes_[tris_[i].v[1]].r;
            deltaY = nodes_[tris_[i].v[2]].z - nodes_[tris_[i].v[1]].z;
            deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
            Cos_trimesh_[i].push_back(deltaX / deltaL);
            Sin_trimesh_[i].push_back(deltaY / deltaL);

            deltaX = nodes_[tris_[i].v[0]].r - nodes_[tris_[i].v[2]].r;
            deltaY = nodes_[tris_[i].v[0]].z - nodes_[tris_[i].v[2]].z;
            deltaL = sqrt(deltaX * deltaX + deltaY * deltaY);
            Cos_trimesh_[i].push_back(deltaX / deltaL);
            Sin_trimesh_[i].push_back(deltaY / deltaL);
        }
        /*ofstream out_temp("doc/lines_info.txt");
        for (int i = 0; i < num_tris_; i++)
        {
            out_temp << i << "\t" << lines_info_[i][0] << "\t" << lines_info_[i][1] << "\t" << lines_info_[i][2] << "\t";
            for (int j = 0; j < 3; j++)
                out_temp << nodes_[tris_[i].v[j]].r << ", " << nodes_[tris_[i].v[j]].z << "\t";
            out_temp << endl;
        }
        out_temp.close();*/
        /*ofstream out_temp("doc/cos_sin.txt");
        for (int i = 0; i < num_tris_; i++)
        {
            out_temp << i << "\t" << Cos_trimesh_[i][0] << "\t" << Cos_trimesh_[i][1] << "\t" << Cos_trimesh_[i][2] << "\t";
            out_temp << Sin_trimesh_[i][0] << "\t" << Sin_trimesh_[i][1] << "\t" << Sin_trimesh_[i][2] << "\t";
            out_temp << endl;
        }
        out_temp.close();*/
    }

    int TriMesh::lines_info(int i, int j)
    {
        return lines_info_[i][j];
    }

    const std::vector<int> &TriMesh::WallCandidatesForTri(int tri) const
    {
        static const std::vector<int> empty;
        if (tri < 0 || tri >= static_cast<int>(wall_candidates_by_tri_.size()))
            return empty;
        return wall_candidates_by_tri_[tri];
    }

    void TriMesh::build_wall_candidates_by_tri()
    {
        wall_candidates_by_tri_.assign(num_tris_, {});
        constexpr double padding = 1.e-8;
        for (int tri = 0; tri < num_tris_; ++tri)
        {
            const Node &a = nodes_[tris_[tri].v[0]];
            const Node &b = nodes_[tris_[tri].v[1]];
            const Node &c = nodes_[tris_[tri].v[2]];
            const double min_x = std::min({a.r, b.r, c.r}) - padding;
            const double max_x = std::max({a.r, b.r, c.r}) + padding;
            const double min_y = std::min({a.z, b.z, c.z}) - padding;
            const double max_y = std::max({a.z, b.z, c.z}) + padding;
            Wall_.CandidateSegments(min_x, max_x, min_y, max_y, wall_candidates_by_tri_[tri]);
        }
    }

    /// @brief if the triangle grid in the plasma region
    /// @param i index of plasma
    /// @return 1: in the plasma region; 2: outside of the plasma region
    int TriMesh::if_in_plasmagrid(int i)
    {
        if (tris_[i].neigh[9] != -1)
            return 1;
        else
            return 0;
    }

    double TriMesh::TritoB2(std::vector<std::vector<double>> &Var, int i)
    {
        int m = tris_[i].neigh[9], n = tris_[i].neigh[10];
        if (m == -2 && n == -1)
            return 0;
        else if (m != -2 && n != -1)
            return Var[m][n];
        else
            throw std::runtime_error("function 1 TritiB2 have error.");
    }
    double TriMesh::TritoB2(double **Var, int i)
    {
        int m = tris_[i].neigh[9], n = tris_[i].neigh[10];
        if (m == -2 && n == -1)
            return 0;
        else if (m != -2 && n != -1)
            return Var[m][n];
        else
            throw std::runtime_error("function 2 TritiB2 have error.");
    }

    void TriMesh::judge_orientation()
    {
        ofstream out_temp("doc/orientation.txt");
        Point P_temp[3];
        double ori;
        for (int i = 0; i < num_tris_; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                P_temp[j].SetX(nodes_[tris_[i].v[j]].r);
                P_temp[j].SetY(nodes_[tris_[i].v[j]].z);
            }
            ori = Tools::orientation(P_temp[0], P_temp[1], P_temp[2]);
            if (ori != 1)
                std::cout << "The direction of " + to_string(i) + " triangle is counterclockwise\n";
            out_temp << i << "\t" << ori << endl;
        }
        out_temp.close();
    }

    WALL::WALL()
    {
    }

    WALL::WALL(string wallpath)
    {
        Wallpath_ = wallpath;
        Wallread();
    }

    WALL::~WALL()
    {
    }

    void WALL::Init(string wallpath)
    {
        Wallpath_ = wallpath;
        Wallread();
    }

    void WALL::Wallread()
    {
        /// Tokamak wall position reading(Clockwise)
        string Path_Wall = Wallpath_;
        std::ifstream fp;
        fp.open(Path_Wall, std::ios::in); // ios::in means read
        if (!fp.is_open())
        {
            std::cerr << "This file READING for wall have some problem!!!\n";
        }
        fp >> Wall_num_;
        Wall_num_ -= 1;
        // std::cout << Wall_num_ << endl;
        double a, b;
        x_.resize(Wall_num_);
        y_.resize(Wall_num_);
        for (int i = 0; i < Wall_num_; i++)
        {
            fp >> x_[i] >> y_[i];
        }
        fp.close();
        x_.push_back(x_[0]);
        y_.push_back(y_[0]);
        segments_.resize(Wall_num_);
        bounds_.resize(Wall_num_);
        for (int i = 0; i < Wall_num_; i++)
        {
            segments_[i] = {x_[i], y_[i], x_[i + 1], y_[i + 1]};
            bounds_[i] = {std::min(segments_[i][0], segments_[i][2]), std::max(segments_[i][0], segments_[i][2]), std::min(segments_[i][1], segments_[i][3]), std::max(segments_[i][1], segments_[i][3])};
            dx_.push_back(x_[i + 1] - x_[i]);
            dy_.push_back(y_[i + 1] - y_[i]);
            Length_sq_.push_back(dx_[i] * dx_[i] + dy_[i] * dy_[i]);
        }
        BuildCandidateBins();
        CalAngleWall();
    }

    int WALL::IfinWall(double x, double y)
    {
        int count = 0;
        double secx, secy;
        for (int i = 0; i < x_.size() - 1; i++)
        {
            if (min(y_[i], y_[i + 1]) < y && max(y_[i], y_[i + 1]) > y)
            {
                if (max(x_[i], x_[i + 1]) < x)
                {
                    count += 1;
                }
                else if (x_[i] < x || x_[i + 1] < x)
                {
                    if (Tools::get_line_intersection(x, y, -1, y, x_[i], y_[i], x_[i + 1], y_[i + 1], &secx, &secy))
                        count += 1;
                }
            }
        }
        if (count % 2 == 1)
            return true;
        else
            return false;
    }

    int WALL::Wall_num()
    {
        return Wall_num_;
    }

    /// @brief position of the wall point
    /// @param num num of the wall point
    /// @param xy 0 for x position, 1 for y position
    /// @return x or y position fo the num_th wall point
    double WALL::P(int num, int xy)
    {
        if (num > Wall_num_ || num < 0)
        {
            std::cout << "Wall index error" << endl;
            exit(0);
        }
        if (xy == 0)
            return x_[num];
        else if (xy == 1)
            return y_[num];
        else
            __throw_runtime_error("Wall xy error!\n");
    }

    const std::array<double, 4> &WALL::Segment(int num) const
    {
        return segments_[num];
    }

    const std::array<double, 4> &WALL::Bounds(int num) const
    {
        return bounds_[num];
    }

    void WALL::BuildCandidateBins()
    {
        Wall_bins_.clear();
        if (bounds_.empty())
            return;

        double min_x = bounds_[0][0];
        double max_x = bounds_[0][1];
        double min_y = bounds_[0][2];
        double max_y = bounds_[0][3];
        for (const auto &bounds : bounds_)
        {
            min_x = std::min(min_x, bounds[0]);
            max_x = std::max(max_x, bounds[1]);
            min_y = std::min(min_y, bounds[2]);
            max_y = std::max(max_y, bounds[3]);
        }

        const int bin_count = std::max(8, static_cast<int>(std::sqrt(static_cast<double>(bounds_.size()))));
        Wall_bin_nx_ = bin_count;
        Wall_bin_ny_ = bin_count;
        Wall_bin_min_x_ = min_x;
        Wall_bin_min_y_ = min_y;
        Wall_bin_dx_ = (max_x - min_x) / Wall_bin_nx_;
        Wall_bin_dy_ = (max_y - min_y) / Wall_bin_ny_;
        if (Wall_bin_dx_ <= 0.0 || Wall_bin_dy_ <= 0.0)
        {
            Wall_bins_.clear();
            return;
        }

        Wall_bins_.assign(Wall_bin_nx_ * Wall_bin_ny_, std::vector<int>());
        for (int i = 0; i < Wall_num_; i++)
        {
            const auto &bounds = bounds_[i];
            int ix0 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((bounds[0] - Wall_bin_min_x_) / Wall_bin_dx_)));
            int ix1 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((bounds[1] - Wall_bin_min_x_) / Wall_bin_dx_)));
            int iy0 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((bounds[2] - Wall_bin_min_y_) / Wall_bin_dy_)));
            int iy1 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((bounds[3] - Wall_bin_min_y_) / Wall_bin_dy_)));
            for (int ix = ix0; ix <= ix1; ix++)
                for (int iy = iy0; iy <= iy1; iy++)
                    Wall_bins_[ix * Wall_bin_ny_ + iy].push_back(i);
        }
    }

    void WALL::CandidateSegments(double min_x, double max_x, double min_y, double max_y, std::vector<int> &candidates) const
    {
        candidates.clear();
        if (Wall_bins_.empty() || Wall_bin_dx_ <= 0.0 || Wall_bin_dy_ <= 0.0)
        {
            candidates.reserve(Wall_num_);
            for (int i = 0; i < Wall_num_; i++)
                candidates.push_back(i);
            return;
        }

        int ix0 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((min_x - Wall_bin_min_x_) / Wall_bin_dx_)));
        int ix1 = std::max(0, std::min(Wall_bin_nx_ - 1, static_cast<int>((max_x - Wall_bin_min_x_) / Wall_bin_dx_)));
        int iy0 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((min_y - Wall_bin_min_y_) / Wall_bin_dy_)));
        int iy1 = std::max(0, std::min(Wall_bin_ny_ - 1, static_cast<int>((max_y - Wall_bin_min_y_) / Wall_bin_dy_)));

        static thread_local std::vector<int> seen;
        static thread_local int stamp = 1;
        if (static_cast<int>(seen.size()) < Wall_num_)
            seen.assign(Wall_num_, 0);
        if (stamp == std::numeric_limits<int>::max())
        {
            std::fill(seen.begin(), seen.end(), 0);
            stamp = 1;
        }
        const int current_stamp = stamp++;

        for (int ix = ix0; ix <= ix1; ix++)
        {
            for (int iy = iy0; iy <= iy1; iy++)
            {
                const auto &bin = Wall_bins_[ix * Wall_bin_ny_ + iy];
                for (int index : bin)
                {
                    if (seen[index] != current_stamp)
                    {
                        seen[index] = current_stamp;
                        candidates.push_back(index);
                    }
                }
            }
        }
    }

    void WALL::CalAngleWall()
    {
        double deltaL;
        for (int i = 0; i < Wall_num(); i++)
        {
            deltaL = sqrt(Length_sq(i));
            Cos_wall_.push_back(-dx(i) / deltaL);
            Sin_wall_.push_back(-dy(i) / deltaL);
        }
    }

    void TriMesh::get_sx()
    {
        double pi = 3.1415926535897;
        double x1, x2, y1, y2, h1, h2, l1, l2, s1, s2;
        double y0;

        // in
        for (int i = 0; i <= Num_Target_ - 1; i = i + 1)
        {
            x1 = Target_[i][0];
            x2 = Target_[i][2];
            y1 = Target_[i][1];
            y2 = Target_[i][3];

            if (x1 == x2)
            {
                Vol_Target_[i] = 2.0 * pi * x1 * fabs(y2 - y1);
                continue;
            }

            if (y1 == y2)
            {
                Vol_Target_[i] = 2.0 * pi * fabs(x1 * x1 - x2 * x2);
                continue;
            }

            get_point_to_y_line(x1, y1, x2, y2, &y0);
            h1 = fabs(y1 - y0);
            h2 = fabs(y2 - y0);
            l1 = sqrt(h1 * h1 + x1 * x1);
            l2 = sqrt(h2 * h2 + x2 * x2);

            s1 = pi * l1 * x1;
            s2 = pi * l2 * x2;

            Vol_Target_[i] = fabs(s1 - s2);
        }
        /*ofstream out("doc/Vol_target.txt");
        for (int i = 0; i < Num_Target_; i++)
            out << i << "\t" << Vol_Target_[i] << endl;
        out.close();*/
    }

    double TriMesh::triVolume(int i)
    {
        return triVolume_[i];
    }
    double WALL::Sin_Wall(int i)
    {
        if (i != -1)
            return Sin_wall_[i];
        else
            throw std::logic_error("wrong Sin_Wall!\n");
    }
    double WALL::Cos_Wall(int i)
    {
        if (i != -1)
            return Cos_wall_[i];
        else
            throw std::logic_error("wrong Cos_Wall!\n");
    }
    double WALL::dx(int i) { return dx_[i]; }
    double WALL::dy(int i) { return dy_[i]; }
    double WALL::Length_sq(int i) { return Length_sq_[i]; }

    /**
     * @brief  终极质量守恒兜底：寻找距离网格碰撞点最近的真实器壁线段，并计算精确投影坐标
     * * @param x_mesh      网格上的粒子碰撞点 X 坐标
     * @param y_mesh        网格上的粒子碰撞点 Y 坐标
     * @param num_walls     真实器壁线段的总数 (wall_x1 数组的长度)
     * @param wall          器壁线段
     * @param nearest_idx   [输出] 距离最近的真实器壁线段的索引 (0 ~ num_walls-1)
     * @param ix            [输出] 粒子最终投射到该真实器壁上的精确 X 坐标
     * @param iy            [输出] 粒子最终投射到该真实器壁上的精确 Y 坐标
     */
    int WALL::findNearestWallFast(double x_mesh, double y_mesh)
    {
        double min_dist_sq = 1e30;
        int nearest_idx = -1;

        int num_walls = Wall_num();

        // 整个循环里没有任何减法算长度，全是乘加运算 (FMA)，CPU 执行极快
        for (int i = 0; i < num_walls; ++i)
        {
            double len_sq = Length_sq(i);
            if (len_sq < 1e-12)
                continue;

            // 利用预存的 dx, dy 算 t
            double t = ((x_mesh - P(i, 0)) * dx(i) + (y_mesh - P(i, 1)) * dy(i)) / len_sq;

            if (t < 0.0)
                t = 0.0;
            if (t > 1.0)
                t = 1.0;

            // 算最近点
            double cx = P(i, 0) + t * dx(i);
            double cy = P(i, 1) + t * dy(i);

            double dist_sq = (x_mesh - cx) * (x_mesh - cx) + (y_mesh - cy) * (y_mesh - cy);

            if (dist_sq < min_dist_sq)
            {
                min_dist_sq = dist_sq;
                nearest_idx = i;
            }
        }
        if (nearest_idx != -1)
            return nearest_idx;
        else
            throw std::logic_error("Wrong answer in findNearestWallFast()");
    }

    double TriMesh::Sin_Target(int i)
    {
        if (i != -1)
            return Sin_Target_[i];
        else
            throw std::logic_error("wrong Sin_Target!\n");
    }
    double TriMesh::Cos_Target(int i)
    {
        if (i != -1)
            return Cos_Target_[i];
        else
            throw std::logic_error("wrong Cos_Target!\n");
    }
    std::pair<double, double> TriMesh::InwardTargetTangent(int i) const
    {
        if (i < 0 || i >= Num_Target_)
            throw std::logic_error("wrong target index for inward tangent\n");
        double target_cos = Cos_Target_[i];
        double target_sin = Sin_Target_[i];
        const int tri_index = TargetIndex_[i][0];
        if (tri_index < 0 || static_cast<std::size_t>(tri_index) >= tris_.size())
            return {target_cos, target_sin};

        const double mid_x = Mid_Target_[i][0];
        const double mid_y = Mid_Target_[i][1];
        const auto &tri = tris_[tri_index];
        const double cx =
            (nodes_[tri.v[0]].r + nodes_[tri.v[1]].r + nodes_[tri.v[2]].r) / 3.;
        const double cy =
            (nodes_[tri.v[0]].z + nodes_[tri.v[1]].z + nodes_[tri.v[2]].z) / 3.;
        const double normal_x = -target_sin;
        const double normal_y = target_cos;
        if ((cx - mid_x) * normal_x + (cy - mid_y) * normal_y < 0.)
        {
            target_cos = -target_cos;
            target_sin = -target_sin;
        }
        return {target_cos, target_sin};
    }
    double TriMesh::Mid_Target(int i, int j) { return Mid_Target_[i][j]; }
    double TriMesh::Vol_Target(int i) { return Vol_Target_[i]; }
    double TriMesh::Sin_trimesh(int i, int j)
    {
        if (i < 0 || i >= static_cast<int>(Sin_trimesh_.size()) ||
            j < 0 || j >= static_cast<int>(Sin_trimesh_[i].size()))
            throw std::out_of_range("Sin_trimesh index out of range: tri=" + std::to_string(i) +
                                    " edge=" + std::to_string(j));
        return Sin_trimesh_[i][j];
    }
    double TriMesh::Cos_trimesh(int i, int j)
    {
        if (i < 0 || i >= static_cast<int>(Cos_trimesh_.size()) ||
            j < 0 || j >= static_cast<int>(Cos_trimesh_[i].size()))
            throw std::out_of_range("Cos_trimesh index out of range: tri=" + std::to_string(i) +
                                    " edge=" + std::to_string(j));
        return Cos_trimesh_[i][j];
    }
    double sign(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
    {
        return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
    }

    bool TriMesh::Ifingrid(int Tri_Index, double px, double py)
    {
        // 2. 获取三个顶点的实际坐标 (A, B, C)
        double x1 = nodes_[tris_[Tri_Index].v[0]].r;
        double y1 = nodes_[tris_[Tri_Index].v[0]].z;
        double x2 = nodes_[tris_[Tri_Index].v[1]].r;
        double y2 = nodes_[tris_[Tri_Index].v[1]].z;
        double x3 = nodes_[tris_[Tri_Index].v[2]].r;
        double y3 = nodes_[tris_[Tri_Index].v[2]].z;

        // 3. 计算点 P 相对于三条边的位置 (叉积符号)
        // 这里的逻辑是：如果点在三角形内，它必须在三条边的“同一侧”
        double d1 = sign(px, py, x1, y1, x2, y2);
        double d2 = sign(px, py, x2, y2, x3, y3);
        double d3 = sign(px, py, x3, y3, x1, y1);

        // 4. 判断逻辑
        bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
        bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

        // 如果同时存在正值和负值，说明点在某些边的左侧，而在另一些边的右侧 -> 也就是在外部
        // 只有当全部为正 或 全部为负时，点才在内部
        // !(has_neg && has_pos) 等价于 (all positive) || (all negative)
        return !(has_neg && has_pos);
    }
    void TriMesh::SetB2Index_find(int i, int m, int n)
    {
        /*b2_index_[i].push_back(m);
        b2_index_[i].push_back(n);*/
        b2_index_[i][0] = m;
        b2_index_[i][1] = n;
    }
}
