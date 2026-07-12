#include "SOLPSField.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "usage: check_d_boundary_flux fnay_Dplus.dat\n";
        return 2;
    }

    try
    {
        const SOLPSIndexedField field =
            ReadSOLPSIndexedField(argv[1], "Radial particle flux");
        int min_i = 1000000;
        int max_i = -1000000;
        int min_j = 1000000;
        int max_j = -1000000;
        double pfr1 = 0.;
        double pfr2 = 0.;
        double outer = 0.;
        for (const auto &row : field.values)
        {
            min_i = std::min(min_i, row.i);
            max_i = std::max(max_i, row.i);
            min_j = std::min(min_j, row.j);
            max_j = std::max(max_j, row.j);
            if (row.j == -1 && row.i >= 0 && row.i <= 23)
                pfr1 -= row.value;
            if (row.j == -1 && row.i >= 72 && row.i <= 95)
                pfr2 -= row.value;
            if (row.j == 36 && row.i >= 0 && row.i <= 95)
                outer += row.value;
        }
        if (field.values.size() != 98u * 38u ||
            min_i != -1 || max_i != 96 || min_j != -1 || max_j != 36)
            throw std::runtime_error("unexpected FNIY index extent");
        if (!(pfr1 > 0.) || !(pfr2 > 0.) || !(outer > 0.))
            throw std::runtime_error("expected outward FNIY is not positive");

        std::cout << std::scientific << std::setprecision(12)
                  << "PFRSide1 " << pfr1 << '\n'
                  << "PFRSide2 " << pfr2 << '\n'
                  << "OuterSide " << outer << '\n'
                  << "Total " << pfr1 + pfr2 + outer << '\n';
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
