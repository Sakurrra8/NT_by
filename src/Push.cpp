#include "Particle.h"

Particle *Q;
Particle ParTemp;
void Push(Particle *PP)
{
    double step_now = 0;
    // while (PP->Weight() != 0)
    while (PP->Weight() > 1e-3)
    {
        if (PP->Charge() != 0)
        {
            exit(5);
        }
        step_now++;
        if (K_log)
        {
            if (PP != &H && PP != &H2 && PP != &D && PP != &D2 && PP != &T && PP != &T2)
            {
                Xlog << PP->name() << '_' << PP->Charge() << ' ' << PP->X(0) << ", " << PP->X(1) << ' ' << PP->X(2) << endl;
                Vlog << PP->name() << '_' << PP->Charge() << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << ' ' << PP->Vreal(7) << endl;
            }
        }
        P = PP;
        P->track();
        PP = P;
        if (PP->Weight() == 0)
            return;
        IfOut = 0;
        if (PP->Charge() != 0)
        {
            if (PP->Zone() >= 6)
            {
                /// search if the charged particle fly out of grid
                int num_intersect = 0;
                for (int i = 0; i < 98; i++)
                {
                    if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[i][1][3], Grid[i][1][7], Grid[i][1][2], Grid[i][1][6], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                    {
                        InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                        InterscePoint[num_intersect++][3] = i;
                    }
                }
                if (num_intersect > 1)
                {
                    for (int i = 1; i < num_intersect - 1; i++)
                    {
                        if (InterscePoint[i][0] < InterscePoint[0][0])
                            for (int j = 0; j < 4; j++)
                                InterscePoint[0][j] = InterscePoint[i][j];
                    }
                }
                if (num_intersect > 0)
                {
                    if (!backGridBoundry)
                    {
                        PP->SetWeight(0.);
                        return;
                    }
                    else if (random() < 0.95)
                    {
                        return;
                        // PP->track();
                    }
                }

                /// search if the charged particle fly into core
                num_intersect = 0;
                for (int i = 25; i < 73; i++)
                {
                    if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[i][36][0], Grid[i][36][4], Grid[i][36][1], Grid[i][36][5], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                    {
                        InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                        InterscePoint[num_intersect++][3] = i;
                    }
                }
                if (num_intersect > 1)
                {
                    for (int i = 1; i < num_intersect - 1; i++)
                    {
                        if (InterscePoint[i][0] < InterscePoint[0][0])
                            for (int j = 0; j < 4; j++)
                                InterscePoint[0][j] = InterscePoint[i][j];
                    }
                }
                if (num_intersect > 0)
                {
                    PP->SetWeight(0.);
                    return;
                }

                /// search if the charged particle fly to target
                num_intersect = 0;
                for (int i = 24; i < 37; i++)
                    if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid1.Wall(i, 0), Grid1.Wall(i, 1), Grid1.Wall(i + 1, 0), Grid1.Wall(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                    {
                        InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                        InterscePoint[num_intersect++][3] = i;
                    }
                if (num_intersect == 0)
                {
                    for (int i = 0; i < 38; i++)
                    {
                        if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[1][i][0], Grid[1][i][4], Grid[1][i][3], Grid[1][i][7], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                        {
                            InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                            InterscePoint[num_intersect++][3] = i;
                        }
                        if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[96][i][1], Grid[96][i][5], Grid[96][i][2], Grid[96][i][6], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                        {
                            InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                            InterscePoint[num_intersect++][3] = i + 38;
                        }
                    }
                }
                else
                {
                    InterscePoint[num_intersect][0] = PP->XY(1);
                }
                if (num_intersect > 1)
                {
                    for (int i = 1; i < num_intersect - 1; i++)
                    {
                        if (InterscePoint[i][0] < InterscePoint[0][0])
                            for (int j = 0; j < 4; j++)
                                InterscePoint[0][j] = InterscePoint[i][j];
                    }
                }

                if (num_intersect > 0)
                {
                    if (PP == &H2 || PP == &D2 || PP == &T2)
                        PP->Init(3, 2);
                    PP->setfate(0, 2, 20);
                    P = PP;
                    P->track();
                    PP = P;
                    if (PP->Weight() == 0)
                        return;
                }
                else
                {
                    PP->SetWeight(0.);
                    if (StepLog)
                    {
                        std::cout << "Particle " + PP->name() + " Charged run out of Grid" << endl;
                        std::cout << PP->NowZone() << endl;
                        std::cout << PP->name() << ' ' << PP->X(0) << ' ' << PP->X(1) << ' ' << PP->X(2) << endl;
                        std::cout << PP->name() << ' ' << PP->X_new(0) << ' ' << PP->X_new(1) << ' ' << PP->X_new(2) << endl;
                        std::cout << PP->name() << ' ' << PP->V(0) << ' ' << PP->V(1) << ' ' << PP->V(2) << endl;
                        std::cout << "V_P" << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << endl;
                        std::cout << "V_R" << ' ' << PP->Vreal(4) << ' ' << PP->Vreal(5) << ' ' << PP->Vreal(6) << ' ' << PP->Vreal(7) << endl;
                        std::cout << B[PP->XY(0)][PP->XY(1)][0] << " " << B[PP->XY(0)][PP->XY(1)][1] << " " << B[PP->XY(0)][PP->XY(1)][2] << endl;
                        std::cout << Num_Reflect_Core << '\t' << PP->Charge() << endl;
                    }
                    return;
                }
            }
        }
        if (PP->Charge() == 0)
        {
            if (K_flight == 1 || K_flight == 3)
            {
                while (PP->Zone() >= 6)
                {
                    if (K_test1)
                    {
                        std::cout << PP->name() << " X1: " << PP->X(0) << ',' << PP->X(1) << '\t' << PP->X_new(0) << ',' << PP->X_new(1) << '\t';
                        std::cout << PP->Zone() << endl;
                    }
                    int num_intersect = 0;
                    if (PP->Zone() == 6)
                    {
                        P = PP;
                        P->track();
                        PP = P;
                        if (PP->Weight() == 0)
                            return;
                    }
                    /// if the test particle flight out from
                    if (PP->IfFlightOut() == 7)
                    {
                        num_intersect = 1;
                        PP->SetIfFlightOut(0);
                        if (InterscePoint[0][3] < N_radial)
                        {
                            PP->SetXY(0, 1);
                            PP->SetXY(1, InterscePoint[0][3]);
                        }
                        else
                        {
                            PP->SetXY(0, 96);
                            PP->SetXY(1, InterscePoint[0][3] - N_radial);
                        }

                        if (PP->XY(1) >= 19 && PP->XY(1) <= 36)
                            PP->SetZone(3);
                        else if (PP->XY(0) <= 24)
                            PP->SetZone(4);
                        else if (PP->XY(0) <= 72)
                            PP->SetZone(2);
                        else
                            PP->SetZone(5);
                    }
                    else
                    {
                        if (K_test1)
                        {
                            std::cout << "Find2: " << '\t';
                        }
                        for (int i = 0; i < Grid1.Wall_num(); i++)
                        {
                            if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid1.Wall(i, 0), Grid1.Wall(i, 1), Grid1.Wall(i + 1, 0), Grid1.Wall(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                            {
                                if (Tools::PointandLine(Grid1.Wall(i, 0), Grid1.Wall(i, 1), Grid1.Wall(i + 1, 0), Grid1.Wall(i + 1, 1), PP->X_new(0), PP->X_new(1)))
                                {
                                    InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                                    InterscePoint[num_intersect][3] = i;
                                    InterscePoint[num_intersect++][4] = 11;
                                }
                                else
                                {
                                    std::cout << "no, not here!" << endl;
                                }
                            }
                        }
                        if (K_test1)
                        {
                            std::cout << num_intersect << endl;
                        }
                        if (K_test1)
                        {
                            std::cout << "Find1: " << '\t';
                        }
                        // std::cout << PP->X(0) << '\t' << PP->X(1) << '\t' << PP->X_new(0) << '\t' << PP->X_new(1) << endl;
                        for (int i = 0; i < Grid1.PLasma_Grid_Boundry_num(); i++)
                        {
                            if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 0), Grid1.PLasma_Grid_Boundry(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                            {
                                if (!Tools::PointandLine(Grid1.PLasma_Grid_Boundry(i, 0), Grid1.PLasma_Grid_Boundry(i, 1), Grid1.PLasma_Grid_Boundry(i + 1, 0), Grid1.PLasma_Grid_Boundry(i + 1, 1), PP->X_new(0), PP->X_new(1)))
                                // if (InterscePoint[num_intersect][1] != PP->X(0) && InterscePoint[num_intersect][2] != PP->X(1))
                                {
                                    InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                                    InterscePoint[num_intersect][3] = i;
                                    InterscePoint[num_intersect++][4] = 19;
                                }
                            }
                        }
                    }
                    if (K_test1)
                    {
                        std::cout << num_intersect << endl;
                    }
                    if (num_intersect > 1)
                    {
                        for (int i = 1; i < num_intersect; i++)
                        {
                            if (InterscePoint[i][0] < InterscePoint[0][0])
                                // if (InterscePoint[i][1] != PP->X(0) && InterscePoint[i][2] != PP->X(1))
                                for (int j = 0; j < 5; j++)
                                    InterscePoint[0][j] = InterscePoint[i][j];
                        }
                    }
                    if (num_intersect < 1)
                    {
                        std::cout << PP->name() << ',' << PP->Charge() << ',' << fatename(PP->fate(0)) << ',' << sourcename(PP->sourcePar(0)) << endl;
                        std::cout << "X," << PP->X(0) << ',' << PP->X(1) << ';' << PP->X_new(0) << ',' << PP->X_new(1) << endl;
                        std::cout << "V," << PP->V(0) << ',' << PP->V(1) << endl;
                        std::cout << "I don't know where the Particle is." << endl;
                        std::cout << endl;
                        // PP->track();
                        // exit(0);
                        pause();
                    }
                    else
                    {
                        if (K_test1)
                        {
                            std::cout << "InterscePoint: " << InterscePoint[0][1] << ',' << InterscePoint[0][2] << '\t' << sourcename((int)InterscePoint[0][4]) << endl;
                        }
                        if (InterscePoint[0][4] == 11 || InterscePoint[0][4] == 1) // wall
                        {
                            PP->SetX(0, InterscePoint[0][1]);
                            PP->SetX(1, InterscePoint[0][2]);
                            P = PP;
                            WallReflect();
                            if (P->Weight() == 0)
                                return;
                            P->track();
                            PP = P;
                            // std::cout << "reflect: " << PP->Zone() << ", ";
                            // std::cout << PP->X(0) << ", " << PP->X(1) << " " << PP->X_new(0) << ", " << PP->X_new(1) << endl;
                            if (PP->Weight() == 0)
                                return;
                            // continue;

                            /*if (!Grid1.IfinIndex1(PP->XY(0), PP->XY(1), PP->X_new(0), PP->X_new(1)))
                            {
                                std::cout << PP->name() + ": " << PP->X(0) << "," << PP->X(1) << " " << PP->X_new(0) << "," << PP->X_new(1) << " " << PP->XY(0) << "," << PP->XY(1) << endl;
                                std::cout << "Push.cpp:328 have some problem.\n";
                            }*/
                        }
                        else if (InterscePoint[0][4] == 19)
                        {
                            // std::cout << PP->X(0) << " " << PP->X(1) << " " << InterscePoint[0][1] << " " << InterscePoint[0][2] << endl;
                            // std::cout << "intersecPoint: " << InterscePoint[0][1] << ", " << InterscePoint[0][2] << endl;
                            PP->SetX(0, InterscePoint[0][1]);
                            PP->SetX(1, InterscePoint[0][2]);
                            // std::cout << InterscePoint[0][3] << "\t";
                            if (InterscePoint[0][3] < N_radial - 2)
                            {
                                PP->SetXY(0, 1);
                                PP->SetXY(1, InterscePoint[0][3] + 1);
                                PP->SetZone(ZonefromXY(PP->XY(0), PP->XY(1)));
                                // std::cout << PP->XY(0) << ", " << PP->XY(1) << endl;
                            }
                            else if (InterscePoint[0][3] < N_poloidal + N_radial - 4)
                            {
                                PP->SetXY(0, InterscePoint[0][3] - N_radial + 3);
                                PP->SetXY(1, 36);
                                PP->SetZone(ZonefromXY(PP->XY(0), PP->XY(1)));
                                // std::cout << PP->XY(0) << ", " << PP->XY(1) << endl;
                            }
                            else if (InterscePoint[0][3] < N_poloidal + 2 * N_radial - 6)
                            {
                                PP->SetXY(0, 96);
                                PP->SetXY(1, N_radial - 2 - (InterscePoint[0][3] - N_poloidal - N_radial + 4));
                                PP->SetZone(ZonefromXY(PP->XY(0), PP->XY(1)));
                                // std::cout << PP->XY(0) << ", " << PP->XY(1) << endl;
                            }
                            else if (InterscePoint[0][3] < N_poloidal + 2 * N_radial + 18)
                            {
                                PP->SetXY(0, N_poloidal - 2 - (InterscePoint[0][3] - (N_poloidal + 2 * N_radial - 6)));
                                PP->SetXY(1, 1);
                                PP->SetZone(ZonefromXY(PP->XY(0), PP->XY(1)));
                                // std::cout << PP->XY(0) << ", " << PP->XY(1) << endl;
                                //  std::cout << PP->name() + ": " << InterscePoint[0][3] << " " << PP->XY(0) << "," << PP->XY(1) << endl;
                            }
                            else if (InterscePoint[0][3] < N_poloidal + 2 * N_radial + 42)
                            {
                                PP->SetXY(0, 25 - 1 - (InterscePoint[0][3] - (N_poloidal + 2 * N_radial + 18)));
                                PP->SetXY(1, 1);
                                PP->SetZone(ZonefromXY(PP->XY(0), PP->XY(1)));
                                // std::cout << PP->XY(0) << ", " << PP->XY(1) << endl;
                            }
                            else
                            {
                                std::cout << "error!" << endl;
                                exit(0);
                            }
                            if (K_test1)
                            {
                                std::cout << "goodbye!" << " Zone = " << PP->Zone() << endl;
                                std::cout << "XY: " << PP->XY(0) << ',' << PP->XY(1) << endl;
                                // pause();}
                            }
                            /*if (!Grid1.IfinIndex1(PP->XY(0), PP->XY(1), PP->X_new(0), PP->X_new(1)))
                            {
                                std::cout << PP->name() + ": " << PP->X(0) << "," << PP->X(1) << " " << PP->X_new(0) << "," << PP->X_new(1) << " " << PP->XY(0) << "," << PP->XY(1) << endl;
                                std::cout << "the Push.cpp:328 have some problem.\n";
                            }*/
                            P = PP;
                            P->track();
                            PP = P;
                            if (PP->Weight() == 0)
                                return;
                        }
                    }
                    /*else
                    {
                        std::cout << "what is the error???" << endl;
                        exit(0);
                    }*/
                }
            }
            else if (K_flight == 2)
            {
                while (PP->Zone() == 7)
                {
                    /*IfOut += 1;
                    if (IfOut > 100)
                    {
                        std::cerr << "Your Particle may run le!" << endl;
                        if (PP->NowZone() > 6)
                            exit(1);
                    }*/
                    int num_intersect = 0;
                    // std::cout << PP->X(0) << '\t' << PP->X(1) << '\t' << PP->X_new(0) << '\t' << PP->X_new(1) << endl;
                    for (int i = 0; i < Grid1.Wall_num(); i++)
                    {
                        if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid1.Wall(i, 0), Grid1.Wall(i, 1), Grid1.Wall(i + 1, 0), Grid1.Wall(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                        {
                            if (InterscePoint[num_intersect][1] != PP->X(0) && InterscePoint[num_intersect][2] != PP->X(1))
                            {
                                InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                                InterscePoint[num_intersect++][3] = i;
                            }
                        }
                    }
                    if (num_intersect > 0)
                    {
                        Num_Reflect++;
                        // std::cout << InterscePoint[num_intersect - 1][3] << endl;
                    }
                    else
                    {
                        std::cerr << PP->name() << ' ' << PP->NowZone() << '\t' << PP->Zone() << '\t' << PP->Charge() << '+' << '\t' << PP->Tn() << '\t' << fatename(PP->fate(0)) << '\t' << sourcename(PP->sourcePar(0)) << '\t' << PP->Weight() << endl;
                        std::cerr << PP->sourceWall(0) << '\t' << PP->sourceWall(1) << '\t' << PP->sourceWall(2) << '\t' << PP->sourceWall(3) << endl;
                        std::cerr << "X: " << PP->X(0) << ',' << PP->X(1) << ' ' << PP->X(2) << endl;
                        std::cerr << "V_P" << ' ' << PP->Vreal(0) << ',' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << endl;
                        std::cerr << "V_R" << ' ' << PP->Vreal(4) << ',' << PP->Vreal(5) << ' ' << PP->Vreal(6) << ' ' << PP->Vreal(7) << endl;
                        std::cerr << "X_new: " << ',' << PP->X_new(0) << ' ' << PP->X_new(1) << ' ' << PP->X_new(2) << endl;
                        std::cerr << Num_Reflect << '\t' << PP->Charge() << '\t' << PP->Tn() << endl;
                        for (int i = 0; i < 3; i++)
                            PP->SetX(i, PP->X_new(i));
                        std::cerr << "Particle " + PP->name() + " have run;\n";
                        if (K_log)
                        {
                            Xlog << PP->Zone() << '\t' << PP->name() << '_' << PP->Charge() << ' ' << PP->X_new(0) << ", " << PP->X_new(1) << ' ' << PP->X_new(2) << endl;
                            Vlog << PP->name() << '_' << PP->Charge() << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << ' ' << PP->Vreal(7) << endl;
                        }
                        exit(1);
                    }
                    if (num_intersect > 1)
                    {
                        for (int i = 1; i < num_intersect - 1; i++)
                        {
                            if (InterscePoint[i][0] < InterscePoint[0][0])
                                for (int j = 0; j < 4; j++)
                                    InterscePoint[0][j] = InterscePoint[i][j];
                        }
                    }
                    PP->FluxCal_Target();
                    // Tn_ = (pow(V_[0], 2) + pow(V_[1], 2) + pow(V_[2], 2)) * mass_ / 3. / qe;
                    // set the recycling coefficient
                    /*if (InterscePoint[0][3] == 27 || InterscePoint[0][3] == 34)
                    {
                        Weight_ = 0;
                        return;
                    }*/
                    P = PP;
                    WallReflect();
                    PP = P;
                    if (PP->Weight() == 0)
                        return;
                }
            }
            if (PP->Zone() == 7)
            {
                std::cerr << PP->name() << ',' << PP->NowZone() << '\t' << PP->Zone() << '\t' << PP->Charge() << '+' << '\t' << PP->Tn() << '\t' << fatename(PP->fate(0)) << '\t' << sourcename(PP->sourcePar(0)) << '\t' << PP->Weight() << endl;
                std::cerr << "InterscePoint" << ',' << InterscePoint[0][1] << '\t' << InterscePoint[0][2] << endl;
                std::cerr << "X" << ' ' << PP->X(0) << ',' << PP->X(1) << ',' << PP->X(2) << endl;
                std::cerr << "V" << ' ' << PP->Vreal(0) << ',' << PP->Vreal(1) << ' ' << PP->Vreal(2) << endl;
                std::cerr << "X_new" << ' ' << PP->X_new(0) << ',' << PP->X_new(1) << ' ' << PP->X_new(2) << endl;
                std::cerr << K_backScatter << '\t' << PP->Charge() << endl;
                std::cerr
                    << "The Paiticle " << PP->name() << " " << "is out of Tokamak." << endl;
                exit(2);
            }
            if (K_log)
            {
                if (PP != &H && PP != &H2 && PP != &D && PP != &D2 && PP != &T && PP != &T2)
                {
                    Xlog << PP->name() << '_' << PP->Charge() << ' ' << PP->X(0) << ", " << PP->X(1) << ' ' << PP->X(2) << endl;
                    Vlog << PP->name() << '_' << PP->Charge() << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << ' ' << PP->Vreal(7) << endl;
                }
            }
        }
        if (K_flight == 2)
        {
            if (PP->Weight() != 0)
                PP->NumParStat();
        }

        P = PP;
        if (K_flight == 1 || K_flight == 3)
        {
            if (P->Zone() < 6 && P->Zone() > 1 && P->IfColl())
            {
                P->Coll();
            }
        }
        else if (K_flight == 2)
        {
            if (P->Zone() < 6 && P->Zone() > 1)
            {
                P->Coll();
            }
        }
        PP = P;
        /*if (Change.Num_ParChan())
        {
            // std::cout << "Coll endl;" << endl;
            Change.ParChange();
        }*/
    }
}

void Push_Tri(Particle *PP)
{
    while (PP->Weight() > 1e-3)
    {
        if (PP->Charge() != 0)
        {
            exit(5);
        }
        if (K_log)
        {
            if (PP != &H && PP != &H2 && PP != &D && PP != &D2 && PP != &T && PP != &T2)
            {
                Xlog << PP->name() << '_' << PP->Charge() << ' ' << PP->X(0) << ", " << PP->X(1) << ' ' << PP->X(2) << endl;
                Vlog << PP->name() << '_' << PP->Charge() << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << ' ' << PP->Vreal(7) << endl;
            }
        }

        P = PP;
        P->track();
        PP = P;

        if (PP->Weight() < 1e-3)
            return;
        else if (PP->IfColl())
        {
            P = PP;
            P->Coll();
            PP = P;
        }
        else if (PP->Zone() > 6)
        {
            int num_intersect = 0; /*
             /// about target reflect
             for (int i = 0; i < 38; i++)
             {
                 if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[1][i][0], Grid[1][i][4], Grid[1][i][3], Grid[1][i][7], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                 {
                     InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                     InterscePoint[num_intersect++][3] = i;
                 }
                 if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid[96][i][1], Grid[96][i][5], Grid[96][i][2], Grid[96][i][6], &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                 {
                     InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                     InterscePoint[num_intersect++][3] = i + 38;
                 }
             }

            /// about wall reflect
            if (InterscePoint[0][4] == 11)
            {
                if (num_intersect == 0)
                {
                    for (int i = 0; i < Grid4.Wall_.Wall_num(); i++)
                    {
                        if (Tools::get_line_intersection(PP->X(0), PP->X(1), PP->X_new(0), PP->X_new(1), Grid1.Wall(i, 0), Grid1.Wall(i, 1), Grid1.Wall(i + 1, 0), Grid1.Wall(i + 1, 1), &InterscePoint[num_intersect][1], &InterscePoint[num_intersect][2]))
                        {
                            InterscePoint[num_intersect][0] = sqrt(pow((InterscePoint[num_intersect][1] - PP->X(0)), 2) + pow((InterscePoint[num_intersect][2] - PP->X(1)), 2));
                            InterscePoint[num_intersect++][3] = i;
                        }
                    }
                }
                if (num_intersect > 1)
                {
                    for (int i = 1; i < num_intersect - 1; i++)
                    {
                        if (InterscePoint[i][0] < InterscePoint[0][0])
                            for (int j = 0; j < 4; j++)
                                InterscePoint[0][j] = InterscePoint[i][j];
                    }
                }
                if (num_intersect > 0)
                {
                    if (PP == &H2 || PP == &D2 || PP == &T2)
                        PP->Init(3, 2);
                    PP->setfate(0, 2, 20);
                    P = PP;
                    P->track();
                    PP = P;
                    if (PP->Weight() == 0)
                        return;
                }
                else
                {
                    if (StepLog)
                    {
                        std::cout << "Particle " + PP->name() + " Charged run out of Grid" << endl;
                        std::cout << PP->NowZone() << endl;
                        std::cout << PP->name() << ' ' << PP->X(0) << ' ' << PP->X(1) << ' ' << PP->X(2) << endl;
                        std::cout << PP->name() << ' ' << PP->X_new(0) << ' ' << PP->X_new(1) << ' ' << PP->X_new(2) << endl;
                        std::cout << PP->name() << ' ' << PP->V(0) << ' ' << PP->V(1) << ' ' << PP->V(2) << endl;
                        std::cout << "V_P" << ' ' << PP->Vreal(0) << ' ' << PP->Vreal(1) << ' ' << PP->Vreal(2) << ' ' << PP->Vreal(3) << endl;
                        std::cout << "V_R" << ' ' << PP->Vreal(4) << ' ' << PP->Vreal(5) << ' ' << PP->Vreal(6) << ' ' << PP->Vreal(7) << endl;
                        std::cout << B[PP->XY(0)][PP->XY(1)][0] << " " << B[PP->XY(0)][PP->XY(1)][1] << " " << B[PP->XY(0)][PP->XY(1)][2] << endl;
                        std::cout << Num_Reflect_Core << '\t' << PP->Charge() << endl;
                    }
                    pause();
                    PP->SetWeight(0.);
                    return;
                }
            }*/

            P = PP;
            WallReflect();
            if (P->Weight() == 0)
                return;
            P->track();
            PP = P;
        }
    }
}