#include "Global.h"
#include "Particle.h"

void Moncar()
{
    if (K_H)
    {
        H.Clear();
        H2.Clear();
    }
    if (K_D)
    {
        D.Clear();
        D2.Clear();
    }
    if (K_T)
    {
        T.Clear();
        T2.Clear();
    }
    if (K_Methane)
    {
        CD4.Clear();
        CD3.Clear();
        CD2.Clear();
        CD1.Clear();
    }
    if (K_Methane || K_C)
    {
        C.Clear();
    }
    if (K_Ar)
    {
        Ar.Clear();
    }

    if (K_Recyc)
        for (int i = 0; i < N_radial * 2; i++)
        {
            if (K_H)
            {
                for (int j = 1; j <= H.Num_Target(i); j++)
                {
                    /// neutral particle transport from target

                    P = &H;
                    P->Init(1, i);
                    flight(i, j);
                    P = &H2;
                    P->Init(1, i);
                    flight(i, j);
                }
            }
            if (K_D)
            {
                for (int j = 1; j <= D.Num_Target(i); j++)
                {
                    P = &D;
                    P->Init(1, i);
                    flight(i, j);
                    P = &D2;
                    P->Init(1, i);
                    flight(i, j);
                }
            }
            if (K_T)
            {
                for (int j = 1; j <= T.Num_Target(i); j++)
                {
                    P = &T;
                    P->Init(1, i);
                    flight(i, j);
                    P = &T2;
                    P->Init(1, i);
                    flight(i, j);
                }
            }
            if (StepLog)
                std::cout << "Target flight " << i << " finished;" << endl;
        }

    /// neutral particle transport from Recom and Grid
    if (K_Rec)
    {
        if (MeshMode == 1)
        {
            for (int i = 1; i < N_poloidal - 1; i++)
            {
                for (int j = 1; j < N_radial - 1; j++)
                {
                    for (int k = 1; k <= numPar_flight; k++)
                        if (K_H)
                        {
                            P = &H;
                            P->Init(4, i * 38 + j);
                            flight(i, j);
                        }
                    for (int k = 1; k <= numPar_flight; k++)
                        if (K_D)
                        {
                            P = &D;
                            P->Init(4, i * 38 + j);
                            flight(i, j);
                        }
                    for (int k = 1; k <= numPar_flight; k++)
                        if (K_T)
                        {
                            P = &T;
                            P->Init(4, i * 38 + j);
                            flight(i, j);
                        }
                    if (StepLog)
                        cout << "D Grid " << i << ' ' << j << " finished;" << endl;
                }
            }
        }
        else if (MeshMode == 3)
        {
            for (int i = 0; i < Grid4.num_tris(); i++)
            {
                if (Grid4.if_in_plasmagrid(i))
                {
                    for (int k = 1; k <= numPar_flight; k++)
                    {
                        if (K_H)
                        {
                            P = &H;
                            P->Init(4, i);
                            flight(i, 0);
                        }
                        if (K_D)
                        {
                            P = &D;
                            P->Init(4, i);
                            flight(i, 0);
                        }
                        if (K_T)
                        {
                            P = &T;
                            P->Init(4, i);
                            flight(i, 0);
                        }
                    }
                    if (StepLog)
                        cout << "D Grid " << i << " finished;" << endl;
                }
            }
        }
    }
    if (K_Pump)
    {
        if (K_D)
            for (int i = 0; i < numPar_flight_D2; i++)
            {
                P = &D2;
                P->Init(8, 0);
                flight(i, 0);
                if (StepLog)
                    cout << "Puff " << i << " finished;" << endl;
            }
        if (K_T)
            for (int i = 0; i < numPar_flight_D2; i++)
            {
                P = &T2;
                P->Init(8, 0);
                flight(i, 0);
                if (StepLog)
                    cout << "Puff " << i << " finished;" << endl;
            }
    }
    if (K_Methane)
    {
        for (int i = 0; i < numPar_flight_CD4; i++)
        {
            /*if (i % 1000 == 0 && i != numPar_flight_CD4)
            {
                D.Stat(2);
                D2.Stat(2);
                CD4.Stat(2);
                CD3.Stat(2);
                CD2.Stat(2);
                CD1.Stat(2);
                C.Stat(2);
                Ar.Stat(2);
                Output();
            }*/
            P = &CD4;
            P->Init(5, 0);
            flight(i, 0);
            if (StepLog)
                cout << "CD4 Puff " << i << " finished;" << endl;
        }
    }
    if (K_Ar)
    {
        for (int i = 0; i < numPar_flight_CD4; i++)
        {
            /*if (i % 1000 == 0 && i != numPar_flight_CD4)
            {
                D.Stat(2);
                D2.Stat(2);
                CD4.Stat(2);
                CD3.Stat(2);
                CD2.Stat(2);
                CD1.Stat(2);
                C.Stat(2);
                Ar.Stat(2);
                Output();
            }*/
            P = &Ar;
            P->Init(5, 0);
            flight(i, 0);
            if (StepLog)
                cout << "Ar Puff " << i << " finished;" << endl;
        }
    }
    if (K_C)
    {
        for (int i = 0; i < numPar_flight_CD4; i++)
        {
            /*if (i % 1000 == 0 && i != numPar_flight_CD4)
            {
                D.Stat(2);
                D2.Stat(2);
                CD4.Stat(2);
                CD3.Stat(2);
                CD2.Stat(2);
                CD1.Stat(2);
                C.Stat(2);
                Ar.Stat(2);
                Output();
            }*/
            P = &C;
            P->Init(5, 0);
            flight(i, 0);
            if (StepLog)
                cout << "C Puff " << i << " finished;" << endl;
            P = &D;
            P->Init(5, 0);
            flight(i, 0);
            if (StepLog)
                cout << "D2 Puff " << i << " finished;" << endl;
        }
    }
    if (K_H)
    {
        H.Stat(2);
        H2.Stat(2);
        H.Stat_Tri(2);
        H2.Stat_Tri(2);
    }
    if (K_D)
    {
        D.Stat(2);
        D2.Stat(2);
        D.Stat_Tri(2);
        D2.Stat_Tri(2);
    }
    if (K_T)
    {
        T.Stat(2);
        T2.Stat(2);
        T.Stat_Tri(2);
        T2.Stat_Tri(2);
    }
    if (K_Methane)
    {
        CD4.Stat(2);
        CD3.Stat(2);
        CD2.Stat(2);
        CD1.Stat(2);
        C.Stat(2);
    }
    if (K_C)
    {
        Ar.Stat(2);
    }
}

void flight(int i, int j)
{
    // out.open("test.txt");
    if (NumPar_now <= 1.)
        P->SetWeight(0.);
    Num_Reflect = 0;
    Num_Reflect_Core = 0;
    if (K_log)
    {
        name_Xlog = "Outputfile/test/log/" + P->name() + "_1_Xlog";
        name_Xlog += "_" + to_string(i) + "_" + to_string(j) + ".txt";
        name_Vlog = "Outputfile/test/log/" + P->name() + "_1_Vlog";
        name_Vlog += "_" + to_string(i) + "_" + to_string(j) + ".txt";
        Xlog.open(name_Xlog);
        Vlog.open(name_Vlog);
    }
    // step_now = 0;
    if (MeshMode == 1)
        Push(P);
    else if (MeshMode == 3)
        Push_Tri(P);
    // PartoPar_Temp.ParChange();
    //  out.close();
    if (K_log)
    {
        Xlog.close();
        Vlog.close();
        if (P != &C)
        {
            // remove(name_Xlog.c_str());
            // remove(name_Vlog.c_str());
        }
    }
}