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
    {
        if (K_H)
        {
            if (H.RecycledSourceSum() > 0.)
            {
                H.BeginDeferredFlightStats(H.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = H.sampleRecycledTarget();
                    P = &H;
                    P->Init(1, target);
                    flight(target, j);
                }
                H.EndDeferredFlightStats();
            }

            if (H2.RecycledSourceSum() > 0.)
            {
                H2.BeginDeferredFlightStats(H2.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = H2.sampleRecycledTarget();
                    P = &H2;
                    P->Init(1, target);
                    flight(target, j);
                }
                H2.EndDeferredFlightStats();
            }
        }
        if (K_D)
        {
            if (D.RecycledSourceSum() > 0.)
            {
                D.BeginDeferredFlightStats(D.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = D.sampleRecycledTarget();
                    P = &D;
                    P->Init(1, target);
                    flight(target, j);
                }
                D.EndDeferredFlightStats();
            }

            if (D2.RecycledSourceSum() > 0.)
            {
                D2.BeginDeferredFlightStats(D2.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = D2.sampleRecycledTarget();
                    P = &D2;
                    P->Init(1, target);
                    flight(target, j);
                }
                D2.EndDeferredFlightStats();
            }
        }
        if (K_T)
        {
            if (T.RecycledSourceSum() > 0.)
            {
                T.BeginDeferredFlightStats(T.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = T.sampleRecycledTarget();
                    P = &T;
                    P->Init(1, target);
                    flight(target, j);
                }
                T.EndDeferredFlightStats();
            }

            if (T2.RecycledSourceSum() > 0.)
            {
                T2.BeginDeferredFlightStats(T2.NumPar_Target());
                for (int j = 1; j <= numPar_flight_Target; j++)
                {
                    const int target = T2.sampleRecycledTarget();
                    P = &T2;
                    P->Init(1, target);
                    flight(target, j);
                }
                T2.EndDeferredFlightStats();
            }
        }
    }

    /// neutral particle transport from Recom and Grid
    if (K_Rec)
    {
        if (K_H && H.RecombinSourceSum() > 0.)
        {
            H.BeginDeferredFlightStats(H.NumPar_Grid());
            for (int k = 1; k <= numPar_flight; k++)
            {
                const int source = H.sampleRecombinCell();
                P = &H;
                P->Init(4, source);
                flight(source, k);
            }
            H.EndDeferredFlightStats();
        }
        if (K_D && D.RecombinSourceSum() > 0.)
        {
            D.BeginDeferredFlightStats(D.NumPar_Grid());
            for (int k = 1; k <= numPar_flight; k++)
            {
                const int source = D.sampleRecombinCell();
                P = &D;
                P->Init(4, source);
                flight(source, k);
            }
            D.EndDeferredFlightStats();
        }
        if (K_T && T.RecombinSourceSum() > 0.)
        {
            T.BeginDeferredFlightStats(T.NumPar_Grid());
            for (int k = 1; k <= numPar_flight; k++)
            {
                const int source = T.sampleRecombinCell();
                P = &T;
                P->Init(4, source);
                flight(source, k);
            }
            T.EndDeferredFlightStats();
        }
    }
    if (K_Pump)
    {
        if (K_D && Num_D2_pump > 0. && numPar_flight_D2 > 0)
        {
            D2.BeginDeferredFlightStats(Num_D2_pump / (double)numPar_flight_D2);
            for (int i = 0; i < numPar_flight_D2; i++)
            {
                P = &D2;
                P->Init(8, 0);
                flight(i, 0);
                if (StepLog)
                    cout << "Puff " << i << " finished;" << endl;
            }
            D2.EndDeferredFlightStats();
        }
        if (K_T && Num_T2_pump > 0. && numPar_flight_T2 > 0)
        {
            T2.BeginDeferredFlightStats(Num_T2_pump / (double)numPar_flight_T2);
            for (int i = 0; i < numPar_flight_D2; i++)
            {
                P = &T2;
                P->Init(8, 0);
                flight(i, 0);
                if (StepLog)
                    cout << "Puff " << i << " finished;" << endl;
            }
            T2.EndDeferredFlightStats();
        }
    }
    if (K_Methane && Num_CD4_pump > 0. && numPar_flight_CD4 > 0)
    {
        CD4.BeginDeferredFlightStats(Num_CD4_pump / (double)numPar_flight_CD4);
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
        CD4.EndDeferredFlightStats();
    }
    if (K_Ar && Num_CD4_pump > 0. && numPar_flight_CD4 > 0)
    {
        Ar.BeginDeferredFlightStats(Num_CD4_pump / (double)numPar_flight_CD4);
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
        Ar.EndDeferredFlightStats();
    }
    if (K_C && Num_CD4_pump > 0. && numPar_flight_CD4 > 0)
    {
        C.BeginDeferredFlightStats(Num_CD4_pump / (double)numPar_flight_CD4);
        D.BeginDeferredFlightStats(4. * Num_CD4_pump / (double)numPar_flight_CD4);
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
        C.EndDeferredFlightStats();
        D.EndDeferredFlightStats();
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