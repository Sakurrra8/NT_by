#include "Global.h"
#include "Particle.h"

void HydrogenOutput(Particle *OutPar1, Particle *OutPar2);
void HydrogenOutput_Tri(Particle *OutPar1, Particle *OutPar2);

string n = "n", E = "E", Tn = "T", Sn = "Sn", SE = "SE", SE_n = "SE_n", Smu = "Smu", Smu1 = "Smu1", Smu2 = "Smu2", Pra = "Pra";
string Ion = "Ion", Rec = "Rec", CX = "CX", Diss1 = "Diss1", Diss2 = "Diss2", Diss3 = "Diss3", Mar = "Mar", Ela = "Ela", CXDT = "CXDT", R_with_H = "R_with_H", R_with_H2 = "R_with_H2";
string DS[4] = {
    "DS1",
    "DS2",
    "DS3",
    "DS4",
};
void Output()
{

    ofstream Out_temp;
    if (K_H)
    {
        HydrogenOutput(&H, &H2);
        if (MeshMode == 3)
        {
            HydrogenOutput_Tri(&H, &H2);
        }

        Out_temp.open(Outputpath + "recycling_H.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_H_recyc[i] << " " << Tn_H_recyc[i] << " " << NumPar_H2_recyc[i] << endl;
        }
        Out_temp.close();
    }
    if (K_D)
    {
        HydrogenOutput(&D, &D2);
        if (MeshMode == 3)
            HydrogenOutput_Tri(&D, &D2);

        Out_temp.open(Outputpath + "recycling_D.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_D_recyc[i] << " " << Tn_D_recyc[i] << " " << NumPar_D2_recyc[i] << endl;
        }
        Out_temp.close();
    }
    if (K_T)
    {
        HydrogenOutput(&T, &T2);
        if (MeshMode == 3)
            HydrogenOutput_Tri(&T, &T2);

        Out_temp.open(Outputpath + "recycling_T.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_T_recyc[i] << " " << Tn_T_recyc[i] << " " << NumPar_T2_recyc[i] << endl;
        }
        Out_temp.close();
    }
    D_W_sputtered.writeH5(Outputpath + "D_W_stats.h5", "/D_W");

    Out_temp.open(Outputpath + "Bx");
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            Out_temp << B[i][j][0] << '\t';
        }
        Out_temp << endl;
    }
    Out_temp.close();
    Out_temp.open(Outputpath + "By");
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            Out_temp << B[i][j][1] << '\t';
        }
        Out_temp << endl;
    }
    Out_temp.close();
    Out_temp.open(Outputpath + "Bz");
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            Out_temp << B[i][j][2] << '\t';
        }
        Out_temp << endl;
    }
    Out_temp.close();
    if (K_Methane)
    {
        CD4.Dump_2D(n, n, 0);
        CD4.Dump_2D(n, n, 1);
        CD4.Dump_2D(Tn, Tn, 0);
        CD4.Dump_2D(Tn, Tn, 1);
        CD4.Dump_2D(Sn, Ion);
        CD4.Dump_2D(Sn, CX);
        CD4.Dump_2D(Sn, Diss1);
        CD4.Dump_2D(Sn, Diss2);
        for (int i = 0; i < 4; i++)
            CD4.Dump_2D(Sn, DS[i], 1);

        CD3.Dump_2D(n, n, 0);
        CD3.Dump_2D(n, n, 1);
        CD3.Dump_2D(Tn, Tn, 0);
        CD3.Dump_2D(Tn, Tn, 1);
        CD3.Dump_2D(Sn, Ion);
        CD3.Dump_2D(Sn, CX);
        CD3.Dump_2D(Sn, Diss1);
        CD3.Dump_2D(Sn, Diss2);
        for (int i = 0; i < 3; i++)
            CD3.Dump_2D(Sn, DS[i], 1);

        CD2.Dump_2D(n, n, 0);
        CD2.Dump_2D(n, n, 1);
        CD2.Dump_2D(Tn, Tn, 0);
        CD2.Dump_2D(Tn, Tn, 1);
        CD2.Dump_2D(Sn, Ion);
        CD2.Dump_2D(Sn, CX);
        CD2.Dump_2D(Sn, Diss1);
        CD2.Dump_2D(Sn, Diss2);
        for (int i = 0; i < 3; i++)
            CD2.Dump_2D(Sn, DS[i], 1);

        CD1.Dump_2D(n, n, 0);
        CD1.Dump_2D(n, n, 1);
        CD4.Dump_2D(Tn, Tn, 0);
        CD4.Dump_2D(Tn, Tn, 1);
        CD1.Dump_2D(Sn, Ion);
        CD1.Dump_2D(Sn, CX);
        CD1.Dump_2D(Sn, Diss1);
        CD1.Dump_2D(Sn, Diss2);
        CD1.Dump_2D(Sn, Diss3);
        for (int i = 0; i < 3; i++)
            CD1.Dump_2D(Sn, DS[i], 1);

        for (int i = 0; i < C.MaxCharge() + 1; i++)
        {
            C.Dump_2D(n, n, i);
            C.Dump_2D(Tn, Tn, i);
            if (i != C.MaxCharge())
            {
                C.Dump_2D(Sn, Ion, i); // ion source output
                C.Dump_2D(Sn, CX, i);  // CX source output
                C.Dump_2D(Pra, Ion, i);
                C.Dump_2D(Pra, CX, i);
            }
            if (i != 0)
            {
                C.Dump_2D(Sn, Rec, i); // rec source output
                C.Dump_2D(Pra, Rec, i);
            }
        }
    }
    if (K_C)
    {
        for (int i = 0; i < C.MaxCharge() + 1; i++)
        {
            C.Dump_2D(n, n, i);
            C.Dump_2D(Tn, Tn, i);
            if (i != C.MaxCharge())
            {
                C.Dump_2D(Sn, Ion, i); // ion source output
                C.Dump_2D(Sn, CX, i);  // CX source output
                C.Dump_2D(Pra, Ion, i);
                C.Dump_2D(Pra, CX, i);
            }
            if (i != 0)
            {
                C.Dump_2D(Sn, Rec, i); // rec source output
                C.Dump_2D(Pra, Rec, i);
            }
        }
    }
    if (K_Ar)
    {
        for (int i = 0; i < Ar.MaxCharge() + 1; i++)
        {
            Ar.Dump_2D(n, n, i);
            Ar.Dump_2D(Tn, Tn, i);
            if (i != Ar.MaxCharge())
            {
                Ar.Dump_2D(Sn, Ion, i);
                Ar.Dump_2D(Sn, CX, i);
                Ar.Dump_2D(Pra, Ion, i);
                Ar.Dump_2D(Pra, CX, i);
            }
            if (i != 0)
            {
                Ar.Dump_2D(Sn, Rec, i);
                Ar.Dump_2D(Pra, Rec, i);
            }
        }
    }

    /*string pathoutofpump = Outputpath + "pump.txt";
    ofstream outpump;
    outpump.open(pathoutofpump);
    outpump << D.name() << '\t' << D.Pump(0) << '\t' << D.Pump(1) << endl;
    outpump << D2.name() << '\t' << D2.Pump(0) << '\t' << D2.Pump(1) << endl;
    if (K_Ar)
        outpump << Ar.name() << '\t' << Ar.Pump(0) << '\t' << Ar.Pump(1) << endl;
    if (K_C || K_Methane)
        outpump << C.name() << '\t' << C.Pump(0) << '\t' << C.Pump(1) << endl;
    if (K_Methane)
    {
        outpump << CD4.name() << '\t' << CD4.Pump(0) << '\t' << CD4.Pump(1) << endl;
        outpump << CD3.name() << '\t' << CD3.Pump(0) << '\t' << CD3.Pump(1) << endl;
        outpump << CD2.name() << '\t' << CD2.Pump(0) << '\t' << CD2.Pump(1) << endl;
        outpump << CD1.name() << '\t' << CD1.Pump(0) << '\t' << CD1.Pump(1) << endl;
    }
    outpump.close();*/

    /*double Var_temp[N_POLOIDAL_GRID][N_RADIAL_GRID] = {0};
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
        {
            Var_temp[i][j] = D2.n_1(i, j);
            std::cout << D2.n_1(i, j) << '\t';
        }
        std::cout << endl;
    }
    Dump_2D_Global(Var_temp, "n_D2_1");*/
    // Dump_2D_Global(D.Ion_rate, "Ion_rate_D");
}

void Dump_2D_Global(double Var[N_POLOIDAL_GRID][N_RADIAL_GRID], string name)
{
    string pathout;
    pathout = Outputpath + "/data/" + name;
    ofstream out;
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << Var[i][j] << " ";
        out << endl;
    }
    out.close();
}

void HydrogenOutput(Particle *OutPar1, Particle *OutPar2)
{
    // Output the density of neutral particle
    OutPar1->Dump_2D(n, n, 0);  // density output
    OutPar1->Dump_2D(Tn, n, 0); // tempture output
    OutPar1->Dump_2D(E, n, 0);  // energy output
    OutPar1->Dump_2D(Sn, Ion);  // ion source output
    OutPar1->Dump_2D(Sn, CX);   // CX source output
    // OutPar1->Dump_array_2D(n, n, 0);

    OutPar1->Dump_2D(Smu, Ion);
    OutPar1->Dump_2D(Smu, CX);
    OutPar1->Dump_2D(Smu1, CX);
    OutPar1->Dump_2D(Smu2, CX);

    if (K_CX_DT == 1)
    {
        OutPar1->Dump_2D(Sn, CXDT); // CX source output
        OutPar1->Dump_2D(Smu, CXDT);
    }

    OutPar1->Dump_2D(SE, Ion);
    OutPar1->Dump_2D(SE, CX);
    OutPar1->Dump_2D(SE_n, CX);

    OutPar1->Dump_2D(Pra, Ion);
    OutPar1->Dump_2D(Pra, CX);

    if (K_Rec)
    {
        OutPar1->Dump_2D(Sn, Rec, 1); // rec source output
        OutPar1->Dump_2D(Smu, Rec, 1);
        OutPar1->Dump_2D(SE, Rec, 1);
        OutPar1->Dump_2D(Pra, Rec, 1);
    }
    OutPar1->Dump_Flux();

    OutPar2->Dump_2D(n, n, 0);
    OutPar2->Dump_2D(n, n, 1);
    // OutPar2->Dump_array_2D(n, n, 0);
    // OutPar2->Dump_array_2D(n, n, 1);

    OutPar2->Dump_2D(Tn, Tn, 0); // density output
    OutPar2->Dump_2D(Sn, Ion);
    OutPar2->Dump_2D(Sn, CX);
    OutPar2->Dump_2D(Sn, Ela);
    OutPar2->Dump_2D(SE, CX);
    OutPar2->Dump_2D(SE_n, CX);
    OutPar2->Dump_2D(SE, Ela);
    OutPar2->Dump_2D(Smu, CX);
    OutPar2->Dump_2D(Smu, Ela);
    OutPar2->Dump_2D(Sn, Diss1);
    OutPar2->Dump_2D(Sn, Diss2);
    OutPar2->Dump_Flux();
    for (int i = 0; i < 4; i++)
        OutPar2->Dump_2D(Sn, DS[i], 1);
    OutPar2->Dump_2D(Pra, DS[2], 1);

    string pathout;
    pathout = Outputpath + OutPar1->name() + "_0_" + "Vx";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + OutPar1->name() + "_0_" + "Vy";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + OutPar2->name() + "_0_" + "Vx";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar2->V_Grid(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + OutPar2->name() + "_0_" + "Vy";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar2->V_Grid(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + OutPar1->name() + "_0_" + "Lambda";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->lambda(i, j, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + OutPar2->name() + "_0_" + "Lambda";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar2->lambda(i, j, 1) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vx_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vy_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vz_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vx_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][0] << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][1] << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][2] << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vx_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();

    OutPar1->OutWallEro(1);
    if (K_Rec)
        OutPar1->OutWallEro(4);
    OutPar1->OutWallEro(-1);

    OutPar1->OutTargetEro(1);
    if (K_Rec)
        OutPar1->OutTargetEro(4);
    OutPar1->OutTargetEro(-1);

    OutPar1->OutPlasmaBoundaryEro(1);
    if (K_Rec)
        OutPar1->OutPlasmaBoundaryEro(4);
    OutPar1->OutPlasmaBoundaryEro(-1);
}

void HydrogenOutput_Tri(Particle *OutPar1, Particle *OutPar2)
{
    // Output the density of neutral particle
    OutPar1->Dump_Tri(n, n, 0);  // density output
    OutPar1->Dump_Tri(Tn, n, 0); // tempture output
    OutPar1->Dump_Tri(E, n, 0);  // energy output
    OutPar1->Dump_Tri(Sn, Ion);  // ion source output
    OutPar1->Dump_Tri(Sn, CX);   // CX source output
    // OutPar1->Dump_array_2D(n, n, 0);

    OutPar1->Dump_Tri(Smu, Ion);
    OutPar1->Dump_Tri(Smu, CX);
    OutPar1->Dump_Tri(Smu1, CX);
    OutPar1->Dump_Tri(Smu2, CX);

    if (K_CX_DT == 1)
    {
        OutPar1->Dump_Tri(Sn, CXDT); // CX source output
        OutPar1->Dump_Tri(Smu, CXDT);
    }

    OutPar1->Dump_Tri(SE, Ion);
    OutPar1->Dump_Tri(SE, CX);
    OutPar1->Dump_Tri(SE_n, CX);

    OutPar1->Dump_Tri(Pra, Ion);
    OutPar1->Dump_Tri(Pra, CX);
    OutPar1->Dump_Tri(Sn, R_with_H);
    OutPar1->Dump_Tri(Sn, R_with_H2);
    OutPar1->Dump_Tri(SE, R_with_H);
    OutPar1->Dump_Tri(SE, R_with_H2);
    OutPar1->Dump_Tri(Smu, R_with_H);
    OutPar1->Dump_Tri(Smu, R_with_H2);

    if (K_Rec)
    {
        OutPar1->Dump_Tri(Sn, Rec, 1); // rec source output
        OutPar1->Dump_Tri(Smu, Rec, 1);
        OutPar1->Dump_Tri(SE, Rec, 1);
        OutPar1->Dump_Tri(Pra, Rec, 1);
    }
    OutPar1->Dump_Flux();

    OutPar2->Dump_Tri(n, n, 0);
    OutPar2->Dump_Tri(n, n, 1);
    OutPar2->Dump_Tri(E, n, 0); // energy output
    // OutPar2->Dump_array_2D(n, n, 0);
    // OutPar2->Dump_array_2D(n, n, 1);

    OutPar2->Dump_Tri(Tn, Tn, 0); // density output
    OutPar2->Dump_Tri(Sn, Ion);
    OutPar2->Dump_Tri(Sn, CX);
    OutPar2->Dump_Tri(Sn, Ela);
    OutPar2->Dump_Tri(Sn, Diss1);
    OutPar2->Dump_Tri(Sn, Diss2);
    OutPar2->Dump_Tri(SE, CX);
    OutPar2->Dump_Tri(SE_n, CX);
    OutPar2->Dump_Tri(SE, Ela);
    OutPar2->Dump_Tri(Smu, CX);
    OutPar2->Dump_Tri(Smu, Ela);
    OutPar2->Dump_Tri(Sn, R_with_H);
    OutPar2->Dump_Tri(Sn, R_with_H2);
    OutPar2->Dump_Tri(SE, R_with_H);
    OutPar2->Dump_Tri(SE, R_with_H2);
    OutPar2->Dump_Tri(Smu, R_with_H);
    OutPar2->Dump_Tri(Smu, R_with_H2);
    // OutPar2->Dump_Flux();
    for (int i = 0; i < 4; i++)
        OutPar2->Dump_Tri(Sn, DS[i], 1);
    OutPar2->Dump_Tri(Pra, DS[2], 1);

    string pathout;
    if (K_NNCs)
    {
        pathout = Casepath + OutPar1->name() + "_0_" + "V_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar1->V_Grid_Tri(i, 0, 0) << "\t" << OutPar1->V_Grid_Tri(i, 1, 0) << "\t" << OutPar1->V_Grid_Tri(i, 2, 0) << endl;
        }
        out.close();

        pathout = Casepath + OutPar2->name() + "_0_" + "V_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar2->V_Grid_Tri(i, 0, 0) << "\t" << OutPar2->V_Grid_Tri(i, 1, 0) << "\t" << OutPar2->V_Grid_Tri(i, 2, 0) << endl;
        }
        out.close();

        pathout = Casepath + OutPar1->name() + "_0_" + "T_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar1->Tn_Tri(i, 0) << endl;
        }
        out.close();

        pathout = Casepath + OutPar2->name() + "_0_" + "T_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar2->Tn_Tri(i, 0) << endl;
        }
        out.close();

        pathout = Casepath + OutPar1->name() + "_0_" + "n_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar1->n_Tri(i, 0) << endl;
        }
        out.close();

        pathout = Casepath + OutPar2->name() + "_0_" + "n_Tri";
        out.open(pathout);
        for (int i = 0; i < Grid4.num_tris(); i++)
        {
            out << OutPar2->n_Tri(i, 0) << endl;
        }
        out.close();
    }

    pathout = Outputpath + OutPar1->name() + "_0_" + "T_Tri";
    out.open(pathout);
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out << OutPar1->Tn_Tri(i, 0) << endl;
    }
    out.close();

    pathout = Outputpath + OutPar2->name() + "_0_" + "T_Tri";
    out.open(pathout);
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out << OutPar2->Tn_Tri(i, 0) << endl;
    }
    out.close();

    pathout = Outputpath + OutPar1->name() + "_0_" + "V_Tri";
    out.open(pathout);
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out << OutPar1->V_Grid_Tri(i, 0, 0) << "\t" << OutPar1->V_Grid_Tri(i, 1, 0) << "\t" << OutPar1->V_Grid_Tri(i, 2, 0) << endl;
    }
    out.close();

    pathout = Outputpath + OutPar2->name() + "_0_" + "V_Tri";
    out.open(pathout);
    for (int i = 0; i < Grid4.num_tris(); i++)
    {
        out << OutPar2->V_Grid_Tri(i, 0, 0) << "\t" << OutPar2->V_Grid_Tri(i, 1, 0) << "\t" << OutPar2->V_Grid_Tri(i, 2, 0) << endl;
    }
    out.close();

    pathout = Outputpath + OutPar1->name() + "_0_" + "Lambda";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->lambda(i, j, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + OutPar2->name() + "_0_" + "Lambda";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar2->lambda(i, j, 1) << "\t";
        out << endl;
    }
    out.close();

    /*pathout = Outputpath + "Vx_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vx_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vy_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vz_D_0_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Be(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1_CX_Be";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Be(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_0_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Neu_Af(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1_CX_Af";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_Grid_CX_Ion_Af(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();*/

    /*pathout = Outputpath + "Vx_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][0] << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][1] << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_Ua";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << ua_D_1[i][j] * B[i][j][2] << "\t";
        out << endl;
    }
    out.close();

    pathout = Outputpath + "Vx_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 0, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vy_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 1, 0) << "\t";
        out << endl;
    }
    out.close();
    pathout = Outputpath + "Vz_D_1";
    out.open(pathout);
    for (int i = 0; i < N_poloidal; i++)
    {
        for (int j = 0; j < N_radial; j++)
            out << OutPar1->V_D_1(i, j, 2, 0) << "\t";
        out << endl;
    }
    out.close();*/
}