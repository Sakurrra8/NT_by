#include "Particle.h"

namespace
{
void accumulateDWallOrTargetImpact()
{
    const int surface_type = static_cast<int>(InterscePoint[0][4]);
    const int surface_index = static_cast<int>(InterscePoint[0][3]);
    int stats_index = -1;
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;

    if (surface_type == 11)
    {
        stats_index = surface_index;
        x1 = Grid1.Wall(surface_index, 0);
        y1 = Grid1.Wall(surface_index, 1);
        x2 = Grid1.Wall(surface_index + 1, 0);
        y2 = Grid1.Wall(surface_index + 1, 1);
    }
    else if (surface_type == 1 && MeshMode == 3)
    {
        stats_index = Grid1.Wall_num() + surface_index;
        x1 = Grid4.Target(surface_index, 0);
        y1 = Grid4.Target(surface_index, 1);
        x2 = Grid4.Target(surface_index, 2);
        y2 = Grid4.Target(surface_index, 3);
    }
    else
    {
        return;
    }

    D_W_sputtered.accumulate(
        stats_index, 1.5 * D.Tn(),
        cosIncidence3D(x1, y1, x2, y2, D.V(0), D.V(1), D.V(2)),
        D.Gamma());
}
}

void WallReflect()
{
    // coeff_recyc is recycling from Wall & target, coeff_reflect is from Trim database
    double coeff_recyc, coeff_reflect;
    coeff_recyc = coeff_ercyc_wall;
    if (InterscePoint[0][4] == 1)
        coeff_recyc = coeff_recyc_target;
    if (K_D2Flight && P->D2pOriginChannel() >= 1 && P->D2pOriginChannel() <= 2)
        D2.AuditD2pSecondaryDBoundaryHit(P->D2pOriginChannel(), (int)InterscePoint[0][4],
                                         P->Weight() * NumPar_now);

    if (MeshMode == 3)
    {
        if (K_Wallelement == 1)
        {
            if (P == &H)
                coeff_reflect = H_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &D)
            {
                const double angle = P->CalAngle((int)InterscePoint[0][3]);
                coeff_reflect = K_DWTrimReflection == 1
                                    ? (1.5 * P->Tn() >= DWTrimERMIN
                                           ? std::min(coeff_recyc, D_W_Trim.ReflectionProbability(1.5 * P->Tn(), angle))
                                           : 0.0)
                                    : D_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), angle);
            }
            else if (P == &T)
                coeff_reflect = T_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &H2 || P == &D2 || P == &T2)
                coeff_reflect = coeff_recyc;
        }
        else if (K_Wallelement == 2)
        {
            if (P == &H)
                coeff_reflect = H_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &D)
                coeff_reflect = D_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &T)
                coeff_reflect = T_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &H2 || P == &D2 || P == &T2)
                coeff_reflect = coeff_recyc;
        }

        if (P == &H)
        {
            if (Tools::Random() > coeff_recyc)
            {
                P->SetWeight(0.);
                return;
            }
            // P->AddWallEro(InterscePoint[0][3]);
            if (Tools::Random() < coeff_reflect / coeff_recyc)
            {
                P->setfate(0, 15, InterscePoint[0][4]);
                P->Init(3, 0);
                // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
            else
            {
                H2.Particlefrom(&H, 0.5, 0);
                P = &H2;
                P->setfate(0, 16, InterscePoint[0][4]);
                P->Init(3, 1);
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
        }
        else if (P == &D)
        {
            if (Tools::Random() > coeff_recyc)
            {
                P->SetWeight(0.);
                return;
            }
            if (InterscePoint[0][4] == 11)
            {
                P->AddWallEro(InterscePoint[0][3]);
                accumulateDWallOrTargetImpact();
            }
            else if (InterscePoint[0][4] == 1)
            {
                P->AddTargetEro(InterscePoint[0][3]);
                accumulateDWallOrTargetImpact();
            }

            if (Tools::Random() < coeff_reflect / coeff_recyc)
            {
                P->setfate(0, 15, InterscePoint[0][4]);
                P->Init(3, 0);
                // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
            else
            {
                D2.Particlefrom(&D, 0.5, 0);
                P = &D2;
                P->setfate(0, 16, InterscePoint[0][4]);
                P->Init(3, 1);
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
        }
        else if (P == &T)
        {
            if (Tools::Random() > coeff_recyc)
            {
                P->SetWeight(0.);
                return;
            }
            // P->AddWallEro(InterscePoint[0][3]);
            if (Tools::Random() < coeff_reflect / coeff_recyc)
            {
                P->setfate(0, 15, InterscePoint[0][4]);
                P->Init(3, 0);
                // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
            else
            {
                P = &T2;
                P->Particlefrom(&T, 0.5, 0);
                P->setfate(0, 16, InterscePoint[0][4]);
                P->Init(3, 1);
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
        }
        else if (P == &H2 || P == &D2 || P == &T2)
        {
            if (Tools::Random() > coeff_recyc)
            {
                P->SetWeight(0.);
                return;
            }
            P->setfate(0, 15, InterscePoint[0][4]);
            P->Init(3, 1);
            // P->track();
            if (P->Weight() == 0)
                return;
        }
        else
        {
            P->setfate(0, 15, InterscePoint[0][4]);
            P->Init(3, 0);
            // P->track();
            if (P->Weight() == 0)
                return;
        }
    }
    else if (MeshMode == 1)
    {
        /*if (K_Puff[0] == 1)
        {
            if (InterscePoint[0][3] == 25 && InterscePoint[0][4] == 11)
                if (InterscePoint[0][1] < Wall[25][0] && InterscePoint[0][1] > 1.73496001)
                {
                    // coeff_recyc = coeff_puff;
                    P->Pump_add(0, (1 - coeff_recyc) * P->Weight() * NumPar_now);
                }
        }
        if (K_Puff[1] == 1)
        {
            if (InterscePoint[0][3] == 35 && InterscePoint[0][4] == 11)
                if (InterscePoint[0][1] < 1.3924100219726538 && InterscePoint[0][1] > Wall[36][0])
                {
                    // coeff_recyc = coeff_puff;
                    P->Pump_add(1, (1 - coeff_recyc) * P->Weight() * NumPar_now);
                }
        }*/
        if (K_Puff[2] == 1)
            if (InterscePoint[0][4] == 11 && (InterscePoint[0][3] == 49 || InterscePoint[0][3] == 58))
            {
                coeff_recyc = coeff_puff;
                if (InterscePoint[0][3] == 58 && InterscePoint[0][4] == 11)
                {
                }
                // std::cout << "coeff_puff: " << coeff_puff << endl;
            }
        if (K_Puff[3] == 1)
        {
            if (InterscePoint[0][3] == 37 && InterscePoint[0][4] == 11)
            {
                coeff_recyc = coeff_puff;
                P->Pump_add(1, (1 - coeff_recyc) * P->Weight() * NumPar_now);
            }
        }
        if (K_Puff[4] == 1)
        {
            if (InterscePoint[0][3] == 53 && InterscePoint[0][4] == 11)
            {
                coeff_recyc = coeff_puff;
                P->Pump_add(1, (1 - coeff_recyc) * P->Weight() * NumPar_now);
            }
        }

        if (K_Wallelement == 1)
        {
            if (P == &H)
                coeff_reflect = H_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &D)
            {
                const double angle = P->CalAngle((int)InterscePoint[0][3]);
                coeff_reflect = K_DWTrimReflection == 1
                                    ? (1.5 * P->Tn() >= DWTrimERMIN
                                           ? std::min(coeff_recyc, D_W_Trim.ReflectionProbability(1.5 * P->Tn(), angle))
                                           : 0.0)
                                    : D_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), angle);
            }
            else if (P == &T)
                coeff_reflect = T_W.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &H2 || P == &D2 || P == &T2)
                coeff_reflect = coeff_recyc;
        }
        else if (K_Wallelement == 2)
        {
            if (P == &H)
                coeff_reflect = H_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &D)
                coeff_reflect = D_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &T)
                coeff_reflect = T_C.n_RefCoeff(K_Reflect, 1.5 * P->Tn(), P->CalAngle((int)InterscePoint[0][3]));
            else if (P == &H2 || P == &D2 || P == &T2)
                coeff_reflect = coeff_recyc;
        }

        // calculate the backscatting and particle recycling
        if (K_backScatter == 1)
        {
            if (P == &H)
            {
                if (Tools::Random() > coeff_recyc)
                {
                    P->SetWeight(0.);
                    return;
                }
                // P->AddWallEro(InterscePoint[0][3]);
                if (Tools::Random() < coeff_reflect / coeff_recyc)
                {
                    P->setfate(0, 15, InterscePoint[0][4]);
                    P->Init(3, 0);
                    // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
                else
                {
                    H2.Particlefrom(&H, 0.5, 0);
                    P = &H2;
                    P->setfate(0, 16, InterscePoint[0][4]);
                    P->Init(3, 1);
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
            }
            else if (P == &D)
            {
                if (Tools::Random() > coeff_recyc)
                {
                    P->SetWeight(0.);
                    return;
                }
                if (InterscePoint[0][4] == 11)
                {
                    P->AddWallEro(InterscePoint[0][3]);
                    accumulateDWallOrTargetImpact();
                }
                else if (InterscePoint[0][4] == 1)
                {
                    P->AddTargetEro(InterscePoint[0][3]);
                    accumulateDWallOrTargetImpact();
                }

                if (Tools::Random() < coeff_reflect / coeff_recyc)
                {
                    P->setfate(0, 15, InterscePoint[0][4]);
                    P->Init(3, 0);
                    // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
                else
                {
                    D2.Particlefrom(&D, 0.5, 0);
                    P = &D2;
                    P->setfate(0, 16, InterscePoint[0][4]);
                    P->Init(3, 1);
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
            }
            else if (P == &T)
            {
                if (Tools::Random() > coeff_recyc)
                {
                    P->SetWeight(0.);
                    return;
                }
                // P->AddWallEro(InterscePoint[0][3]);
                if (Tools::Random() < coeff_reflect / coeff_recyc)
                {
                    P->setfate(0, 15, InterscePoint[0][4]);
                    P->Init(3, 0);
                    // std::cout << P->V(0) << '\t' << P->V(1) << '\t' << P->V(2) << '\t' << P->ChargeTag() << endl;
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
                else
                {
                    P = &T2;
                    P->Particlefrom(&T, 0.5, 0);
                    P->setfate(0, 16, InterscePoint[0][4]);
                    P->Init(3, 1);
                    // P->track();
                    if (P->Weight() == 0)
                        return;
                }
            }
            else if (P == &H2 || P == &D2 || P == &T2)
            {
                if (Tools::Random() > coeff_recyc)
                {
                    P->SetWeight(0.);
                    return;
                }
                P->setfate(0, 15, InterscePoint[0][4]);
                P->Init(3, 1);
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
            else
            {
                P->setfate(0, 15, InterscePoint[0][4]);
                P->Init(3, 0);
                // P->track();
                if (P->Weight() == 0)
                    return;
            }
        }
        if (K_log)
        {
            if (P != &H && P != &H2 && P != &D && P != &D2 && P != &T && P != &T2)
            {
                Xlog << P->name() << '_' << P->Charge() << ' ' << P->X(0) << ", " << P->X(1) << ' ' << P->X(2) << endl;
                Vlog << P->name() << '_' << P->Charge() << ' ' << P->Vreal(0) << ' ' << P->Vreal(1) << ' ' << P->Vreal(2) << ' ' << P->Vreal(3) << ' ' << P->Vreal(7) << endl;
            }
        }
    }
}
