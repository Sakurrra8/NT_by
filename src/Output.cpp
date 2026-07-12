#include "Global.h"
#include "DBoundarySource.h"
#include "Particle.h"
#include <algorithm>
#include <map>
#include <set>

void HydrogenOutput(Particle *OutPar1, Particle *OutPar2);
void HydrogenOutput_Tri(Particle *OutPar1, Particle *OutPar2);
void SourceStratumSummaryOutput();
void EireneD2pReactionAuditOutput();
void D2pMesh3FlightAuditOutput();

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
        H.Dump_Flux();
        H2.Dump_Flux();

        Out_temp.open(Outputpath + "recycling_H.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_H_recyc[i] << " " << Tn_H_recyc[i] << " " << NumPar_H2_recyc[i] << endl;
        }
        Out_temp.close();
    }
    if (K_D)
    {
        D2.DumpD2pBalance_B2();
        D2.DumpD2pBalance_Tri();
        D2.DumpD2pPhysicsDecomposition_B2();
        if (MeshMode == 3 && K_D2Flight == 1)
        {
            D2.DumpD2pTrackLengthTri();
            D2.UseD2pTransportDensityForOutput();
        }
        HydrogenOutput(&D, &D2);
        if (MeshMode == 3)
            HydrogenOutput_Tri(&D, &D2);
        D.Dump_Flux();
        D2.Dump_Flux();
        EireneD2pReactionAuditOutput();
        D2pMesh3FlightAuditOutput();

        Out_temp.open(Outputpath + "recycling_D.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_D_recyc[i] << " " << Tn_D_recyc[i] << " " << NumPar_D2_recyc[i] << endl;
        }
        Out_temp.close();

        Out_temp.open(Outputpath + "target_incident_D.csv");
        Out_temp << "target,mean_incident_energy_eV,mean_incident_angle_deg,"
                    "mean_fast_probability,mean_reflected_energy_eV,D_source_s-1,D2_source_s-1\n";
        const int target_count = std::min(
            2 * N_radial, static_cast<int>(DTargetIncidentAngle.size()));
        for (int i = 0; i < target_count; ++i)
        {
            Out_temp << i << ','
                     << DTargetMeanIncidentEnergy[i] << ','
                     << DTargetIncidentAngle[i] << ','
                     << DTargetFastProbability[i] << ','
                     << DTargetMeanReflectedEnergy[i] << ','
                     << NumPar_D_recyc[i] << ','
                     << NumPar_D2_recyc[i] << '\n';
        }
        Out_temp.close();
        if (K_DBoundarySource)
            D_BoundarySource.WriteSummary(Outputpath + "D_boundary_source.csv");
    }
    if (K_T)
    {
        HydrogenOutput(&T, &T2);
        if (MeshMode == 3)
            HydrogenOutput_Tri(&T, &T2);
        T.Dump_Flux();
        T2.Dump_Flux();

        Out_temp.open(Outputpath + "recycling_T.txt");
        for (int i = 0; i < 76; i++)
        {
            Out_temp << NumPar_T_recyc[i] << " " << Tn_T_recyc[i] << " " << NumPar_T2_recyc[i] << endl;
        }
        Out_temp.close();
    }
    SourceStratumSummaryOutput();
    if (K_H5Output)
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

namespace
{
struct EireneReactionAuditRow
{
    int reaction_id = 0;
    string database;
    string table;
    string reaction_number;
    string eirene_type;
};

bool ParseEireneD2pReactionAudit(const string &path,
                                std::map<int, EireneReactionAuditRow> &reactions,
                                std::set<int> &d2p_reaction_ids)
{
    ifstream input(path);
    if (!input.is_open())
        return false;

    bool in_reaction_list = false;
    bool in_d2p_species = false;
    string line;
    while (std::getline(input, line))
    {
        if (line.find("**  Reactions") != string::npos)
        {
            in_reaction_list = true;
            continue;
        }
        if (in_reaction_list && line.find("** 4a.") != string::npos)
        {
            in_reaction_list = false;
            continue;
        }
        if (line.find("D2+") != string::npos && line.find("**") == string::npos)
        {
            in_d2p_species = true;
            continue;
        }
        if (in_d2p_species && line.find("*** 5.") != string::npos)
        {
            in_d2p_species = false;
            continue;
        }

        if (in_reaction_list)
        {
            EireneReactionAuditRow row;
            std::istringstream fields(line);
            if (fields >> row.reaction_id >> row.database >> row.table >>
                    row.reaction_number >> row.eirene_type)
            {
                reactions[row.reaction_id] = row;
            }
        }
        else if (in_d2p_species)
        {
            int reaction_id = 0;
            int code1 = 0;
            int code2 = 0;
            std::istringstream fields(line);
            if (fields >> reaction_id >> code1 >> code2)
                d2p_reaction_ids.insert(reaction_id);
        }
    }
    return true;
}
}

void EireneD2pReactionAuditOutput()
{
    const std::vector<string> candidates = {
        Casepath + "input.dat",
        Casepath + "solps_output/input.dat",
        Outputpath + "input.dat",
        "input.dat"};

    std::map<int, EireneReactionAuditRow> reactions;
    std::set<int> d2p_reaction_ids;
    string input_path;
    for (const string &candidate : candidates)
    {
        reactions.clear();
        d2p_reaction_ids.clear();
        if (ParseEireneD2pReactionAudit(candidate, reactions, d2p_reaction_ids))
        {
            input_path = candidate;
            break;
        }
    }

    ofstream readme(Outputpath + "EIRENE_D2p_reaction_audit_readme.txt");
    if (input_path.empty())
    {
        readme << "No EIRENE input.dat was found. No reaction audit CSV was generated.\n"
               << "Searched Casepath/input.dat, Casepath/solps_output/input.dat, Outputpath/input.dat, and ./input.dat.\n"
               << "Copy or link the current case input.dat to one of those read-only locations to enable the audit.\n"
               << "Do not infer that AMJUEL H.4 2.2.12 (DS2) is absent from this missing-file condition.\n";
        return;
    }

    ofstream out(Outputpath + "EIRENE_D2p_reaction_audit.csv");
    out << "reaction_id,database,table,reaction_number,eirene_type,attached_to_D2p\n";
    for (int reaction_id : d2p_reaction_ids)
    {
        const auto reaction = reactions.find(reaction_id);
        if (reaction == reactions.end())
            continue;
        const EireneReactionAuditRow &row = reaction->second;
        out << row.reaction_id << ',' << row.database << ',' << row.table << ','
            << row.reaction_number << ',' << row.eirene_type << ",1\n";
    }

    const auto reaction13 = reactions.find(13);
    const bool reaction13_attached = d2p_reaction_ids.count(13) != 0;
    const bool reaction13_is_2p2p12 =
        reaction13 != reactions.end() &&
        reaction13->second.database == "AMJUEL" &&
        reaction13->second.table == "H.4" &&
        reaction13->second.reaction_number == "2.2.12";
    readme << "Parsed EIRENE input: " << input_path << '\n'
           << "DS2 is present in EIRENE input as reaction 13: "
           << ((reaction13_attached && reaction13_is_2p2p12) ? "confirmed" : "not confirmed")
           << ".\n"
           << "The audit is read-only; EIRENE was not run and input.dat was not modified.\n";
}

void D2pMesh3FlightAuditOutput()
{
    ofstream audit(Outputpath + "D2p_mesh3_flight_audit.csv");
    const bool transport_enabled = MeshMode == 3 && K_D2Flight == 1;
    audit << "MeshMode,K_D2Flight,K_flight,K_Prob,K_NNCs,EireneHeavyEnergyScale,"
          << "push_tri_allows_charged,track_tri_has_charge_gt_0_branch,"
          << "coll_immediate_charge1_branch_exists,official_Tri_D2p_density_estimator,"
          << "mesh3_D2p_flight_status,D2p_max_allowed_speed_m_s\n";
    audit << MeshMode << ',' << K_D2Flight << ',' << K_flight << ',' << K_Prob << ','
          << K_NNCs << ',' << EireneHeavyEnergyScale << ','
          << "1,1,1,"
          << (transport_enabled ? "transport_track_length" : "local_balance_P_over_L") << ','
          << (transport_enabled ? "transport_enabled" : "transport_disabled")
          << ',' << D2pMaxAllowedSpeed << '\n';

    ofstream interpretation(Outputpath + "D2p_tri_density_interpretation.txt");
    interpretation
        << "Triangular-grid D2+ local-balance definition:\n"
        << "Tri_n_[i][1] = (Ion_D2 + CX_D2 + optional CXDT_D2) / "
        << "(triVolume * ne * (k_2.2.11 + k_2.2.12 + k_2.2.14))\n\n"
        << "When K_D2Flight=0, n_D2_1 and n_D2_1_Tri use this local balance estimate.\n"
        << "When K_D2Flight=1, n_D2_1 and n_D2_1_Tri are overwritten with the "
        << "D2+ track-length transport density before output.\n"
        << "The local-balance density remains available in D2p_balance_Tri.csv and "
        << "D2p_track_length_Tri.csv.\n";

    ofstream design(Outputpath + "D2p_mesh3_test_ion_tracking_design.txt");
    design
        << "K_D2Flight controls charged D2+ transport in MeshMode=3.\n"
        << "K_D2Flight=0: D2+ is not transported as a test ion; the formal output remains the local-balance density P/L.\n"
        << "K_D2Flight=1: D2 ionization/CX creates a charged D2+ test ion and follows it with fixed-time steps.\n\n"
        << "Transport principle:\n"
        << "1. New D2+ inherits the parent D2 neutral velocity projected onto the local magnetic field.\n"
        << "2. Each charged step applies the parallel ion friction/thermal relaxation and optional electric/drift terms.\n"
        << "3. The step is limited by the current triangle size so a fixed-time charged step cannot skip many mesh cells.\n"
        << "4. Triangle crossing updates Tri_Index_/XY_ before the next force evaluation.\n"
        << "5. DS1/DS2/DS3 are sampled after finite D2+ residence time using the local D2+ reaction rates.\n"
        << "6. The sampled DS event is tallied as a source/fate contribution and terminates the current D2+.\n"
        << "7. DS1 creates two ion products and no neutral test particle is launched.\n"
        << "8. DS2 launches one neutral D, and DS3 launches neutral D with twice the event weight.\n"
        << "9. DS2 secondary D uses 4.3 eV and DS3 secondary D uses 16 eV, matching the EIRENE reaction cards.\n"
        << "10. The secondary neutral D energy/temperature is recalculated with the D mass and then tracked by the normal neutral-D chain until ionization or another terminal fate.\n"
        << "11. D+ products are accounted in the source/fate tallies; they are not launched as separate test ions.\n"
        << "12. Leaving the plasma triangle mesh or producing non-finite/unphysical speed terminates that test ion as a boundary/numerical loss.\n\n"
        << "Speed guard:\n"
        << "D2pMaxAllowedSpeed=" << D2pMaxAllowedSpeed << " m/s.\n"
        << "This is a numerical sanity limit, not a physical velocity clamp for normal D2+ ions.\n"
        << "For D2+, 1e6 m/s is about 20 keV, far above the tens-of-eV edge-ion scale expected here.\n"
        << "Normal D2+ speeds in this case should be O(1e3-1e4) m/s, well below the guard.\n"
        << "The guard only catches corrupted states such as invalid mesh-field access or runaway arithmetic.\n";

    ofstream status(Outputpath + "D2p_transport_status.txt");
    status
        << "MeshMode=" << MeshMode << '\n'
        << "K_D2Flight=" << K_D2Flight << '\n'
        << "K_NNCs=" << K_NNCs << '\n'
        << "EireneHeavyEnergyScale=" << EireneHeavyEnergyScale << '\n'
        << "hydrogen_isotope_scaling=automatic_by_background_species\n"
        << "D2p_transport=" << (transport_enabled ? "enabled" : "disabled") << '\n'
        << "D2p_secondary_product_transport=neutral_D_from_DS2_DS3_tracked\n"
        << "D2p_max_allowed_speed_m_s=" << D2pMaxAllowedSpeed << '\n'
        << "density_reference="
        << (transport_enabled ? "transport_track_length" : "local_balance_P_over_L")
        << '\n'
        << "track_length_outputs="
        << (transport_enabled ? "enabled" : "disabled") << '\n';

    if (transport_enabled)
    {
        ofstream warning(Outputpath + "D2p_mesh3_flight_warning.txt");
        warning
            << "MeshMode=3 with K_D2Flight=1 enables D2+ test-ion tracking.\n"
            << "Particle::track() follows charged D2+ with fixed-time steps and triangle-index updates.\n"
            << "n_D2_1 and n_D2_1_Tri are D2+ track-length transport densities in this mode.\n"
            << "Use D2p_balance_Tri.csv for the local P/L balance diagnostic.\n";
    }
}

void SourceStratumSummaryOutput()
{
    ofstream out(Outputpath + "source_stratum_summary.csv");
    out << std::setprecision(17);
    out << "species,charge,source_stratum,launched_weight_s-1,launched_events,b2_track_length_inventory\n";
    if (K_H)
    {
        H.AppendSourceStratumSummary(out);
        H2.AppendSourceStratumSummary(out);
    }
    if (K_D)
    {
        D.AppendSourceStratumSummary(out);
        D2.AppendSourceStratumSummary(out);
    }
    if (K_T)
    {
        T.AppendSourceStratumSummary(out);
        T2.AppendSourceStratumSummary(out);
    }
    if (K_Methane)
    {
        CD4.AppendSourceStratumSummary(out);
        CD3.AppendSourceStratumSummary(out);
        CD2.AppendSourceStratumSummary(out);
        CD1.AppendSourceStratumSummary(out);
    }
    if (K_C || K_Methane)
        C.AppendSourceStratumSummary(out);
    if (K_Ar)
        Ar.AppendSourceStratumSummary(out);

    ofstream readme(Outputpath + "source_stratum_readme.txt");
    readme << "source_stratum_summary.csv reports launched/re-emitted particle weights, event counts, and B2 track-length inventory.\n"
           << "b2_track_length_inventory is the volume integral of density attributed to each primary source stratum.\n"
           << "It can be compared with EIRENE pdena_int_b2 and pdenm_int_b2 after matching source strata.\n"
           << "The launch columns are not collisional source-term decompositions; those require Sn_by_stratum diagnostics.\n";
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

    if (K_H5Output)
    {
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
    for (int i = 0; i < 4; i++)
        OutPar2->Dump_Tri(Sn, DS[i], 1);
    OutPar2->Dump_Tri(Pra, DS[2], 1);

    string pathout;
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

    OutPar1->WriteTargetImpactSummary(
        Outputpath + "target_impact_" + OutPar1->name() + ".csv");
    OutPar2->WriteTargetImpactSummary(
        Outputpath + "target_impact_" + OutPar2->name() + ".csv");

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
