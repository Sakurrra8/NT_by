#include "Global.h"

void Initialize(int Input, char *settingfile[])
{
  string Inputstring, line;
  double number_temp;
  std::ifstream In;
  In.open(settingfile[1], std::ios::in); // ios::in means read
  if (!In.is_open())
  {
    std::cerr << "This file READING for setting.log have some problem!!!\n";
  }
  // std::cout << settingfile << " " << "Successfully read" << endl;
  std::getline(In, line);
  In >> Inputstring >> Inputpath;
  std::getline(In, line);
  In >> Inputstring >> Casepath;
  std::getline(In, line);
  In >> Inputstring >> Outputpath;

  if (Input == 3)
  {
    Outputpath = settingfile[2];
  }

  Databasepath = Inputpath + "database/";
  std::getline(In, line);
  std::streampos after_output = In.tellg();
  if (In >> Inputstring)
  {
    if (Inputstring == "Databasepath")
    {
      In >> Databasepath;
      std::getline(In, line);
    }
    else
    {
      In.clear();
      In.seekg(after_output);
    }
  }
  In >> Inputstring >> N_poloidal;
  std::getline(In, line);
  In >> Inputstring >> N_radial;
  /// read the N_Trimesh

  const int P = N_poloidal;
  const int R = N_radial;
  const int T = 8010;

  // 三维：Grid[P][R][8]，B[P][R][4]
  Grid.assign(P,
              std::vector<std::vector<double>>(R, std::vector<double>(8, 0.0)));
  B.assign(P, std::vector<std::vector<double>>(R, std::vector<double>(4, 0.0)));

  // 二维：每个都是 [P][R]
  auto make2D = [&](std::vector<std::vector<double>> &a)
  {
    a.assign(P, std::vector<double>(R, 0.0));
  };

  // 一维：每个都是 [T]
  auto make1D = [&](std::vector<double> &a)
  { a.assign(T, double(0.0)); };

  auto make1D_V = [&](std::vector<std::vector<double>> &a)
  {
    a.assign(T, std::vector<double>(3, 0.0));
  };

  make2D(ne);
  make2D(n_D_1);
  make2D(Te);
  make2D(Ti);
  // make2D(Te_coll);
  // make2D(Ti_coll);

  make2D(ua_D_1);
  make2D(Ti_D_thermal);

  make1D(n_D_0_Tri);
  make1D(T_D2_0_Tri);
  make1D(T_D_0_Tri);
  make1D(n_D2_0_Tri);
  make1D_V(ua_D_0_Tri);
  make1D_V(ua_D2_0_Tri);
  make2D(n_D_0);
  make2D(T_D2_0);
  make2D(T_D_0);
  make2D(n_D2_0);
  make2D(n_D2_1);
  make2D(Vi_D);

  make2D(Volume);
  make2D(Area);
  make2D(Btot);
  make2D(Btor);
  make2D(Brad);
  make2D(Bpol);
  make2D(Erad);
  make2D(Epol);
  make2D(Etor);
  make2D(vdEp);
  make2D(vdEr);

  NumPar_H_recyc.resize(N_radial * 2);
  NumPar_H2_recyc.resize(N_radial * 2);
  Tn_H_recyc.resize(N_radial * 2);
  NumPar_D_recyc.resize(N_radial * 2);
  NumPar_D2_recyc.resize(N_radial * 2);
  Tn_D_recyc.resize(N_radial * 2);
  DTargetIncidentAngle.resize(N_radial * 2, 60.0);
  DTargetFastProbability.resize(N_radial * 2, 0.0);
  DTargetMeanIncidentEnergy.resize(N_radial * 2, 0.0);
  DTargetMeanReflectedEnergy.resize(N_radial * 2, 0.0);
  NumPar_T_recyc.resize(N_radial * 2);
  NumPar_T2_recyc.resize(N_radial * 2);
  Tn_T_recyc.resize(N_radial * 2);

  std::getline(In, line);
  In >> Inputstring >> K_log;
  // std::cout << K_log << endl;
  std::getline(In, line);
  In >> Inputstring >> StepLog;
  // std::cout << StepLog << endl;
  std::getline(In, line);
  In >> Inputstring >> K_H2_elastic;
  // std::cout << K_H2_elastic << endl;
  std::getline(In, line);
  In >> Inputstring >> K_EcrossBDrift;
  // std::cout << K_EcrossBDrift << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Recyc;
  // std::cout << K_Recyc << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Rec;
  // std::cout << K_Rec << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Methane;
  // std::cout << K_Methane << endl;
  std::getline(In, line);
  In >> Inputstring >> K_C;
  // std::cout << K_C << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Ar;
  // std::cout << K_Ar << endl;
  std::getline(In, line);
  In >> Inputstring >> K_database;
  // std::cout << K_database << endl;
  std::getline(In, line);
  In >> Inputstring >> K_database_Pra;
  // std::cout << K_database_Pra << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Maxwell;
  // std::cout << K_Maxwell << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Ei;
  // std::cout << K_Ei << endl;
  std::getline(In, line);
  In >> Inputstring >> K_WallRefl;
  // std::cout << K_WallRefl << endl;
  std::getline(In, line);
  In >> Inputstring >> K_backScatter;
  // std::cout << K_backScatter << endl;
  std::getline(In, line);
  In >> Inputstring >> K_MarColl;
  // std::cout << K_MarColl << endl;
  std::getline(In, line);
  In >> Inputstring >> K_flight;
  // std::cout << K_flight << endl;
  std::getline(In, line);
  In >> Inputstring >> K_CX_impurity;
  // std::cout << K_CX_impurity << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Puff[0] >> K_Puff[1] >> K_Puff[2] >> K_Puff[3] >>
      K_Puff[4];
  // std::cout << K_Puff[0] << '\t' << K_Puff[1] << endl;
  std::getline(In, line);
  In >> Inputstring >> K_D2Flight;
  // std::cout << K_D2Flight << endl;
  std::getline(In, line);
  In >> Inputstring >> K_abnormal_transport;
  // std::cout << K_abnormal_transport << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Reflect;
  // std::cout << K_Reflect << endl;
  std::getline(In, line);
  In >> Inputstring >> K_ReflectDirection;
  // std::cout << K_ReflectDirection << endl;
  std::getline(In, line);
  In >> Inputstring >> backGridBoundry;
  std::getline(In, line);
  In >> Inputstring >> K_Prob;
  // std::cout << K_Prob << endl;
  std::getline(In, line);
  In >> Inputstring >> K_H;
  // std::cout << K_H << endl;
  std::getline(In, line);
  In >> Inputstring >> K_D;
  // std::cout << K_D << endl;
  std::getline(In, line);
  In >> Inputstring >> K_T;
  // std::cout << K_T << endl;
  std::getline(In, line);
  In >> Inputstring >> K_back;
  // std::cout << K_back << endl;

  // Optional named switch. A second pass keeps all existing positional setting
  // files backward compatible when this key is absent.
  {
    std::ifstream option_stream(settingfile[1]);
    std::string option_name;
    std::string option_line;
    while (std::getline(option_stream, option_line))
    {
      std::istringstream option_input(option_line);
      if (!(option_input >> option_name))
        continue;
      if (option_name == "K_DWTrimReflection")
        option_input >> K_DWTrimReflection;
      else if (option_name == "SurfaceTemperature_eV")
        option_input >> T_wall;
      else if (option_name == "DWTrimERMIN_eV")
        option_input >> DWTrimERMIN;
      else if (option_name == "DWTrimDatabase")
        option_input >> DWTrimDatabase;
      else if (option_name == "K_DWTargetActualAngle")
        option_input >> K_DWTargetActualAngle;
      else if (option_name == "K_DTargetSourceMode")
        option_input >> K_DTargetSourceMode;
      else if (option_name == "DTargetIncidentModel")
        option_input >> DTargetIncidentModel;
      else if (option_name == "DTargetSheathFactor")
        option_input >> DTargetSheathFactor;
      else if (option_name == "DTargetIncidentSamples")
        option_input >> DTargetIncidentSamples;
      else if (option_name == "K_DBoundarySource")
        option_input >> K_DBoundarySource;
      else if (option_name == "DBoundaryLaunchModel")
        option_input >> DBoundaryLaunchModel;
      else if (option_name == "numPar_flight_DBoundary")
        option_input >> numPar_flight_DBoundary;
      else if (option_name == "DBoundaryFluxFile")
        option_input >> DBoundaryFluxFile;
      else if (option_name == "K_Roulette")
        option_input >> K_Roulette;
      else if (option_name == "K_Splitting")
        option_input >> K_Splitting;
      else if (option_name == "K_H5Output")
        option_input >> K_H5Output;
      else if (option_name == "K_NNCs")
        option_input >> K_NNCs;
      else if (option_name == "D2ElasticModel")
        option_input >> D2ElasticModel;
      else if (option_name == "EireneRateArgumentScale")
      {
        option_input >> EireneRateArgumentScale;
        EireneElectronTemperatureScale = EireneRateArgumentScale;
        EireneHeavyEnergyScale = EireneRateArgumentScale;
      }
      else if (option_name == "EireneElectronTemperatureScale")
        option_input >> EireneElectronTemperatureScale;
      else if (option_name == "EireneHeavyEnergyScale")
        option_input >> EireneHeavyEnergyScale;
      else if (option_name == "W_RouletteMin")
        option_input >> W_RouletteMin;
      else if (option_name == "W_RouletteTarget")
        option_input >> W_RouletteTarget;
      else if (option_name == "W_SplitMax")
        option_input >> W_SplitMax;
      else if (option_name == "W_SplitTarget")
        option_input >> W_SplitTarget;
      else if (option_name == "W_SplitMin")
        option_input >> W_SplitMin;
      else if (option_name == "MaxSplit")
        option_input >> MaxSplit;
      else if (option_name == "MaxSplitDepth")
        option_input >> MaxSplitDepth;
      else if (option_name == "ImportanceOutside")
        option_input >> RegionImportance[0];
      else if (option_name == "ImportanceDivertor")
        option_input >> RegionImportance[1];
      else if (option_name == "ImportanceMainChamber")
        option_input >> RegionImportance[2];
      else if (option_name == "ImportanceMainPoloidalBegin")
        option_input >> ImportanceMainPoloidalBegin;
      else if (option_name == "ImportanceMainPoloidalEnd")
        option_input >> ImportanceMainPoloidalEnd;
    }
  }
  D2ElasticModel = std::clamp(D2ElasticModel, 0, 3);
  DTargetIncidentModel = std::clamp(DTargetIncidentModel, 0, 2);
  DTargetIncidentSamples = std::max(256, DTargetIncidentSamples);
  K_DBoundarySource = std::clamp(K_DBoundarySource, 0, 1);
  DBoundaryLaunchModel = std::clamp(DBoundaryLaunchModel, 0, 1);
  numPar_flight_DBoundary = std::max(0, numPar_flight_DBoundary);
  if (!K_H2_elastic)
    D2ElasticModel = 0;
  std::cout << "D2ElasticModel: " << D2ElasticModel << std::endl;
  std::cout << "DTargetIncidentModel: " << DTargetIncidentModel
            << (DTargetIncidentModel == 1
                    ? " (EIRENE NEMODS=7 drifting flux plus sheath"
                    : DTargetIncidentModel == 2
                          ? " (normal monoenergetic sensitivity, target 3Ti+0.5Te+sheath"
                          : " (legacy mean energy/angle")
            << ", sheath="
            << (DTargetSheathFactor <= 0. ? "EIRENE dynamic" : "fixed")
            << ", factor=" << DTargetSheathFactor
            << ", averaging samples=" << DTargetIncidentSamples << ")"
            << std::endl;
  std::cout << "D boundary FNIY source: "
            << (K_DBoundarySource ? "on" : "off")
            << " (launch model=" << DBoundaryLaunchModel
            << (DBoundaryLaunchModel == 0 ? ", EIRENE surface recycling" :
                                             ", direct outward D sensitivity")
            << ", histories=" << numPar_flight_DBoundary
            << ", file=" << DBoundaryFluxFile << ")" << std::endl;
  std::getline(In, line);
  In >> Inputstring >> K_dn;
  // std::cout << K_dn << endl;
  std::getline(In, line);
  In >> Inputstring >> Ratio_T;
  // std::cout << Ratio_T << endl;
  std::getline(In, line);
  In >> Inputstring >> K_CX_DT;
  // std::cout << K_CX_DT << endl;
  std::getline(In, line);
  In >> Inputstring >> K_DT;
  // std::cout << K_DT << endl;
  std::getline(In, line);
  In >> Inputstring >> K_Wallelement;
  // std::cout << K_Wallelement << endl;

  std::getline(In, line);
  In >> Inputstring >> dt;
  // std::cout << dt << endl;
  std::getline(In, line);
  In >> Inputstring >> t_max;
  // std::cout << t_max << endl;
  std::getline(In, line);
  In >> Inputstring >> numPar_flight;
  // std::cout << numPar_flight << endl;
  std::getline(In, line);
  In >> Inputstring >> numPar_flight_D2;
  // std::cout << numPar_flight_D2 << endl;
  std::getline(In, line);
  In >> Inputstring >> numPar_flight_T2;
  // std::cout << numPar_flight_T2 << endl;
  std::getline(In, line);
  In >> Inputstring >> numPar_flight_CD4;
  // std::cout << numPar_flight_CD4 << endl;
  std::getline(In, line);
  In >> Inputstring >> numPar_flight_Target;
  // std::cout << numPar_flight_Target << endl;
  std::getline(In, line);
  In >> Inputstring >> Num_D2_pump;
  // std::cout << Num_D2_pump << endl;
  std::getline(In, line);
  In >> Inputstring >> Num_T2_pump;
  // std::cout << Num_T2_pump << endl;
  std::getline(In, line);
  In >> Inputstring >> Num_CD4_pump;
  // std::cout << Num_CD4_pump << endl;
  std::getline(In, line);
  K_Pump = ((K_D && Num_D2_pump > 0. && numPar_flight_D2 > 0) ||
            (K_T && Num_T2_pump > 0. && numPar_flight_T2 > 0) ||
            (K_Methane && Num_CD4_pump > 0. && numPar_flight_CD4 > 0));
  In >> Inputstring >> coeff_recyc_target;
  // std::cout << coeff_recyc_target << endl;
  std::getline(In, line);
  In >> Inputstring >> coeff_ercyc_wall;
  // std::cout << coeff_ercyc_wall << endl;
  std::getline(In, line);
  In >> Inputstring >> coeff_puff;
  // std::cout << coeff_puff << endl;

  // Optional switches. Old setting files can omit these lines.
  bool mesh_mode_read = false;
  while (std::getline(In, line))
  {
    if (line.empty())
      continue;
    std::istringstream iss(line);
    iss >> Inputstring;
    if (Inputstring.empty() || Inputstring[0] == '#')
      continue;
    if (Inputstring == "MeshMode")
    {
      iss >> MeshMode;
      mesh_mode_read = true;
      break;
    }
    else if (Inputstring == "K_Roulette")
      iss >> K_Roulette;
    else if (Inputstring == "K_Splitting")
      iss >> K_Splitting;
    else if (Inputstring == "K_H5Output")
      iss >> K_H5Output;
    else if (Inputstring == "K_NNCs")
      iss >> K_NNCs;
    else if (Inputstring == "D2ElasticModel")
      iss >> D2ElasticModel;
    else if (Inputstring == "DTargetIncidentModel")
      iss >> DTargetIncidentModel;
    else if (Inputstring == "DTargetSheathFactor")
      iss >> DTargetSheathFactor;
    else if (Inputstring == "DTargetIncidentSamples")
      iss >> DTargetIncidentSamples;
    else if (Inputstring == "K_DBoundarySource")
      iss >> K_DBoundarySource;
    else if (Inputstring == "DBoundaryLaunchModel")
      iss >> DBoundaryLaunchModel;
    else if (Inputstring == "numPar_flight_DBoundary")
      iss >> numPar_flight_DBoundary;
    else if (Inputstring == "DBoundaryFluxFile")
      iss >> DBoundaryFluxFile;
    else if (Inputstring == "W_RouletteMin")
      iss >> W_RouletteMin;
    else if (Inputstring == "W_RouletteTarget")
      iss >> W_RouletteTarget;
    else if (Inputstring == "W_SplitMax")
      iss >> W_SplitMax;
    else if (Inputstring == "W_SplitTarget")
      iss >> W_SplitTarget;
    else if (Inputstring == "W_SplitMin")
      iss >> W_SplitMin;
    else if (Inputstring == "MaxSplit")
      iss >> MaxSplit;
    else if (Inputstring == "MaxSplitDepth")
      iss >> MaxSplitDepth;
    else if (Inputstring == "ImportanceOutside")
      iss >> RegionImportance[0];
    else if (Inputstring == "ImportanceDivertor")
      iss >> RegionImportance[1];
    else if (Inputstring == "ImportanceMainChamber")
      iss >> RegionImportance[2];
    else if (Inputstring == "ImportanceMainPoloidalBegin")
      iss >> ImportanceMainPoloidalBegin;
    else if (Inputstring == "ImportanceMainPoloidalEnd")
      iss >> ImportanceMainPoloidalEnd;
    else if (Inputstring == "EireneRateArgumentScale")
    {
      iss >> EireneRateArgumentScale;
      EireneElectronTemperatureScale = EireneRateArgumentScale;
      EireneHeavyEnergyScale = EireneRateArgumentScale;
    }
    else if (Inputstring == "EireneElectronTemperatureScale")
      iss >> EireneElectronTemperatureScale;
    else if (Inputstring == "EireneHeavyEnergyScale")
      iss >> EireneHeavyEnergyScale;
  }
  if (!mesh_mode_read)
    MeshMode = 1;

  /// Extern Mode
  // std::cout << MeshMode << endl;
  if (MeshMode == 2)
  {
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_nx(number_temp);
    // std::cout << Grid3.nx() << endl;
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_ny(number_temp);
    // std::cout << Grid3.ny() << endl;
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_nCv(number_temp);
    // std::cout << Grid3.nCv() << endl;
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_nVx(number_temp);
    // std::cout << Grid3.nVx() << endl;
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_nFc(number_temp);
    // std::cout << Grid3.nFc() << endl;
    std::getline(In, line);
    In >> Inputstring >> number_temp;
    Grid3.Set_nCmxVx(number_temp);
    // std::cout << Grid3.nCmxVx() << endl;
  }
  else
  {
    Grid3.Set_nx(0);
    Grid3.Set_ny(0);
    Grid3.Set_nCv(0);
    Grid3.Set_nVx(0);
    Grid3.Set_nFc(0);
    Grid3.Set_nCmxVx(0);
  }

  In.close();

  if (K_H)
  {
    make2D(n_H_1);
    make2D(ua_H_1);
    make2D(Ti_H_thermal);
    make1D(n_H_0_Tri);
    make1D(T_H2_0_Tri);
    make1D(T_H_0_Tri);
    make1D(n_H2_0_Tri);
    make1D_V(ua_H_0_Tri);
    make1D_V(ua_H2_0_Tri);
    make2D(n_H_0);
    make2D(T_H2_0);
    make2D(T_H_0);
    make2D(n_H2_0);
    make2D(n_H2_1);
    make2D(Vi_H);
  }

  if (K_T)
  {
    make2D(n_T_1);
    make2D(ua_T_1);
    make2D(Ti_T_thermal);
    make1D(n_T_0_Tri);
    make1D(T_T2_0_Tri);
    make1D(T_T_0_Tri);
    make1D(n_T2_0_Tri);
    make1D_V(ua_T_0_Tri);
    make1D_V(ua_T2_0_Tri);
    make2D(n_T_0);
    make2D(T_T2_0);
    make2D(T_T_0);
    make2D(n_T2_0);
    make2D(n_T2_1);
    make2D(Vi_T);
  }

  V_1.assign(3, 0.0);
  V_2.assign(3, 0.0);
}
