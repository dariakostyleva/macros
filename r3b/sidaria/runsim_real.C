void runsim_real(Int_t nEvents = 0)
{

  // =========== Configuration area =============================

  //TString OutFile = "/data.local2/G4_sim_momenta/real_case/sim_out_real_air_500k_041019.root";  // Output file for data
 // TString OutFile = "sim_out_real_thesis_10k.root";
  TString OutFile = "sim_out_real.root";
  TString ParFile = "sim_par.root";  // Output file for params

  Bool_t fVis = true;                // Store tracks for visualization
  Bool_t fUserPList= false;          // Use of R3B special physics list
  Bool_t fR3BMagnet = false;          // Magnetic field definition
  Bool_t fCalifaHitFinder = false;    // Apply hit finder task

  TString fMC = "TGeant4";           // MonteCarlo engine: TGeant3, TGeant4, TFluka
  TString fGenerator = "ascii";        // Event generator type: box, gammas, r3b, ion, ascii
  TString fEventFile = "evt_gen3.dat";           // Input event file in the case of ascii generator
//  TString fGenerator = "box";        // Event generator type: box, gammas, r3b, ion, ascii
//  TString fEventFile = "";           // Input event file in the case of ascii generator
  Int_t randomSeed = 0; // 0 for time-dependent random numbers, 335566 for time-independent


  // ---------------  Detector selection: true - false -------------------------------

  Bool_t  fTarget = false;           // Target
  TString fTargetType = "LiH";       // Target selection: LeadTarget, Para, Para45, LiH


  Bool_t  fTracker = true;          // Tracker
  //TString fTrackerGeo = "ams_s444.geo.root"; //geometry for kinematic case
  TString fTrackerGeo = "ams_s444_real.geo.root"; //geometry for realistic case
  cout << "Geometry " << fTrackerGeo << " is chosen!" << endl;


  // ========= End of Configuration area =======================


  TString dir = gSystem->Getenv("VMCWORKDIR");
  TString r3bdir = dir + "/macros/";
  r3bdir.ReplaceAll("//","/");

  TString r3b_geomdir = dir + "/geometry/";
  gSystem->Setenv("GEOMPATH",r3b_geomdir.Data());
  r3b_geomdir.ReplaceAll("//","/");

  TString r3b_confdir = dir + "/gconfig/";
  gSystem->Setenv("CONFIG_DIR",r3b_confdir.Data());
  r3b_confdir.ReplaceAll("//","/");


  // ----    Debug option   -------------------------------------------------
  gDebug = 0;
  // ------------------------------------------------------------------------

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------


  // -----   Create simulation run   ----------------------------------------
  FairRunSim* run = new FairRunSim();
  run->SetName(fMC);                  // Transport engine
  run->SetOutputFile(OutFile);        // Output file

  //  R3B Special Physics List in G4 case
  if ( (fUserPList) && (fMC.CompareTo("TGeant4") == 0) ) {
       run->SetUserConfig("g4R3bConfig.C");
       run->SetUserCuts("SetCuts.C");
   }

  // -----   Create media   -------------------------------------------------
  run->SetMaterials("media_r3b.geo");       // Materials

  // -----   Create R3B geometry --------------------------------------------

  //Cave definition
  FairModule* cave= new R3BCave("CAVE");
  cave->SetGeometryFileName("r3b_cave.geo");
 // cave->SetGeometryFileName("r3b_cave_vacuum.geo");
  run->AddModule(cave);

  //Target definition
  if (fTarget) {
    run->AddModule(new R3BTarget(fTargetType, "target_" + fTargetType + ".geo.root"));
  }



  // Tracker
  if (fTracker) {
    run->AddModule(new R3BTra(fTrackerGeo));
  }


  // -----   Create PrimaryGenerator   --------------------------------------

  // 1 - Create the Main API class for the Generator
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();

  if (fGenerator.CompareTo("box") == 0  ) {
	  // 2- Define the BOX generator
	/*  Double_t pdgId=211; // pion beam
	  Double_t theta1= 30.;  // polar angle distribution
	  Double_t theta2= 160.;
	  Double_t momentum=.8; // 10 GeV/c
	  FairBoxGenerator* boxGen = new FairBoxGenerator(pdgId, 10);
	  boxGen->SetThetaRange (   theta1,   theta2);
	  boxGen->SetPRange     (momentum,momentum*2.);
	  boxGen->SetPhiRange   (0.,360.);
	  boxGen->SetXYZ(0.0,0.0,0.0);
  */
    //FairIonGenerator* boxGen = new FairIonGenerator(Z, A, Q, multiplicity, Px, Py, Pz, vertex_x, vertex_y, vertex_z);
    FairIonGenerator* boxGen = new FairIonGenerator(17, 30, 17, 1, 0., 0., 1.4, 0., 0., 0.);
   // primGen->AddGenerator(boxGen);
	  // add the box generator
    boxGen->SetMass(27.950444);
	  primGen->AddGenerator(boxGen);
  }


  if (fGenerator.CompareTo("ascii") == 0  ) {
    R3BAsciiGenerator* gen = new R3BAsciiGenerator((dir+"/input/"+fEventFile).Data());
    primGen->AddGenerator(gen);
  }

  run->SetGenerator(primGen);


  //-------Set visualisation flag to true------------------------------------
  run->SetStoreTraj(fVis);

  FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");


  // -----   Initialize simulation run   ------------------------------------

  TRandom3 random(randomSeed);
  gRandom = &random;
  //TVirtualMC::GetMC()->SetRandom(new TRandom3(randomSeed));
  run->Init();

  // -----   Runtime database   ---------------------------------------------
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  Bool_t kParameterMerged = kTRUE;
  FairParRootFileIo* parOut = new FairParRootFileIo(kParameterMerged);
  parOut->open(ParFile.Data());
  rtdb->setOutput(parOut);
  rtdb->saveOutput();
  //rtdb->print();

  // -----   Start run   ----------------------------------------------------
  if (nEvents>0) run->Run(nEvents);

  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is "    << OutFile << endl;
  cout << "Parameter file is " << ParFile << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime
       << "s" << endl << endl;
  // ------------------------------------------------------------------------

  cout << " Test passed" << endl;
  cout << " All ok " << endl;

}


