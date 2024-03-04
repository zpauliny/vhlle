/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include <queue>
#include "fld.h"
#include "hdo.h"
#include "ic.h"
#include "ickw.h"
#include "icPartUrqmd.h"
#include "icPartSMASH.h"
#include "icGlauber.h"
#include "icGubser.h"
#include "icGlissando.h"
#include "icTrento.h"
#include "icTest.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoAZH.h"
#include "eoHadron.h"
#include "eoSmash.h"
#include "trancoeff.h"
#include "particle.h"

using namespace std;

// program parameters, to be read from file
int nx, ny, nz, eosType, etaSparam = 0, zetaSparam = 0;
int eosTypeHadron = 0;
// in Cartesian frame etamin, etamax = zmin, zmax
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, tauResize, dtau;
string collSystem, outputDir, isInputFile;
double etaS, zetaS, eCrit, eEtaSMin, al, ah, aRho, T0, etaSMin;
int icModel,glauberVariable =1;  // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, Rgt, Rgz, impactPar, s0ScaleFactor;
bool freezeoutOnly {false};  // freezoutOnly 1 for true, 0 for false

void setDefaultParameters() {
  // specifically for dynamical initialization, do not resize
  tauResize = 100.0;
}

void readParameters(char *parFile) {
 char parName[255], parValue[255];
 ifstream fin(parFile);
 if (!fin.is_open()) {
  cout << "cannot open parameters file " << parFile << endl;
  exit(1);
 }
 cout << "vhlle: reading parameters from " << parFile << endl;
 while (fin.good()) {
  string line;
  getline(fin, line);
  istringstream sline(line);
  sline >> parName >> parValue;
  if (strcmp(parName, "eosType") == 0)
   eosType = atoi(parValue);
  else if (strcmp(parName, "eosTypeHadron") == 0)
   eosTypeHadron = atoi(parValue);
  else if (strcmp(parName, "nx") == 0)
   nx = atoi(parValue);
  else if (strcmp(parName, "ny") == 0)
   ny = atoi(parValue);
  else if (strcmp(parName, "nz") == 0)
   nz = atoi(parValue);
  else if (strcmp(parName, "icModel") == 0)
   icModel = atoi(parValue);
  else if (strcmp(parName, "glauberVar") == 0)
   glauberVariable = atoi(parValue);
  else if (strcmp(parName, "xmin") == 0)
   xmin = atof(parValue);
  else if (strcmp(parName, "xmax") == 0)
   xmax = atof(parValue);
  else if (strcmp(parName, "ymin") == 0)
   ymin = atof(parValue);
  else if (strcmp(parName, "ymax") == 0)
   ymax = atof(parValue);
  else if (strcmp(parName, "etamin") == 0)
   etamin = atof(parValue);
  else if (strcmp(parName, "etamax") == 0)
   etamax = atof(parValue);
  else if (strcmp(parName, "tau0") == 0)
   tau0 = atof(parValue);
  else if (strcmp(parName, "tauMax") == 0)
   tauMax = atof(parValue);
  else if (strcmp(parName, "tauGridResize") == 0)
   tauResize = atof(parValue);
  else if (strcmp(parName, "dtau") == 0)
   dtau = atof(parValue);
  else if (strcmp(parName, "e_crit") == 0)
   eCrit = atof(parValue);
  else if (strcmp(parName, "etaS") == 0)
   etaS = atof(parValue);
  else if (strcmp(parName, "zetaS") == 0)
   zetaS = atof(parValue);
  else if (strcmp(parName, "zetaSparam") == 0)
   zetaSparam = atoi(parValue);
  else if (strcmp(parName, "epsilon0") == 0)
   epsilon0 = atof(parValue);
  else if (strcmp(parName, "Rg") == 0)
   Rgt = atof(parValue);
  else if (strcmp(parName, "Rgz") == 0)
   Rgz = atof(parValue);
  else if (strcmp(parName, "impactPar") == 0)
   impactPar = atof(parValue);
  else if (strcmp(parName, "s0ScaleFactor") == 0)
   s0ScaleFactor = atof(parValue);
  else if (strcmp(parName, "etaSparam") == 0)
   etaSparam = atoi(parValue);
  else if (strcmp(parName, "aRho") == 0)
   aRho = atof(parValue);
  else if (strcmp(parName, "ah") == 0)
   ah = atof(parValue);
  else if (strcmp(parName, "al") == 0)
   al = atof(parValue);
  else if (strcmp(parName, "T0") == 0)
   T0 = atof(parValue);
  else if (strcmp(parName, "eEtaSMin") == 0)
   eEtaSMin = atof(parValue);
  else if (strcmp(parName, "etaSMin") == 0)
   etaSMin = atof(parValue);
  else if (strcmp(parName, "freezeoutOnly") == 0)
   freezeoutOnly = atoi(parValue);
  else if (parName[0] == '!')
   cout << "CCC " << sline.str() << endl;
  else
   cout << "UUU " << sline.str() << endl;
 }
}

void printParameters() {
 cout << "====== parameters ======\n";
 cout << "outputDir = " << outputDir << endl;
 cout << "freezeoutOnly = " << freezeoutOnly << endl;
 cout << "eosType = " << eosType << endl;
 cout << "eosTypeHadron = " << eosTypeHadron << endl;
 cout << "nx = " << nx << endl;
 cout << "ny = " << ny << endl;
 cout << "nz = " << nz << endl;
 cout << "icModel = " << icModel << endl;
 cout << "glauberVar = " << glauberVariable << "   ! 0=epsilon,1=entropy"
      << endl;
 cout << "xmin = " << xmin << endl;
 cout << "xmax = " << xmax << endl;
 cout << "ymin = " << ymin << endl;
 cout << "ymax = " << ymax << endl;
 cout << "etamin = " << etamin << endl;
 cout << "etamax = " << etamax << endl;
 cout << "tau0 = " << tau0 << endl;
 cout << "tauMax = " << tauMax << endl;
 cout << "tauGridResize = " << tauResize << endl;
 cout << "dtau = " << dtau << endl;
 cout << "e_crit = " << eCrit << endl;
 cout << "zeta/s param : " << zetaSparam << endl;
 cout << "etaSparam = " << etaSparam << endl;
 if (etaSparam == 0){
    cout << "eta/s = " << etaS << endl;
 }
 else if (etaSparam == 1){
    cout << "al = " << al << endl;
    cout << "ah = " << ah << endl;
    cout << "etaSMin = " << etaSMin << endl;
    cout << "T0 = " << T0 << endl;
 }
 else if (etaSparam == 2){
    cout << "al = " << al << endl;
    cout << "ah = " << ah << endl;
    cout << "aRho = " << aRho << endl;
    cout << "etaSMin = " << etaSMin << endl;
    cout << "eEtaSMin = " << eEtaSMin << endl;
 }
 cout << "zeta/s = " << zetaS << endl;
 cout << "epsilon0 = " << epsilon0 << endl;
 cout << "Rgt = " << Rgt << "  Rgz = " << Rgz << endl;
 cout << "impactPar = " << impactPar << endl;
 cout << "s0ScaleFactor = " << s0ScaleFactor << endl;
 cout << "======= end parameters =======\n";
}

void readCommandLine(int argc, char** argv)
{
  if(argc==1){
  cout << "no CL params - exiting.\n"; exit(1) ;
 }
 else{
  for(int iarg=1; iarg<argc-1; iarg++){
   if(strcmp(argv[iarg],"-system")==0) collSystem = argv[iarg+1];
   if(strcmp(argv[iarg],"-params")==0) readParameters(argv[iarg+1]);
   if(strcmp(argv[iarg],"-ISinput")==0) isInputFile = argv[iarg+1];
   if(strcmp(argv[iarg],"-outputDir")==0) outputDir = argv[iarg+1];
  }
  cout << "vhlle: command line parameters are:\n";
  cout << "collision system:  " << collSystem << endl;
  cout << "ini.state input:  " << isInputFile << endl;
  cout << "output directory:  " << outputDir << endl;
 }
}


Fluid* expandGrid2x(Hydro* h, EoS* eos, EoS* eosH, TransportCoeff *trcoeff) {
 Fluid* f = h->getFluid();
 if(f->getX(0) + f->getX(f->getNX()-1)>0.001
  || f->getY(0) + f->getY(f->getNY()-1)>0.001) {
   cout << "this grid expansion works only with symmetric min/max ranges\n";
   return f;
 }
 // creating a new fluid, twice the transverse size of the current one
 Fluid* fnew = new Fluid(eos, eosH, trcoeff, f->getNX(), f->getNY(), f->getNZ(),
   2.0*f->getX(0), 2.0*f->getX(f->getNX()-1), 2.0*f->getY(0), 2.0*f->getY(f->getNY()-1),
   f->getZ(0), f->getZ(f->getNZ()-1), 2.0*h->getDtau(), f->geteCrit());
 // filling the new fluid
 for(int ix=0; ix<f->getNX(); ix++)
  for(int iy=0; iy<f->getNY(); iy++)
   for(int iz=0; iz<f->getNZ(); iz++) {
    fnew->getCell(ix, iy, iz)->importVars(f->getCell(2*(ix - f->getNX()/2) + f->getNX()/2,
      2*(iy - f->getNY()/2) + f->getNY()/2, iz));
   }
 h->setFluid(fnew);  // now Hydro object operates on the new fluid
 h->setDtau(2.0*h->getDtau());
 delete f;
 return fnew;
}

void output_e_nb(double t, Fluid* f, ofstream& file_e, ofstream& file_nb) {
 double x, y, z;
 double e, p, nb, nq, ns, vx, vy, vz;
 Cell *c;
 file_e << t << "\n";
 file_nb << t << "\n";
 for (int ix = 0; ix < f->getNX(); ix++) 
   for (int iy = 0; iy < f->getNY(); iy++) 
    for (int iz = 0; iz < f->getNZ(); iz++) {
     c = f->getCell(ix, iy, iz);
     x = f->getX(ix);
     y = f->getY(iy);
     z = f->getZ(iz);
     f->getCMFvariables(c, 1.0, e, nb, nq, ns, vx, vy, vz);
     file_e << e << " ";
     file_nb << nb << " "; 
  }
  file_e << "\n";
  file_nb << "\n";
}

// program parameters, to be read from file
// int nx, ny, nz, eosType ;
// double xmin, xmax, ymin, ymax, zmin, zmax, tau0, tauMax, dtau ;
// char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
// char icInputFile [255] ;
// double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume ;
// int icModel, NPART, glauberVariable=1 ;
// double epsilon0, alpha, impactPar, s0ScaleFactor ;

int main(int argc, char **argv) {
 // pointers to all the main objects
 EoS *eos;
 EoS *eosH;
 TransportCoeff *trcoeff;
 Fluid *f;
 Hydro *h;
 queue<Particle>* particles = new queue<Particle>();
 time_t start = 0, end;
 time(&start);
 Particle one_particle;

 // read parameters from file
 setDefaultParameters();
 readCommandLine(argc, argv);
 printParameters();

 // EoS for hydro evolution
 if (eosType == 0)
  eos = new EoSs("eos/Laine_nf3.dat", 3);
 else if (eosType == 1)
  eos = new EoSChiral();
 else if (eosType == 2)
  eos = new EoSAZH();
 else {
  cout << "eosType != 0,1,2\n";
  return 0;
 }

 // hadronic EoS for hypersurface creation
 if (eosTypeHadron == 0) {
   eosH = new EoSHadron((char*)"eos/eosHadronLog.dat"); //PDG hadronic EoS
 } else if (eosTypeHadron == 1) {
   eosH = new EoSSmash((char*)"eos/hadgas_eos_SMASH.dat", 101, 51, 51); //SMASH hadronic EoS
 } else {
   cout << "Unknown haronic EoS type for hypersurface creation.\n" <<
           "eosTypeHadron should be either \"0\" (PDG hadronic EoS) or " <<
           "\"1\" (SMASH hadronic EoS).\n";
   return 0;
 }


 // transport coefficients
 trcoeff = new TransportCoeff(etaS, zetaS, zetaSparam, eos, etaSparam, ah, al, aRho, T0, etaSMin, eEtaSMin);

 f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit);
 cout << "fluid allocation done\n";

 double* timeInit = new double; // current time, tau or t depending on the coordinate frame

 // initial conditions
 if (icModel == 1) {  // optical Glauber
  ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 2) {  // Glauber_table + parametrized rapidity dependence
  IC *ic = new IC(isInputFile.c_str(), s0ScaleFactor);
  ic->setIC(f, eos, tau0);
  delete ic;
 } else if (icModel == 3) {  // UrQMD IC
  IcPartUrqmd *ic = new IcPartUrqmd(f, isInputFile.c_str(), Rgt, Rgz, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 4) {  // analytical Gubser solution
  ICGubser *ic = new ICGubser();
  ic->setIC(f, eos, tau0);
  delete ic;
  }else if(icModel==5){ // IC from GLISSANDO + rapidity dependence
   IcGlissando *ic = new IcGlissando(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else if (icModel == 6){ // SMASH IC
   IcPartSMASH *ic = new IcPartSMASH(f, isInputFile.c_str(), Rgt, Rgz, tau0);
   ic->setIC(f, eos);
   delete ic;
 } else if(icModel==7){ // IC from Trento
   IcTrento *ic = new IcTrento(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else if(icModel==8){ // IC for testing purposes
   ICTest *ic = new ICTest();
   ic->setIC(f, eos, 1);
   delete ic;
 } else if(icModel==9) { // SMASH dynamical IC
   IcPartSMASH *ic = new IcPartSMASH(f, isInputFile.c_str(), Rgt, Rgz, particles);
   ic->setIC(f,eos,particles,timeInit);
   delete ic;
 } else {
  cout << "icModel = " << icModel << " not implemented\n";
 }
 cout << "IC done\n";

 // For calculating initial anisotropy without running full hydro, uncomment following line
 //f->InitialAnisotropies(tau0) ;

 time_t tinit = 0;
 time(&tinit);
 float diff = difftime(tinit, start);
 cout << "Init time = " << diff << " [sec]" << endl;

 // hydro init
 double ctime;
 #ifdef CARTESIAN
 h = new Hydro(f, eos, trcoeff, *timeInit, dtau);
 ctime = h->time();
 cout << "Hydro starting at " << ctime << endl;
 double ftime = ctime + 2.0;
 #else
 h = new Hydro(f, eos, trcoeff, tau0, dtau);
 #endif
 time(&start);
 // h->setNSvalues() ; // initialize viscous terms
 f->initOutput(outputDir.c_str(), tau0, freezeoutOnly);
 //f->outputCorona(tau0);

 // initialize energy density output
 string outfile_e = outputDir.c_str();
 outfile_e.append("/energy_density.dat");
 ofstream file_e(outfile_e.c_str());
 file_e << "# block at t, energy_density: iterating over z, y, x \n";
 file_e << "# grid: xmin, xmax, ymin, ymax, zmin, zmax: " << "\n";
 file_e << "# " << xmin << " " << xmax << " " <<
   ymin << " " << ymax << " " << etamin << " " <<
   etamax;
 file_e << "# Nx, Ny, Nz: \n";
 file_e << "# " << nx << " " << ny << " " << nz << "\n";

 // initialize baryon number density output
 string outfile_nb = outputDir.c_str();
 outfile_nb.append("/baryon_number.dat");
 ofstream file_nb(outfile_nb.c_str());
 file_nb << "# block at t, energy_density: iterating over z, y, x \n";
 file_nb << "# grid: xmin, xmax, ymin, ymax, zmin, zmax: " << "\n";
 file_nb << "# " << xmin << " " << xmax << " " <<
   ymin << " " << ymax << " " << etamin << " " <<
   etamax;
 file_nb << "# Nx, Ny, Nz: \n";
 file_nb << "# " << nx << " " << ny << " " << nz << "\n";


 bool resized = false; // flag if the grid has been resized
 
 int timestep = 0;
 int nelements = 1;  
 do {
  // small tau: decrease timestep by making substeps, in order
  // to avoid instabilities in eta direction (signal velocity ~1/tau)
  int nSubSteps = 1;
  #ifdef CARTESIAN
  ctime = h->time();
  #else
  ctime = h->getTau();
  #endif
  while (dtau / nSubSteps >
         1.0 * ctime * (etamax - etamin) / (nz - 1)) {
   nSubSteps *= 2;  // 0.02 in "old" coordinates
  }
  if(nSubSteps>1) {
   h->setDtau(h->getDtau() / nSubSteps);
   for (int j = 0; j < nSubSteps; j++){
    h->performStep();
   }
   h->setDtau(h->getDtau() * nSubSteps);
   cout << "timestep reduced by " << nSubSteps << endl;
  } else 
   h->performStep();
   
  if (particles->size() > 0) h->addParticles(particles);
  
  // freeze-out only after ftime and up until nelemens is 0
  if ((ctime > ftime) && (nelements>0))
   nelements = f->outputSurface(ctime);
  if (!freezeoutOnly)
   f->outputGnuplot(ctime);

  // when nelements is 0: print out the particles still in the queue
  // file an smash oscar format
  if (nelements == 0) {
    int particle_number = particles->size();
    int n_part = 0;
    if (particle_number > 0) {
      string filename = outputDir.c_str();
      filename.append("/particle_lists.oscar");
      ofstream outfile(filename.c_str());
      outfile << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge \n";
      outfile << "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e \n";
      while (particle_number > 0)
      {
        Particle particle_to_dump = particles->front();
        double t = particle_to_dump.getT();
        double x = particle_to_dump.getX();
        double y = particle_to_dump.getY();
        double z = particle_to_dump.getZ();
        double mass = particle_to_dump.getM();
        double p0 = particle_to_dump.getE();
        double px = particle_to_dump.getPx();
        double py = particle_to_dump.getPy();
        double pz = particle_to_dump.getPz();
        int pdg = particle_to_dump.getPdg();
        int id = 0;
        int charge = particle_to_dump.getQ();
        outfile << t << " " << x << " " << y << " " << z << " " << mass << " " << 
                   p0 << " " << px << " " << py << " " << pz << " " << pdg << 
                   " " << id << " " << charge << "\n";
        n_part += 1;
        particles->pop();
        particle_number = particles->size();
      }
    }
    nelements = -1;
    // check-in about the procedure
    std::cout << n_part << " particles were left after hydro evolution, " <<
                 "nelements = " << nelements << ", particle queue is empty: " <<
                 particles->empty() << "\n";
  }
  
  // output energy density at every 10th timestep
  if (timestep%10 == 0) output_e_nb(ctime, f, file_e, file_nb); 
    
  if(ctime>=tauResize and resized==false) {
   cout << "grid resize\n";
   f = expandGrid2x(h, eos, eosH, trcoeff);
   resized = true;
  }
  timestep++;
 } while(ctime<tauMax+0.0001);

 end = 0;
 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;
 
 cout << "Is particle queue empty: " << particles->empty() << endl;

 delete f;
 delete h;
 delete eos;
 delete eosH;
 delete particles;
 file_e.close();
 file_nb.close();
}
