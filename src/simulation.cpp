#include <cstring>
#include <ctime>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include <vhlle/colour.h>
#include <vhlle/eoAZH.h>
#include <vhlle/eoChiral.h>
#include <vhlle/eoCMF.h>
#include <vhlle/eoCMFe.h>
#include <vhlle/eoHadron.h>
#include <vhlle/eos.h>
#include <vhlle/eoSmash.h>
#include <vhlle/eo1.h>
#include <vhlle/eo3.h>
#include <vhlle/fld.h>
#include <vhlle/hdo.h>
#include <vhlle/ic.h>
#include <vhlle/icDynFlu.h>
#include <vhlle/icGlauber.h>
#include <vhlle/icGlissando.h>
#include <vhlle/icGubser.h>
#include <vhlle/ickw.h>
#include <vhlle/icPartSMASH.h>
#include <vhlle/icPartUrqmd.h>
#include <vhlle/icTest.h>
#include <vhlle/icTrento.h>
#include <vhlle/icTrento3d.h>
#include <vhlle/icSuperMC.h>
#include <vhlle/particle.h>
#include <vhlle/simulation.h>
#include <vhlle/trancoeff.h>
#include <vhlle/vtk.h>

// ============================================================
// constructor / destructor
// ============================================================

Simulation::Simulation(int argc, char** argv) {
    particles = new std::deque<Particle>();
    time(&tstart);
    readCommandLine(argc, argv);
}

Simulation::~Simulation() {
    delete f;
    delete h;
    delete eos;
    delete eosH;
    delete trcoeff;
    delete vtk_out;
    delete particles;
}

// ============================================================
// 2. setup
// ============================================================

// setup eos, setup IC, setup hydro
void Simulation::setup() {
    // EoS for hydro evolution
    if (eosType == 0)
        eos = new EoSs("eos/Laine_nf3.dat", 3);
    else if (eosType == 1)
        eos = new EoSChiral();
    else if (eosType == 2)
        eos = new EoSAZH();
    else if (eosType == 3)
        eos = new EoSCMF();
    else if (eosType == 4)
        eos = new EoSCMFe();
    else {
        std::cerr << red << "FATAL: eosType != 0,1,2,3,4\n" << reset;
        exit(1);
    }

    // hadronic EoS for hypersurface creation
    if (eosTypeHadron == 0)
        eosH = new EoSHadron((char*)"eos/eosHadronLog.dat");
    else if (eosTypeHadron == 1)
        eosH = new EoSSmash((char*)"eos/hadgas_eos_SMASH.dat", 101, 51, 51);
    else {
        std::cerr << red << "FATAL: Unknown hadronic EoS type.\n" << reset;
        exit(1);
    }

    // transport coefficients
    trcoeff = new TransportCoeff(etaS, zetaS, ah, al, aRho, T0, etaSMin,
        etaSEpsilonMin, etaSShiftMuB, etaSScaleMuB, zetaSPeakEpsilon,
        zetaSScaleBeta, zetaSSigmaMinus, zetaSSigmaPlus, eos, etaSparam,
        zetaSparam);

    // fluid grid
    f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax,
                  etamin, etamax, dtau, eCrit, cartesian);
    std::cout << "fluid allocation done\n";

    // initial conditions
    if (icModel == 1) {
        ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0);
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 2) {
        IC *ic = new IC(isInputFile.c_str(), s0ScaleFactor, glauberVariable);
        ic->setIC(f, eos, tau0);
        delete ic;
    } else if (icModel == 3) {
        IcPartUrqmd *ic = new IcPartUrqmd(f, isInputFile.c_str(), Rgt, Rgz, tau0);
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 4) {
        ICGubser *ic = new ICGubser();
        ic->setIC(f, eos, tau0);
        delete ic;
    } else if (icModel == 5) {
        IcGlissando *ic = new IcGlissando(f, isInputFile.c_str(), tau0, collSystem.c_str());
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 6) {
        IcPartSMASH *ic = new IcPartSMASH(f, isInputFile.c_str(), Rgt, Rgz, smoothingType);
        tau0 = ic->getTau0();
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 7) {
        IcTrento *ic = new IcTrento(f, isInputFile.c_str(), tau0, collSystem.c_str());
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 8) {
        IcTrento3d *ic = new IcTrento3d(f, isInputFile.c_str(), tau0, collSystem.c_str());
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 9) {
        ICSuperMC *ic = new ICSuperMC(isInputFile, tau0, collSystem);
        ic->setIC(f, eos);
        delete ic;
    } else if (icModel == 10) {
        if (!cartesian) {
            std::cerr << red << "FATAL: IC model = 10: Dynamical Fluidization\n"
                 << "Cartesian coordinate system must be set.\n" << reset;
            exit(1);
        }
        IcDynFlu *ic = new IcDynFlu(f, isInputFile.c_str(), gaussian_sigma, particles);
        ic->setIC(f, eos, particles, ctime, minParticlesFO, timeInitFO);
        delete ic;
    } else if (icModel == 11) {
        ICTest *ic = new ICTest();
        ic->setIC(f, eos, 1);
        delete ic;
    } else {
        std::cerr << red << "FATAL: icModel = " << icModel << " not implemented.\n"
                  << reset;
        exit(1);
    }
    std::cout << "IC done\n";

    time(&tinit);
    std::cout << "Init time = " << difftime(tinit, tstart) << " [sec]\n";

    // hydro init
    if (cartesian) {
        h = new Hydro(f, eos, trcoeff, ctime, dtau, cartesian);
        ctime = h->getTime();
    } else {
        h = new Hydro(f, eos, trcoeff, tau0, dtau, cartesian);
        ctime = h->getTau();
    }

    if (vorticityOn)
        h->enableVorticity();

    time(&tstart);  // reset wall clock for hydro loop timing

    f->initOutput(outputDir.c_str(), tau0, freezeoutOnly);
    if (icModel != 10) {
        f->outputCorona(tau0, freezeoutExtend);
        if (vorticityOn) f->printDbetaHeader();
    }

    vtk_out = new VtkOutput(outputDir, eos, xmin, ymin, etamin, cartesian);
}

// ============================================================
// 3. step
// ============================================================

bool Simulation::step() {
    if (!vtk_values.empty())
        vtk_out->write(*h, vtk_values);

    // substep logic to avoid instabilities at small tau
    int nSubSteps = 1;
    while (dtau / nSubSteps > 1.0 * ctime * (etamax - etamin) / (nz - 1))
        nSubSteps *= 2;

    if (nSubSteps > 1) {
        h->setDtau(h->getDtau() / nSubSteps);
        for (int j = 0; j < nSubSteps; j++)
            h->performStep();
        h->setDtau(h->getDtau() * nSubSteps);
        std::cout << "timestep reduced by " << nSubSteps << "\n";
    } else {
        h->performStep();
    }

    ctime = cartesian ? h->getTime() : h->getTau();

    if (icModel == 10) {
        if (particles->size() > 0) h->addParticles(particles);
        if ((ctime > timeInitFO) && (nelements > 0))
            nelements = f->outputSurface(ctime, freezeoutExtend);
    } else {
        nelements = f->outputSurface(ctime, freezeoutExtend);
    }

    if (!freezeoutOnly)
        f->outputGnuplot(ctime);

    if (nelements == 0) {
        if (icModel == 10)
            outputCoronaParticles(particles, outputDir);
        return false;
    }

    if (ctime >= tauResize && !resized) {
        std::cout << "grid resize\n";
        expandGrid2x();
        resized = true;
    }

    return ctime < tauMax + 0.0001;
}

// ============================================================
// 4. finalize
// ============================================================

void Simulation::finalize() {
    time_t tend = 0;
    time(&tend);
    std::cout << "Execution time = " << difftime(tend, tstart) << " [sec]\n";
    f->renameOutput(outputDir.c_str());
}

// ============================================================
// private helpers
// ============================================================

void Simulation::expandGrid2x() {
    if (f->getX(0) + f->getX(f->getNX()-1) > 0.001
     || f->getY(0) + f->getY(f->getNY()-1) > 0.001) {
        std::cout << "grid expansion works only with symmetric min/max ranges\n";
        return;
    }
    Fluid* fnew = new Fluid(eos, eosH, trcoeff,
        f->getNX(), f->getNY(), f->getNZ(),
        2.0*f->getX(0), 2.0*f->getX(f->getNX()-1),
        2.0*f->getY(0), 2.0*f->getY(f->getNY()-1),
        f->getZ(0), f->getZ(f->getNZ()-1),
        2.0*h->getDtau(), f->geteCrit(), cartesian);
    if (vorticityOn) fnew->enableVorticity();
    for (int ix = 0; ix < f->getNX(); ix++)
        for (int iy = 0; iy < f->getNY(); iy++)
            for (int iz = 0; iz < f->getNZ(); iz++) {
                fnew->getCell(ix, iy, iz)->importVars(
                    f->getCell(2*(ix - f->getNX()/2) + f->getNX()/2,
                               2*(iy - f->getNY()/2) + f->getNY()/2, iz));
                if (vorticityOn) fnew->getCell(ix, iy, iz)->enableVorticity();
            }
    h->setFluid(fnew);
    h->setDtau(2.0 * h->getDtau());
    delete f;
    f = fnew;
}

void Simulation::checkGridDimension(int n, char axis) {
    if (n < 5) {
        std::cerr << red << "FATAL: grid too small in " << axis << " direction\n" << reset;
        exit(1);
    }
}

void Simulation::checkGridBorders(double min, double max, const std::string& axis) {
    if (min >= max) {
        std::cerr << red << "FATAL: " << axis << "min >= " << axis << "max\n" << reset;
        exit(1);
    }
}

bool Simulation::parse_bool(const std::string& value) {
    if (value == "1" || value == "true")  return true;
    if (value == "0" || value == "false") return false;
    throw std::runtime_error("Invalid boolean in config: " + value);
}

void Simulation::readCommandLine(int argc, char** argv) {
    if (argc == 1) {
        std::cerr << "FATAL: no CL params.\n" << reset;
        exit(1);
    }
    for (int i = 1; i < argc - 1; i++) {
        if (strcmp(argv[i], "-system")    == 0) collSystem  = argv[i+1];
        if (strcmp(argv[i], "-params")    == 0) readParameters(argv[i+1]);
        if (strcmp(argv[i], "-ISinput")   == 0) isInputFile = argv[i+1];
        if (strcmp(argv[i], "-outputDir") == 0) outputDir   = argv[i+1];
    }
    std::cout << "collision system: " << collSystem << "\n"
              << "ini.state input:  " << isInputFile << "\n"
              << "output directory: " << outputDir << "\n";
}

void Simulation::readParameters(const char* parFile) {
    std::ifstream fin(parFile);
    if (!fin.is_open()) {
        std::cerr << red << "FATAL: cannot open parameters file " << parFile << "\n"
                  << reset;
        exit(1);
    }
    std::cout << "vhlle: reading parameters from " << parFile << "\n";

    std::map<std::string, std::function<void(const std::string&)>> handlers = {
        {"eosType",          [this](const std::string& v) { eosType          = stoi(v); }},
        {"eosTypeHadron",    [this](const std::string& v) { eosTypeHadron    = stoi(v); }},
        {"nx",               [this](const std::string& v) { nx               = stoi(v); }},
        {"ny",               [this](const std::string& v) { ny               = stoi(v); }},
        {"nz",               [this](const std::string& v) { nz               = stoi(v); }},
        {"icModel",          [this](const std::string& v) { icModel          = stoi(v); }},
        {"glauberVar",       [this](const std::string& v) { glauberVariable  = stoi(v); }},
        {"xmin",             [this](const std::string& v) { xmin             = stod(v); }},
        {"xmax",             [this](const std::string& v) { xmax             = stod(v); }},
        {"ymin",             [this](const std::string& v) { ymin             = stod(v); }},
        {"ymax",             [this](const std::string& v) { ymax             = stod(v); }},
        {"etamin",           [this](const std::string& v) { etamin           = stod(v); }},
        {"etamax",           [this](const std::string& v) { etamax           = stod(v); }},
        {"tau0",             [this](const std::string& v) { tau0             = stod(v); }},
        {"tauMax",           [this](const std::string& v) { tauMax           = stod(v); }},
        {"tauGridResize",    [this](const std::string& v) { tauResize        = stod(v); }},
        {"dtau",             [this](const std::string& v) { dtau             = stod(v); }},
        {"e_crit",           [this](const std::string& v) { eCrit            = stod(v); }},
        {"etaS",             [this](const std::string& v) { etaS             = stod(v); }},
        {"zetaS",            [this](const std::string& v) { zetaS            = stod(v); }},
        {"etaSparam",        [this](const std::string& v) { etaSparam        = stoi(v); }},
        {"zetaSparam",       [this](const std::string& v) { zetaSparam       = stoi(v); }},
        {"zetaSScaleBeta",   [this](const std::string& v) { zetaSScaleBeta   = stod(v); }},
        {"zetaSPeakEpsilon", [this](const std::string& v) { zetaSPeakEpsilon = stod(v); }},
        {"zetaSSigmaMinus",  [this](const std::string& v) { zetaSSigmaMinus  = stod(v); }},
        {"zetaSSigmaPlus",   [this](const std::string& v) { zetaSSigmaPlus   = stod(v); }},
        {"epsilon0",         [this](const std::string& v) { epsilon0         = stod(v); }},
        {"Rg",               [this](const std::string& v) { Rgt              = stod(v); }},
        {"Rgz",              [this](const std::string& v) { Rgz              = stod(v); }},
        {"impactPar",        [this](const std::string& v) { impactPar        = stod(v); }},
        {"s0ScaleFactor",    [this](const std::string& v) { s0ScaleFactor    = stod(v); }},
        {"VTK_output_values",[this](const std::string& v) { vtk_values       = v;       }},
        {"aRho",             [this](const std::string& v) { aRho             = stod(v); }},
        {"ah",               [this](const std::string& v) { ah               = stod(v); }},
        {"al",               [this](const std::string& v) { al               = stod(v); }},
        {"T0",               [this](const std::string& v) { T0               = stod(v); }},
        {"etaSEpsilonMin",   [this](const std::string& v) { etaSEpsilonMin   = stod(v); }},
        {"etaSMin",          [this](const std::string& v) { etaSMin          = stod(v); }},
        {"etaSShiftMuB",     [this](const std::string& v) { etaSShiftMuB     = stod(v); }},
        {"etaSScaleMuB",     [this](const std::string& v) { etaSScaleMuB     = stod(v); }},
        {"freezeoutOnly",    [this](const std::string& v) { freezeoutOnly    = parse_bool(v); }},
        {"freezeoutExtend",  [this](const std::string& v) { freezeoutExtend  = parse_bool(v); }},
        {"vorticity",        [this](const std::string& v) { vorticityOn      = stoi(v); }},
        {"smoothingType",    [this](const std::string& v) { smoothingType    = stoi(v); }},
        {"Gaussian_Sigma",   [this](const std::string& v) { gaussian_sigma   = stod(v); }},
        {"minParticlesFO",   [this](const std::string& v) { minParticlesFO   = stoi(v); }},
        {"cartesian",        [this](const std::string& v) { cartesian        = parse_bool(v); }},
    };

    char parName[255], parValue[255];
    while (fin.good()) {
        std::string line;
        getline(fin, line);
        std::istringstream sline(line);
        sline >> parName >> parValue;
        auto handler = handlers.find(parName);
        if (handler != handlers.end())
            handler->second(parValue);
        else if (parName[0] == '!')
            std::cout << "CCC " << sline.str() << "\n";
        else
            std::cout << "UUU " << sline.str() << "\n";
    }

    checkGridBorders(xmin, xmax, "x");
    checkGridBorders(ymin, ymax, "y");
    checkGridBorders(etamin, etamax, "eta");

    if (icModel == 10)
        tauResize = 100.0;  // do not resize grid in dynIC
}

void Simulation::printParameters() {
    std::cout << "====== parameters ======\n"
         << "outputDir = "        << outputDir      << "\n"
         << "freezeoutOnly = "    << freezeoutOnly  << "\n"
         << "freezeoutExtend = "  << freezeoutExtend<< "\n"
         << "vorticity = "        << vorticityOn    << "\n"
         << "eosType = "          << eosType        << "\n"
         << "eosTypeHadron = "    << eosTypeHadron  << "\n"
         << "nx = "               << nx             << "\n"
         << "ny = "               << ny             << "\n"
         << "nz = "               << nz             << "\n"
         << "icModel = "          << icModel        << "\n"
         << "cartesian = "        << cartesian      << "\n"
         << "xmin = "             << xmin           << "\n"
         << "xmax = "             << xmax           << "\n"
         << "ymin = "             << ymin           << "\n"
         << "ymax = "             << ymax           << "\n"
         << "etamin = "           << etamin         << "\n"
         << "etamax = "           << etamax         << "\n"
         << "tau0 = "             << tau0           << "\n"
         << "tauMax = "           << tauMax         << "\n"
         << "tauGridResize = "    << tauResize      << "\n"
         << "dtau = "             << dtau           << "\n"
         << "e_crit = "           << eCrit          << "\n"
         << "etaS = "             << etaS           << "\n"
         << "zetaS = "            << zetaS          << "\n"
         << "epsilon0 = "         << epsilon0       << "\n"
         << "impactPar = "        << impactPar      << "\n"
         << "s0ScaleFactor = "    << s0ScaleFactor  << "\n"
         << "======= end parameters =======\n";
}