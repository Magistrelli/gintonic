/*
Script used to simulate the dynamics of simple superfluid configurations in
    order to check the GPE solvers
It can simulate superfluids either starting with simple vortices configurations,
    initialized with uniform or random densiy and uniform or random phase or
    excited with a couple of opposite solitons
To choose your simulation, comment or decomment the specific lines between 43-51

Usage: ./gin [<vBgT> <vBgEvol>]
  :arg: <vBgT>      time scale for the variation of the bg velocity (if < 0, constant vBg)
  :arg: <vBgEv>     type of prescribed evolution for the background velocity
*/


#define VBGT -1.
#define VBGEV ""

#include "./PDE.h"

int main(int argc, char *argv[]) {
// define variables
GPE * test;
ifstream inputPar;
string buffer;
// time variables
double tmax;                // final time
int tout, tcons;            // how often (time steps) to print output and conserved quatities
// initial configuration
unsigned int initConf;      // initial configuration flag
bool ifActiveParticles;
string vortexFile;          // file with vortices initial positions and charges
bool rndN, rndPhase;        // flags for random phase and density
double x01, x02, nfact;     // solitons parameters
// simulation parameters
unsigned int method;        // PDE solver
bool ifKelvin;              // wheter to check Kelvin circulation theorem
int nrows, ncols;           // number of cell printed with QPBC
double fracCellOut;         // fraction of the fundamental cell printed in outputs
unsigned int outGrid;       // fraction of the points of the spatial grid printed in outputs
// background velocity parameters
double vBgT = VBGT;
string vBgEv = VBGEV;
if (argc != 3) {
    cout << "Using default values: vBgT = " \
            << VBGT << ", vBgEv = " << VBGEV << endl;
} else {
    vBgT = atof(argv[1]);
    vBgEv = argv[2];
    cout << "Using imput values: vBgT = " \
            << vBgT << ", vBgEv = " << vBgEv << endl;
}

// read parameter file
inputPar.open(PARAM_file);
for (unsigned int ibuff=0; ibuff<12; ibuff++) {ReadPar(inputPar, buffer);}
ReadPar(inputPar, tmax);
ReadPar(inputPar, tout);
for (unsigned int ibuff=0; ibuff<2; ibuff++) {ReadPar(inputPar, buffer);}
ReadPar(inputPar, ifActiveParticles);
for (unsigned int ibuff=0; ibuff<7; ibuff++) {ReadPar(inputPar, buffer);}
ReadPar(inputPar, initConf);
ReadPar(inputPar, buffer);
ReadPar(inputPar, vortexFile);
ReadPar(inputPar, buffer);
ReadPar(inputPar, rndN);
ReadPar(inputPar, rndPhase);
ReadPar(inputPar, buffer);
ReadPar(inputPar, x01);
ReadPar(inputPar, x02);
ReadPar(inputPar, nfact);
ReadPar(inputPar, buffer);
ReadPar(inputPar, method);
for (unsigned int ibuff=0; ibuff<4; ibuff++) {ReadPar(inputPar, buffer);}
ReadPar(inputPar, ifKelvin);
ReadPar(inputPar, nrows);
ReadPar(inputPar, ncols);
ReadPar(inputPar, fracCellOut);
ReadPar(inputPar, outGrid);
inputPar.close();

// how often to print the conserved quantities
if (tmax < 0.1) {tcons = 1;}
else {tcons = 10;}

// initialize the simulation
if (initConf == 0) {test = new GPE(vortexFile);}
else if (initConf == 1) {test = new GPE(rndN, rndPhase);}
else if (initConf == 2) {test = new GPE(x01, x02, nfact);}
else {
    cout << endl << "ERROR! Initial configuration not valid!" << endl;
    exit(-10);
}

test->SetOutCell(fracCellOut, outGrid);
test->SetTSkip(tcons);
test->SetVbgT(vBgT);
test->SetVbgEvol(vBgEv);

// initialize support arrays (used only if Kelvin is on)
if (ifKelvin) {test->SetKelvin();}

// print out the initial condition
test->OutputPhase("phase_ini.out", 0);
test->OutputPin("pin_ini.out");
if (ifActiveParticles) {
    test->OutputParticles("particles_ini.out");
}
if (fracCellOut == 1.) {
    test->PrintAllGrid(nrows, ncols, "iniPbc.out");
    test->OutCs("Cs_ini.out");
}
cout << endl;

// solves the GPE (applies <method> a until <tmax> is reached)
test->Solver(tmax, tout, method);
// FIXME something during the evolution with active particles does NOT work (compare with old active particles results)

// print out the final condition
if (fracCellOut == 1.) {
    test->PrintAllGrid(nrows, ncols, "endPbc.out");
}

cout << endl;
delete test;
return 0;
}
