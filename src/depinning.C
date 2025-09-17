/*
Script used to simulate the depinning of a vortex from a specified pinning potential due
    to a vertical background flow
The critical unpinning velocity is looked for by running the simulation many times with
    an increasing background flow
The script also estimates the vortex mean velocity and saves it in 'out/vortexVel.out'

Usage: ./Test.x <t_max> <t_out>
  :arg: <t_max>     total number of temporal steps (in natural units)
  :arg: <t_out>     number of steps after which the intermediate ouput files are printed

FIXME this script is outdated. The structure is correct, but e.g. the reading
  of the parameter file (here called input.dat) is not up to date
*/


#include "./PDE.h"

int main(int argc, char *argv[]) {
if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <t_max> <t_out>" << endl;
    return -1;
}

int tmax = static_cast<int>(atof(argv[1]));
int tout = static_cast<int>(atof(argv[2]));
ifstream input;
ofstream outMeanV;
GPE* sf;
double vscale;          // order of the lag velocity
double vxBg, vyBg;
double xvinit, yvinit;  // initial vortex position
double xvend, yvend;    // final vortex position
double vmeanx, vmeany;  // mean vortex velocity
// compute the total evolution time (n.u.)
double time = 0.0005*1e5;

input.open("input.dat");
// burn input until vyBg and save it in vscale
for (int in = 0; in < 4; in++) {input >> vscale;}
input.close();

outMeanV.open("out/vortexVel.out");
vxBg = 0.;
for (int iv = 0; iv < 10; iv++) {
    vyBg = (iv+1)*2.*vscale;
    cout << endl << "\tRunning the "+to_string(iv)+"-th simulation, vxBg = " \
                 << vxBg << ", vyBg = " << vyBg << endl;

    sf = new GPE("inVort.dat");
    sf->SetVbg(vxBg, vyBg);
    if (iv == 0) {
        sf->GetIniPosV(xvinit, yvinit, 0);
        sf->OutputPhase("ini.out");
        sf->OutputPin("pin.out");
        //sf->OutCs("Cs_0.out");
       }

    cout << endl;
    sf->Solver(tmax, tout, 0, "_"+to_string(iv));
    sf->OutputPhase("end_"+to_string(iv)+".out");

    sf->GetPosV(xvend, yvend, 0);
    vmeanx = (xvend-xvinit)/time;
    vmeany = (yvend-yvinit)/time;

    outMeanV << left << setw(15) << vyBg << setw(15) << vmeanx << setw(15) \
                                                     << vmeany << endl;
    delete sf;
}
outMeanV.close();

cout << endl;
return 0;
}
