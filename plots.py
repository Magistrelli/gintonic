#!/usr/bin/env python3

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# import modules, define constants and functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import os.path
import argparse
import matplotlib.pyplot as plt
import tqdm
from multiprocessing import Pool

import matplotlib
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

HBAR = 6.58212e-22      # Planck reduced constant in units of Mev s
PMASS = 1.87913e3       # 2.*neutron_mass*c^2 in units of Mev
C_LIGHT = 2.99792e10    # speed of light in units of cm/s

std_figsize = [20, 9]                   # default size of the figures
figsize = std_figsize                   # default size of the figures
plt.rcParams.update({'font.size': 18})  # default size of writings
splot_stc = [1, 2]                      # default suplot structure


def vortScatter():
    '''
    Show the initial position of the vortices with red points
    This method prints all the vortices of n x m cells centered on the
        fundamental one

    The position of the fundamental cell's vortices, <xv0> and <yv0>, and the
        dimensions of the cell's sides, <L1> and <L2>, must be specified as
        global variables before using this method!
    <n> and <m> are global variables
    '''

    for i in range(n):
        for j in range(m):
            if (hasattr(vcharge, "__len__")):    # more vortices
                for k in range(len(xv0)):
                    plt.scatter(xv0[k]-(m-1-2*i)*L1/2, yv0[k]-(n-1-2*j)*L2/2,
                                s=10, c='red')
            else:   # single vortex
                plt.scatter(xv0-(m-1-2*i)*L1/2, yv0-(n-1-2*j)*L2/2,
                            s=10, c='red')


def axis_limlab(xlim, ylim):
    '''
    Set the plot's axes limits and labels
    '''

    plt.xlim(xlim)
    plt.ylim(ylim)
    units = ""
    if (args.physical_units):
        units=" [cm]"
    plt.xlabel("x"+units)
    plt.ylabel("y"+units)


def plot_common_ops(xlim, ylim):
    '''
    Operations used in all plots
    '''

    plt.colorbar()
    if (initWithVortices):
        vortScatter()
    axis_limlab(xlim, ylim)


def plotField(Field, title, levels=np.linspace(0., 2.*np.pi, 500), ifvel=1):
    '''
    Plots as a color map the values assumed by a superfluid's field on the
        fundamental cell
    It optionally represents as an arrow map also the velocity field

    The coordinate and velocity arrays, <X>, <Y>, <Velx> and <Vely>, and the
        axes limits, <xlim> and <ylim>, must be specified as global variables
        before using this method!

    :arg: <Field>   array which stores the values of field to be plotted
    :arg: <title>   title of the plot
    :arg: <levels>  'levels' argument of pyplot.contourf
    :arg: <ifvel>   = 1 to plot the velocity field
                    = 0 elsewhere
    '''

    plt.contourf(X, Y, Field, levels=levels)
    plot_common_ops(xlim, ylim)
    if ((ifvel == 1) and (args.narr>0) and
        (np.nanmax(Velx)**2. + np.nanmax(Vely)**2. > 1e-10)):
        plt.quiver(X[::args.narr], Y[::args.narr],
                   Velx[::args.narr, ::args.narr],
                   Vely[::args.narr, ::args.narr],
                   color='k')
    plt.title(title)


def plotQPBCField(Field, title,
                  levels=np.linspace(0., 2.*np.pi, 500), ifvel=1):
    '''
    Plots as a color map the values assumed by a superfluid's field on a
        n x m grid of unit cells
    It optionally represents as an arrow map also the velocity field

    The coordinate and velocity arrays, <XPbc>, <YPbc>, <VelxPbc> and
        <VelyPbc>, and the axes limits, <xlimPbc> and <ylimPbc>, must be
        specified as global variables before using this method!
    <n> and <m> are global variables

    :arg: <Field>   array which stores the values of field to be plotted
    :arg: <title>   title of the plot
    :arg: <levels>  'levels' argument of pyplot.contourf
    :arg: <ifvel>   = 1 to plot the velocity field
                    = 0 elsewhere
    '''

    plt.contourf(XPbc, YPbc, Field, levels=levels)
    plot_common_ops(xlimPbc, ylimPbc)
    if ((ifvel == 1) and (args.narr>0) and
        (np.nanmax(VelxPbc)**2. + np.nanmax(VelyPbc)**2. > 1e-30)):
        plt.quiver(XPbc[::m*args.narr], YPbc[::n*args.narr],
                   VelxPbc[::m*args.narr, ::n*args.narr],
                   VelyPbc[::m*args.narr, ::n*args.narr],
                   color='k')
    plt.title(title)


def mod_average(vals, weight=[-1.,], mod=2.*np.pi):
    '''
    Calculate the modular average of an array of elements

    :arg: <vals>    array to average
    :arg: <mod>     modulus of the values>
    '''

    if (np.max(vals) - np.min(vals) > 0.9*mod):
       for i,val in enumerate(vals):
            if val > mod/2.: vals[i] -= mod

    val_sum = 0.
    if (np.any(weight) < 0.):
        weight=np.full(len(vals), 1.)
    for i,val in enumerate(vals):
        val_sum += val * weight[i]

    return (val_sum/sum(weight)) %mod


def plotPhaseDens(Phase, Dens, time, savename,
        ifKelvin=False, xCirc=None, yCirc=None, ifQpbc=0, ifShow=0, it_cs=-1):
    '''
    Represents a snapshot of the superfluid taken at a certain time
    Plots as color maps, in two different subplots, the phase field and the
        superfluid density
    Plots as an arrow map the velocity field in both subplots

    <args.nlim> is a global variable

    :arg: <time>        time when the snapshot is taken
    :arg: <savename>    the figure will be saved as <savename>
    :arg: <ifQpbc>      = 1 if the grid of <n> x <m> cell is shown
                        = 0 if only the fundamental cell is shown
    :arg: <ifShow>      = 1 to show the plot during the run
                        = 0 elsewhere
    '''

    fig = plt.figure(figsize=figsize)

    plt.subplot(splot_stc[0], splot_stc[1], 1)
    if (ifQpbc == 0):
        plotField(Phase, "Phase")
        if (ifKelvin): plt.scatter(xCirc, yCirc, s=10, c='white')
    else: plotQPBCField(Phase, "Phase Field")

    units = ""
    if (args.physical_units):
        units=" [fm^-2]"

    plt.subplot(splot_stc[0], splot_stc[1], 2)
    if (ifQpbc == 0):
        plotField(Dens, "Superfluid Density"+units,
                  levels=np.linspace(0, args.nlim, 100))
        if (ifKelvin): plt.scatter(xCirc, yCirc, s=10, c='white')
    else:
        plotQPBCField(Dens, "Superfluid Density"+units,
                      levels=np.linspace(0., args.nlim, 100))

    units = ""
    if (args.physical_units):
        units=" s"

    if ((it_cs >= 0) and (initWithVortices)):
        with open( os.path.join(
                    args.input, "Cs_"+str(it_cs)+".out"), 'r' ) as ffile:
            lines = ffile.readlines()
            c1 = [float(x) for x in lines[0].split()]
            c2 = [float(x) for x in lines[1].split()]
        # average of cs: carefull with mod2pi
        c1 = mod_average(c1, weight=Dens[:,0])
        c2 = mod_average(c2, weight=Dens[0,:])
        # Wood+2019 eq 14
        yvort = - L2 * (c1/np.pi - Nv)/ (2.*Nv)
        xvort = L1 * (c2/np.pi - Nv) / (2.*Nv)
        plt.scatter(xvort, yvort, s=10, c='black')

    plt.suptitle("t = " + str("%.2e" % (time)) + units)
    plt.subplots_adjust(right = 1., left = 0.05, wspace=0.2, hspace=0.5)
    plt.savefig(savename)     # , dpi = 200) for more resolution
    if (ifShow == 1):
        plt.show()
    plt.close(fig)


def plotCs(fnumber):
    '''
    Plots the values assumed by the QPBC integration constants along the sides
        of the fundamental cell

    The coordinate arrays, <X> and <Y>, must be specified as global variables
        before using this method!
    <std_figsize> is a global variable

    :arg: <fnumber>     the input file storing the values of the c_i must be
                        named 'os.path.join(args.input, 'Cs_<fnumber>.out')
    '''

    with open( os.path.join(
                args.input, "Cs_"+str(fnumber)+".out"), 'r' ) as ffile:
        lines = ffile.readlines()
        c1 = [float(x) for x in lines[0].split()]
        c2 = [float(x) for x in lines[1].split()]
    fig = plt.figure(figsize=std_figsize)

    units = ""
    if (args.physical_units):
        units=" [cm]"

    plt.subplot(2, 1, 1)
    plt.plot(Y, c1, marker='.')
    plt.xlabel("y"+units)
    plt.ylabel("c1 (y)"+units)

    plt.subplot(2, 1, 2)
    plt.plot(X, c2, marker='.')
    plt.xlabel("x"+units)
    plt.ylabel("c2 (x)"+units)

    plt.subplots_adjust(right = 0.95, left = 0.05, wspace=0.2, hspace=0.5)
    try:
        plt.savefig(os.path.join(args.output, "ci_"+f'{fnumber:03d}'+".jpg"))
    except:
        plt.savefig(os.path.join(args.output, "ci_"+fnumber+".jpg"))
    plt.close(fig)


def plot_f_of_t(time, cons, label, title, res=-1):
    '''
    Plots the time evolution of a superfluid's conserved quantities

    :arg: <cons>    the conserved quantity whose evolution is plotted
    :arg: <label>   label of the y-axis
    :arg: <title>   title of the plot
    :arg: <res>     if >0, the y-axis is centered on the mean value of <cons>
                    and it covers an interval equals to 2.*<res>
    :arg: <time>    time vector; by default, this is specified as a global
                    variable and must thus be specified before using
                    this method!
    '''

    plt.plot(time, cons)
    plt.title(title)
    units = ""
    if (args.physical_units):
        units=" [s]"
    plt.xlabel("t"+units)
    plt.ylabel(label)
    if (res > 0):
        cMean = np.mean(cons)
        plt.ylim(cMean-cMean/res, cMean+cMean/res)


def plotvBgAndNStdDev(vbg_file, ifShow=False):
    '''
    Plots module of vBg against the standard deviation of the density
    '''

    time, vbgx, vbgy, nStdDev = \
            np.genfromtxt(vbg_file, unpack=True)
    vBgMod = np.sqrt(vbgx**2 + vbgy**2)

    fig = plt.figure(figsize=(20, 20))
    plt.plot(vBgMod, nStdDev)
    plt.title("nStdDev")
    units_v, units_n = "", ""
    if (args.physical_units):
        units_v=" [fm/s]"
        units_n=" [1/fm]"
    plt.xlabel("|vBg|"+units_v)
    plt.ylabel("nStdDev"+units_n)
    plt.savefig(os.path.join(args.output, "vBgVsNStdDev.jpg"))
    if ifShow: plt.show()
    plt.close(fig)

    fig = plt.figure(figsize=(20, 20))
    plt.subplot(1, 2, 1)
    plot_f_of_t(time, vbgx, "vBg_x(t)", "")
    plt.subplot(1, 2, 2)
    plot_f_of_t(time, vbgy, "vBg_x(t)", "")
    plt.savefig(os.path.join(args.output, "vbg_t.jpg"))
    plt.close(fig)


def plot_screenshot(args_screenshot):
    i, fig, nfiles_burn, ifKelvin, if_physical_units, \
            lenx, leny, dtout, if_cs, fracCellOut, \
            if_potential, if_particles, vxBg, vyBg, minPin, \
            maxPin, L1, L2, Pin = args_screenshot

    for idname in ("", "burn_-"):
        if (idname=="burn_-"):
            idx_time = -i
            idx = nfiles_burn - i
        else:
            idx_time = i + 1
            idx = i + nfiles_burn + 1
        try:
            x, y, dens, phase, velx, vely = \
                    np.genfromtxt( os.path.join(args.input,
                                        "phase_"+idname+str(i+1)+".out"),
                                    skip_header=3, unpack=True,
                                    usecols=(0, 1, 2, 3, 4, 5) )
        except: continue
        if (ifKelvin):
            xCirc, yCirc = \
                    np.genfromtxt( os.path.join(
                            args.input, "Kelvin_"+idname+str(i+1)+".out"),
                            unpack=True )
            if (if_physical_units):
                xCirc *= CSI
                yCirc *= CSI
        else: xCirc, yCirc = None, None
        Phase = np.transpose(np.reshape(phase, (lenx, leny)))
        Dens = np.transpose(np.reshape(dens, (lenx, leny)))
        Velx = np.transpose(np.reshape(velx, (lenx, leny)))
        Vely = np.transpose(np.reshape(vely, (lenx, leny)))

        # plot phase, density and velocity fields
        plotPhaseDens(Phase, Dens, dtout*idx_time,
                    os.path.join(args.output, f"phase_{idx:04d}.jpg"),
                    ifKelvin=ifKelvin, xCirc=xCirc, yCirc=yCirc, it_cs=-1)

        # plot the c_i constants
        if ((if_cs) and (fracCellOut == 1.)): plotCs(idx)

        # plot the landscape potential
        if if_potential:
            if ((if_particles) or
                    ((vxBg**2+vyBg**2 != 0.) and (minPin != maxPin))):
                x, y, pin = np.genfromtxt( os.path.join(
                                args.input, "pin_"+idname+str(i+1)+".out"),
                            skip_header=1, unpack=True )
                Pin = np.transpose(np.reshape(pin, (lenx, leny)))
                fig = plt.figure(figsize=(figsize[0]*L1/L2*0.5, figsize[1]*0.5))
                if (minPin != maxPin):
                    plotField(Pin, "Pinning Potential [MeV]",
                        levels=np.linspace(minPin, maxPin, 100), ifvel=0)
                else: plotField(Pin, "Pinning Potential [MeV]", ifvel=0)
                plt.savefig(os.path.join(args.output, f"pinning_{idx:04d}.jpg"))
                plt.close(fig)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Plots results from gintonic simulations')
    parser.add_argument('-i', '--input', default="out",
                        help='input folder storing the results of the main code')
    parser.add_argument('-o', '--output', default="jpg",
                        help='output folder storing the plots')
    parser.add_argument('--narr', type=int, default=0,
                        help='prints velocity every narr points of the grid')
    parser.add_argument('--nlim', type=float, default=1.1,
                        help='upper limit of the density axis in plots')
    parser.add_argument('--cs', action='store_true',
                        help='plots QPBC constants')
    parser.add_argument('--potential', action='store_true',
                        help='plots background potential')
    parser.add_argument('--physical_units', action='store_true',
                        help='plots in physical units')
    parser.add_argument('--show', action='store_true',
                        help='show (some) plots during production')
    parser.add_argument('-p', '--nprocs', type=int, default=20,
                        help='Number of processes')
    args = parser.parse_args()

    # %%%%%%%%%%%%%%%%
    # read input files
    # %%%%%%%%%%%%%%%%

    # read the global input parameters
    nbulk, mu, vxBg, vyBg, \
      L1, L2, lenx, leny, dt, tmax, ntout, \
      ptype, ifActivePar, pStr, pSig, px0, py0, npx, npy, \
      initConfig, vort_file, _, _, _, _, _, \
      _, tburn, _, _, ifKelvin, n, m, fracCellOut, outGrid \
            = np.genfromtxt( os.path.join(args.input, "parameters.dat"),
                             dtype=str )
    nbulk = float(nbulk)
    mu = float(mu)
    vxBg = float(vxBg)
    vyBg = float(vyBg)
    L1 = float(L1)
    L2 = float(L2)
    tmax = float(tmax)
    ntout = int(ntout)
    tburn = float(tburn)
    pStr = float(pStr)
    dt = float(dt)
    lenx = int(lenx)
    leny = int(leny)
    ifActivePar = int(ifActivePar)
    initConfig = int(initConfig)
    ifKelvin = int(ifKelvin)
    m = int(m)
    n = int(n)
    fracCellOut = float(fracCellOut)
    outGrid = int(outGrid)

    # if dt < 0 the code uses the CFL condition with Courant number = |dt|
    if (dt < 0.): dt = np.abs(dt) / np.sqrt(mu) / (lenx/L1 + leny/L2)
    ntmax = int(tmax / dt)
    ntburn = int(tburn / dt)

    CSI = HBAR*C_LIGHT/np.sqrt(PMASS*mu)
    T0 = HBAR/mu
    if ((pStr < 0) or (ptype == "cos")):
        args.nlim += 0.5
    elif (ptype == "rod"):
        args.nlim += 0.1
    elif (ifActivePar):
        args.nlim += 0.3
    if (initConfig == 0): initWithVortices = True
    else: initWithVortices = False
    nchuncks = int((ntmax+ntburn)/ntout)
    if (L1 > 2.*L2):
        buffer = splot_stc[0]
        splot_stc[0] = splot_stc[1]
        splot_stc[1] = buffer
        figsize = [20,16]
    else: figsize[0] *= L1/L2

    # compute the output cell dimensions
    lenx = int(int(lenx*fracCellOut)/outGrid)
    leny = int(int(leny*fracCellOut)/outGrid)
    # physical units
    if (args.physical_units):
        L1 *= CSI
        L2 *= CSI
        dt *= T0
        args.nlim = 1.1e-26/(CSI**2)
    # grid limits
    xlim = (-0.5*L1*fracCellOut, 0.5*L1*fracCellOut)
    ylim = (-0.5*L2*fracCellOut, 0.5*L2*fracCellOut)
    # n x m cells grid limits
    if (fracCellOut == 1.):
        xlimPbc = (xlim[0]-(m-1)*L1/2, xlim[1]+(m-1)*L1/2)
        ylimPbc = (ylim[0]-(n-1)*L2/2, ylim[1]+(n-1)*L2/2)

    print(nbulk, mu, vxBg, vyBg, dt, lenx, leny, L1, L2, ptype, ifActivePar, ntburn)

    if (initWithVortices):
        # read the vortex input parameters
        xv0, yv0, vcharge = np.genfromtxt(
                    os.path.join(args.input, vort_file),
                    skip_header=1, unpack=True)
        # physical units
        if (args.physical_units):
            xv0 *= CSI
            yv0 *= CSI
        # compute the total vortices charge
        if (hasattr(vcharge, "__len__")):
            Nv = int(sum(vcharge))  # more vortices
        else:
            Nv = int(vcharge)       # single vortex
        print(xv0, yv0, Nv, "\n")
    else :
        print("No vortices\n")


    # %%%%%%%%%%%%%%%%%%%%%%%
    # plot initial condition
    # %%%%%%%%%%%%%%%%%%%%%%%

    # create jpg folder
    if not os.path.exists(args.output): os.makedirs(args.output)

    # read fundamental cell's density, phase and velocity
    x, y, dens, phase, velx, vely = \
                np.genfromtxt(os.path.join(args.input, "phase_ini.out"),
                      skip_header=3, unpack=True, usecols=(0, 1, 2, 3, 4, 5))
    # if Kelvin is checked, read the Kelvin circuit
    if (ifKelvin):
        xCirc, yCirc = np.genfromtxt(
                    os.path.join(args.input, "Kelvin_0.out"), unpack=True)
        if (args.physical_units):
            xCirc *= CSI
            yCirc *= CSI
    else: xCirc, yCirc = None, None
    # reshape inputs to match a 2D grid
    Dens = np.transpose(np.reshape(dens, (lenx, leny)))
    Phase = np.transpose(np.reshape(phase, (lenx, leny)))
    Velx = np.transpose(np.reshape(velx, (lenx, leny)))
    Vely = np.transpose(np.reshape(vely, (lenx, leny)))
    # remove redundant coordinates
    X = []
    Y = []
    [X.append(xx) for xx in x if xx not in X]
    [Y.append(yy) for yy in y if yy not in Y]
    X = np.array(X)
    Y = np.array(Y)

    # plot fundamental cell's density, phase and velocity
    plotPhaseDens(Phase, Dens, -ntburn*dt,
                os.path.join(args.output, "phase_ini.jpg"),
                ifKelvin=ifKelvin, xCirc=xCirc, yCirc=yCirc, ifShow=args.show)


    maxPin, minPin, Pin = None, None, None
    if args.potential:
        # read the landscape potential
        x, y, pin = \
                np.genfromtxt( os.path.join(args.input, "pin_ini.out"),
                                skip_header=1, unpack=True )
        maxPin = max(pin)   # upper limit of the pinning potential (color) axis
        minPin = min(pin)   # lower limit of the pinning potential (color) axis
        # reshape input to match a 2D grid
        Pin = np.transpose(np.reshape(pin, (lenx, leny)))
        # plot the landscape potential
        fig = plt.figure(figsize=(figsize[0]*L1/L2*0.5, figsize[1]*0.5))
        if (minPin != maxPin):
            plotField(Pin, "Pinning Potential [MeV]",
                levels=np.linspace(minPin, maxPin, 100), ifvel=0)
        else: plotField(Pin, "Pinning Potential [MeV]", ifvel=0)
        plt.savefig(os.path.join(args.output, "pinning_ini.jpg"))  #, dpi = 200)
        if args.show: plt.show()
        plt.close(fig)


    # read and plot the QPBC integration constants
    if ((args.cs) and (fracCellOut == 1.)):
        plotCs("ini")


    # read n x m fundamental cells' density, phase and velocity in QPBC
    if (fracCellOut == 1.):
        xPbc, yPbc, nPbc, phasePbc, rPbc, iPbc, velxPbc, velyPbc = \
                np.genfromtxt( os.path.join(args.input, "iniPbc.out"),
                                skip_header=1, unpack=True )
        # reshape inputs to match a 2D grid
        NPbc = np.transpose(np.reshape(nPbc, (m*lenx, n*leny)))
        PhasePbc = np.transpose(np.reshape(phasePbc, (m*lenx, n*leny)))
        RPbc = np.transpose(np.reshape(rPbc, (m*lenx, n*leny)))
        IPbc = np.transpose(np.reshape(iPbc, (m*lenx, n*leny)))
        VelxPbc = np.transpose(np.reshape(velxPbc, (m*lenx, n*leny)))
        VelyPbc = np.transpose(np.reshape(velyPbc, (m*lenx, n*leny)))
        # remove redundant coordinates
        XPbc = []
        YPbc = []
        [XPbc.append(xx) for xx in xPbc if xx not in XPbc]
        [YPbc.append(yy) for yy in yPbc if yy not in YPbc]
        XPbc = np.array(XPbc)
        YPbc = np.array(YPbc)

        # plot n x m fundamental cells' density, phase and velocity
        plotPhaseDens(PhasePbc, NPbc, -ntburn*dt,
                    os.path.join(args.output, "QPBC_IniPhase.jpg"),
                    ifKelvin=ifKelvin, xCirc=xCirc, yCirc=yCirc,
                    ifQpbc=1, ifShow=args.show)

        units = ""
        if (args.physical_units):
            units=" [fm^-1]"

        # and Real and Imaginary parts of the wavefunction
        fig = plt.figure(figsize=figsize)
        plt.subplot(splot_stc[0], splot_stc[1], 1)
        plotQPBCField(RPbc, "Wavefunction's Real Part"+units,
                      levels=100, ifvel=0)
        plt.subplot(splot_stc[0], splot_stc[1], 2)
        plotQPBCField(IPbc, "Wavefunction's Imaginary Part"+units,
                      levels=100, ifvel=0)
        plt.savefig(os.path.join(args.output, "QPBC_IniRI.jpg"))  # , dpi = 200)
        plt.close(fig)


    # %%%%%%%%%%%%%%%%%%%%%%%%%
    # plot conserved quantities
    # %%%%%%%%%%%%%%%%%%%%%%%%%
    t, Ntot, Etot, vxmean, vymean, circ = \
                np.genfromtxt( os.path.join(args.input, "cons.out"),
                                skip_header=1, unpack=True )

    fig = plt.figure(figsize=(20, 20))

    plt.subplot(3, 2, 1)
    plot_f_of_t(t, Ntot, "N(t)", "Fundamental cell's total particle number",
                  res=1e3)

    plt.subplot(3, 2, 2)
    plot_f_of_t(t, Etot, "E(t) [MeV]", "Mean particle total energy")

    plt.subplot(3, 2, 3)
    plot_f_of_t(t, vxmean, "Px(t) [MeV/c]", "Superfluid's mean x momentum")

    plt.subplot(3, 2, 4)
    plot_f_of_t(t, vymean, "Py(t) [MeV/c]", "Superfluid's mean y momentum")

    plt.subplot(3, 2, 5)
    plot_f_of_t(t, circ, "circ(t) [2pi]", "circulation (cell's boundaries)")

    plt.subplots_adjust(right = 0.8, left = 0.1, wspace=0.2, hspace=0.5)
    plt.savefig(os.path.join(args.output, "conserved.jpg"))
    if args.show: plt.show()
    plt.close(fig)


    # %%%%%%%%%%%%%%%%%%%%%%%%%
    # plot |vBg| vs nStdDev to look for critical velocities
    # %%%%%%%%%%%%%%%%%%%%%%%%%
    if (vxBg**2 + vyBg**2 > 0.):
        plotvBgAndNStdDev( os.path.join(args.input, "vel_and_nStdDev.out"),
                            ifShow=args.show )


    # %%%%%%%%%%%%%%%%%%%
    # plot time evolution
    # %%%%%%%%%%%%%%%%%%%

    nfiles_burn = 0
    for i in range(nchuncks):
      if os.path.isfile( os.path.join(
                        args.input, "phase_burn_-"+str(i+1)+".out") ):
        nfiles_burn += 1
      else: break

    args_screenshot = [(i, fig, nfiles_burn, ifKelvin, args.physical_units,
                        lenx, leny, dt*ntout, args.cs, fracCellOut,
                        args.potential, ifActivePar, vxBg, vyBg, minPin,
                        maxPin, L1, L2, Pin) for i in range(nchuncks)]
    pool = Pool(processes=args.nprocs)
    plots = list( tqdm.tqdm(pool.imap(plot_screenshot, args_screenshot),
                            total=nchuncks, desc="Plotting") )


    # %%%%%%%%%%%%%%%%%%%%%
    # plot final condition
    # %%%%%%%%%%%%%%%%%%%%%

    if (fracCellOut == 1.):
        # read n x m fundamental cells' density, phase and velocity in QPBC
        xPbc, yPbc, nPbc, phasePbc, rPbc, iPbc, velxPbc, velyPbc = \
                    np.genfromtxt( os.path.join(args.input, "endPbc.out"),
                                    skip_header=1, unpack=True )
        # reshape inputs to match a 2D grid
        NPbc = np.transpose(np.reshape(nPbc, (m*lenx, n*leny)))
        PhasePbc = np.transpose(np.reshape(phasePbc, (m*lenx, n*leny)))
        RPbc = np.transpose(np.reshape(rPbc, (m*lenx, n*leny)))
        IPbc = np.transpose(np.reshape(iPbc, (m*lenx, n*leny)))
        VelxPbc = np.transpose(np.reshape(velxPbc, (m*lenx, n*leny)))
        VelyPbc = np.transpose(np.reshape(velyPbc, (m*lenx, n*leny)))
        # plot n x m fundamental cells' density, phase and velocity
        plotPhaseDens(PhasePbc, NPbc, 0,
                os.path.join(args.output, "QPBC_EndPhase.jpg"),
                ifKelvin=ifKelvin, xCirc=xCirc, yCirc=yCirc, ifQpbc=1)
