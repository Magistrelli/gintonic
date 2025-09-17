#!/usr/bin/env python3

import numpy as np
import argparse, os

parser = argparse.ArgumentParser(
    description='Generates a series of inVort.dat parameter files for a specified array of x and y position (single charged vortices)')
parser.add_argument('--L1', type=float, default=40.,
    help='x-dimension of the fundamental cell')
parser.add_argument('--L2', type=float, default=40.,
    help='y-dimension of the fundamental cell')
parser.add_argument('--npos', type=int, default=6,
    help='lenght of the position arrays (will define a grid on a triangle)')
parser.add_argument('--charge', type=float, default=1.,
    help='vortex charge')
parser.add_argument('-o', '--outfold', default="noDiss",
    help='output folder')
parser.add_argument('--simmetry', action="store_true",
    help='consider only the bottom-left part of the cell with y<=x<=0')
args = parser.parse_args()

if args.simmetry:
    xvs = np.linspace(-args.L1/2., 0., args.npos)
    yvs = np.linspace(-args.L2/2., 0., args.npos)
    xyvs = [(x, y) for x in xvs for y in yvs if y<=x]
else:
    xvs = np.linspace(-args.L1/2., args.L1/2., args.npos*2-1)[:-1]
    yvs = np.linspace(-args.L2/2., args.L2/2., args.npos*2-1)[:-1]
    xyvs = [(x, y) for x in xvs for y in yvs]

if not os.path.isdir(args.outfold): os.mkdir(args.outfold)
for ip,xy in enumerate(xyvs):
    fold = os.path.join(args.outfold, f"sim_{ip:03d}")
    if not os.path.isdir(fold): os.mkdir(fold)
    outfold, jpgfold = os.path.join(fold, "out"), os.path.join(fold, "jpg")
    if not os.path.isdir(outfold): os.mkdir(outfold)
    if not os.path.isdir(jpgfold): os.mkdir(jpgfold)
    fout = os.path.join(fold, f"inVort.dat")
    with open(fout, "w") as f:
        f.write("xv0[csi_units]\tyv0[csi_units]\tvortex_charge\n\n")
        f.write(str(xy[0])+"\t\t\t"+str(xy[1])+"\t\t\t"+str(args.charge)+"\n")
