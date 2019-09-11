#!/usr/bin/env python3
from __future__ import absolute_import, division, print_function, unicode_literals
#import matplotlib
import adios2
import argparse
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import decomp
import time
import os

#matplotlib.use('Agg')
is_init = 0

def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--instream", "-i", help="Name of the input stream", required=True)
    parser.add_argument("--outfile", "-o", help="Name of the output file", default="screen")
    parser.add_argument("--varname", "-v", help="Name of variable read", default="U")
    parser.add_argument("--nompi", "-nompi", help="ADIOS was installed without MPI", action="store_true")
    parser.add_argument("--displaysec", "-dsec", help="Float representing gap between plot window refresh", default=0.1)
    parser.add_argument("--nx", "-nx", help="Integer representing process decomposition in the x direction",default=1)
    parser.add_argument("--ny", "-ny", help="Integer representing process decomposition in the y direction",default=1)
    parser.add_argument("--nz", "-nz", help="Integer representing process decomposition in the z direction",default=1)
    parser.add_argument("--plane", "-plane", help="The 2D plane to be displayed/stored xy/yz/xz/all", default='yz')
    args = parser.parse_args()

    args.displaysec = float(args.displaysec)
    args.nx = int(args.nx)
    args.ny = int(args.ny)
    args.nz = int(args.nz)

    if args.plane not in ('xz', 'yz', 'xy', 'all'):
        raise "Input argument --plane must be one of xz/yz/xy/all"

    return args


def Plot2D(plane_direction, data, args, fullshape, step, fontsize):
    # Plotting part
    displaysec = args.displaysec
    gs = gridspec.GridSpec(1, 1)
    fig = plt.figure(1, figsize=(8,10))
    ax = fig.add_subplot(gs[0, 0])
    colorax = ax.imshow(data, origin='lower', interpolation='quadric',extent=[0, fullshape[1], 0, fullshape[0]], cmap=plt.get_cmap('gist_ncar'))
    fig.colorbar(colorax, orientation='horizontal')

    for i in range(args.ny):
        y = fullshape[0] / args.ny * i
        ax.plot([0, fullshape[1]], [y, y], color='black')

    for i in range(args.nx):
        x = fullshape[1] / args.nx * i
        ax.plot([x, x], [0, fullshape[0]], color='black')

    ax.set_title("{0} plane, Timestep {1}".format(plane_direction, step), fontsize=fontsize)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.ion()
    global is_init
    print ("Step: {0}".format(step))
    if (args.outfile == "screen"):
        plt.show()
        plt.pause(displaysec)
    elif args.outfile.endswith(".bp"):
        if is_init == 0:
            global adios
            global ioWriter
            global var
            global writer
            is_init = 1 
            adios = adios2.ADIOS(mpi.comm_app)
            ioWriter = adios.DeclareIO("VizOutput")
            var = ioWriter.DefineVariable(args.varname, data.shape, [0,0], data.shape, adios2.ConstantDims, data)
            writer = ioWriter.Open(args.outfile, adios2.Mode.Write)

        writer.BeginStep()
        writer.Put(var, data) #, adios2.Mode.Sync)
        writer.EndStep()
    else:
        imgfile = args.outfile+"{0:0>3}".format(step)+"_" + plane_direction + ".png"
        fig.savefig(imgfile)

    plt.clf()


def read_data(args, fr, start_coord, size_dims):
    
    var1 = args.varname
    data= fr.read(var1, start_coord, size_dims )
    data = np.squeeze(data)
    return data


if __name__ == "__main__":
    # fontsize on plot
    fontsize = 22

    args = SetupArgs()
#    print(args)

    # Setup up 2D communicators if MPI is installed
    mpi = decomp.MPISetup(args, 3)

    # Read the data from this object
    fr = adios2.open(args.instream, "r", mpi.comm_app,"adios2.xml", "VizInput")
#    vars_info = fr.availablevariables()


    # Get the ADIOS selections -- equally partition the data if parallelization is requested
 

    # Read through the steps, one at a time
    k = 0
    k_start = 5
    k_end = 10
    for fr_step in fr:
        k = k + 1
        #if fr_step.current_step()
        start, size, fullshape = mpi.Partition_3D_3D(fr, args)
        cur_step= fr_step.current_step()
        vars_info = fr.available_variables()
        print (vars_info)
        shape3_str = vars_info[args.varname]["Shape"].split(',')
        shape3 = list(map(int,shape3_str))
        data = read_data (args, fr_step, [0,0,0], [shape3[0],shape3[1],shape3[2]])
        fr_step.end_step() 

        #if k < k_start:
        #   continue 

        #if k > k_end:
        #    k_start = k_end + 5
        #    k_end = k_start + 5
        #    continue

        data_fft = np.fft.fftn(data) 
        data_shft = np.abs(np.fft.fftshift(data_fft))
        if args.plane in ('xy', 'all'):
            data_xy = data_shft[:,:,int(shape3[2]/2):int(shape3[2]/2)+1]
            data_xy = np.squeeze(data_xy)
            Plot2D ('xy', data_xy, args, fullshape, cur_step, fontsize)
        
        if args.plane in ('xz', 'all'):
            data_xz = data_shft[:,int(shape3[1]/2):int(shape3[1]/2)+1,:]
            data_xz = np.squeeze(data_xz)
            Plot2D ('xz', data_xz, args, fullshape, cur_step, fontsize)
        
        if args.plane in ('yz', 'all'):
            data_yz = data_shft[int(shape3[0]/2): int(shape3[0]/2)+ 1, :, :]
            data_yz = np.squeeze(data_yz)
            Plot2D ('yz',  data_yz, args, fullshape, cur_step, fontsize)

    fr.close()

