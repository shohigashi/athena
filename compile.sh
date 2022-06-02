#!/bin/bash
make clean &> /dev/null
### without Magnetic field
#python configure.py --prob=collapse_primordial --grav=mg  --flux=hllc -hdf5 -mpi -h5double  --mpiccmd=mpicxx --nghost=4 -fft &> complog.d
### with Magnetic field
python configure.py --prob=collapse_primordial --grav=mg -b --flux=hlld -hdf5 -mpi -h5double  --mpiccmd=mpicxx --nghost=4 -fft &> complog.d
make -j &>> complog.d
