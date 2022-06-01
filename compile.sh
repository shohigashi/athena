#!/bin/bash
make clean &> /dev/null
##python configure.py --prob=collapse_primordial --grav=mg -b --flux=hlld -hdf5 -mpi -h5double --nghost=4 --mpiccmd=mpiicpc -fft --cflag=-static_mpi &> complog.d
python configure.py --prob=collapse_primordial --grav=mg --flux=hllc -hdf5 -mpi -h5double  --mpiccmd=mpicxx --nghost=4 -fft --cflag=static_mpi &> complog.d
#python configure.py --prob=collapse_primordial  --flux=hllc --grav=mg -hdf5 -mpi --mpiccmd=mpiicpc --cflag=-static_mpi &> complog.d
make -j &>> complog.d
