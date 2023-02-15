#module purge
#module unload PrgEnv-cray/6.0.4
module clear
module refresh
module load gnu/pythonpackages/2.7.9
module load PrgEnv-intel/6.0.4
module load intel/odbserver/0.16.2.omp.1
module load gnu/anaconda2/5.3.0
module load gnu/python_dateutil/2.7.3
module list
which python

export MONITOBS=${MONITOBS:-${0%/*}}
export OBSLIB=${MONITOBS}/pylib
export NMLDIR=${MONITOBS}/nml

echo ${OBSLIB}

python ${MONITOBS}/monitobs_gui.py
