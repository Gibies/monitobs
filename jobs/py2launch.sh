PYSCRYPT=$(realpath $1)
shift
ARGS=$@

#gnu/python/3.9.1
module load pbs
module load craype-broadwell
module load cray-snplauncher
module load gnu/pythonpackages/2.7.9
module load gnu/packagesuite/1
module load gnu/pandas/0.18.1
module load gnu/matplotlib/2.2.2
module load gnu/basemap/1.1.0
module load gnu/python_eccodes/eccodes-2.13.1_utility
module load gnu/py-ncepbufr/1.1.1
module load gnu/pathlib/1.0.1
module load gnu/cython/0.28
module load gnu/iris/2.0.0  
module load gnu/toolz/0.9.0
module load gnu/lib/udunits/2.2.26
module load gnu/cartopy/0.16.0
module load gnu/utility/pynio/1.5.5
module load gnu/pyngl/1.6.1
module load gnu/pyke/1.1.1
module load gnu/xarray/0.11.3

module list

python ${PYSCRYPT} ${ARGS}

#display /home/gibies/plots/research/osse_buoy/outfile/surface.png
