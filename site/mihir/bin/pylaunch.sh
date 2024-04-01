PYSCRYPT=$(realpath $1)
shift
ARGS=$@

module load gnu/python/3.9.1
module load gnu/pandas/0.18.1
module load gnu/utility/pynio/1.5.5
module load gnu/pyngl/1.6.1
module load gnu/netcdf4/1.5.3-netcdf4.6.0hdf51.10.0
module load gnu/matplotlib/2.2.2
module load gnu/basemap/1.1.0
module load gnu/python_eccodes/eccodes-2.13.1_utility
module load gnu/py-ncepbufr/1.1.1
module load gnu/pathlib/1.0.1
module load gnu/cython/0.28

python ${PYSCRYPT} ${ARGS}

#display /home/gibies/plots/research/osse_buoy/outfile/surface.png
