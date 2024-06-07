# pyiron HDF5 format
Currently the three major simulation codes in pyiron (LAMMPS, SPHInX and VASP) all implment their own HDF5 format. To address this challenge we aim to develop a minimalistic interface to access the HDF5 files created by pyiron. 

This minimalistic interface should be able to read and write pyiron compatible HDF5 files using the [h5io_browser](https://github.com/h5io/h5io_browser) package.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pyiron-dev/pyiron-hdf5-format/HEAD)
