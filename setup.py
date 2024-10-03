from setuptools import setup, find_packages

setup(
    name='pyiron_io',
    version='0.0.1',
    url='https://github.com/pyiron-dev/pyiron-hdf5-format',
    description='Compare how different pyiron job objects are serialised in HDF5',
    packages=find_packages(),
    install_requires=[
        "pint",
        "numpy",
        "h5io_browser",
        "pyiron_dataclasses"
    ],
)
