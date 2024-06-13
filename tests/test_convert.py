import os
from unittest import TestCase

from h5io_browser.base import read_nested_dict_from_hdf
from pint import UnitRegistry

from pyiron_io.convert import convert_sphinx_job_dict, convert_lammps_job_dict, convert_vasp_job_dict


class TestConvert(TestCase):
    def test_sphinx(self):
        ureg = UnitRegistry()
        job_dict = read_nested_dict_from_hdf(
            file_name=os.path.join(os.path.dirname(__file__), "static", "sx.h5"),
            h5_path="/sx",
            recursive=True,
            slash='ignore',
        )
        job_sphinx = convert_sphinx_job_dict(job_dict=job_dict)
        self.assertEqual(job_sphinx.calculation_output.generic.energy_tot[-1], -228.7831594379917 * ureg.eV)

    def test_lammps(self):
        ureg = UnitRegistry()
        job_dict = read_nested_dict_from_hdf(
            file_name=os.path.join(os.path.dirname(__file__), "static", "lmp.h5"),
            h5_path="/lmp",
            recursive=True,
            slash='ignore',
        )
        job_lammps = convert_lammps_job_dict(job_dict=job_dict)
        self.assertEqual(job_lammps.calculation_output.generic.energy_tot[-1], -9428.45286561574 * ureg.eV)

    def test_vasp(self):
        ureg = UnitRegistry()
        job_dict = read_nested_dict_from_hdf(
            file_name=os.path.join(os.path.dirname(__file__), "static", "vasp.h5"),
            h5_path="/vasp",
            recursive=True,
            slash='ignore',
        )
        job_vasp = convert_vasp_job_dict(job_dict=job_dict)
        self.assertEqual(job_vasp.calculation_output.generic.energy_tot[-1], -14.7459202 * ureg.eV)
