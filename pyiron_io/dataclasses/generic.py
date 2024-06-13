from dataclasses import dataclass
import numpy as np
from typing import Optional, List


@dataclass
class OutputGenericDFT:
    energy_free: np.ndarray
    n_valence: dict
    bands_k_weights: np.ndarray
    kpoints_cartesian: np.ndarray
    bands_e_fermi: np.ndarray
    bands_occ: np.ndarray
    bands_eigen_values: np.ndarray
    scf_convergence: List[bool]
    scf_energy_int: np.ndarray
    scf_energy_free: np.ndarray
    scf_computation_time: np.ndarray
    scf_energy_zero: np.ndarray
    scf_energy_band: np.ndarray
    scf_electronic_entropy: np.ndarray
    scf_residue: np.ndarray
    energy_int: np.ndarray
    computation_time: np.ndarray
    energy_zero: np.ndarray
    energy_band: np.ndarray
    electronic_entropy: np.ndarray
    residue: np.ndarray


@dataclass
class GenericOutput:
    cells: np.ndarray  # N_steps * 3 *3  [Angstrom]
    energy_pot: np.ndarray  # N_steps  [eV]
    energy_tot: np.ndarray  # N_steps  [eV]
    forces: np.ndarray  # N_steps * N_atoms * 3  [eV/Angstrom]
    indices: Optional[np.ndarray]  # N_steps * N_atoms
    natoms: Optional[np.ndarray]  # N_steps
    positions: np.ndarray  # N_steps * N_atoms * 3  [Angstrom]
    pressures: Optional[np.ndarray]  # N_steps * 3 * 3
    steps: Optional[np.ndarray]  # N_steps
    temperature: Optional[np.ndarray]  # N_steps
    unwrapped_positions: Optional[np.ndarray]  # N_steps * N_atoms * 3  [Angstrom]
    velocities: Optional[np.ndarray]  # N_steps * N_atoms * 3  [Angstrom/fs]
    volume: np.ndarray  # N_steps
    dft: Optional[OutputGenericDFT]


@dataclass
class Server:
    user: str
    host: str
    run_mode: str
    queue: Optional[str]
    qid: Optional[int]
    cores: int
    threads: int
    new_h5: bool
    structure_id: Optional[int]
    run_time: int  # [seconds]
    memory_limit: Optional[str]
    accept_crash: bool


@dataclass
class Executable:
    version: str
    name: str
    operation_system_nt: bool
    executable: Optional[str]
    mpi: bool
    accepted_return_codes: List[int]


@dataclass
class GenericDict:
    restart_file_list: list
    restart_file_dict: dict
    exclude_nodes_hdf: list
    exclude_groups_hdf: list


@dataclass
class Interactive:
    interactive_flush_frequency: int
    interactive_write_frequency: int


@dataclass
class GenericInput:
    calc_mode: str
    structure: str
    fix_symmetry: Optional[bool]
    k_mesh_spacing: Optional[float]
    k_mesh_center_shift: Optional[np.ndarray]
    reduce_kpoint_symmetry: Optional[bool]
    restart_for_band_structure: Optional[bool]
    path_name: Optional[str]
    n_path: Optional[str]
    fix_spin_constraint: Optional[bool]
    max_iter: Optional[int]
    temperature: Optional[float]
    n_ionic_steps: Optional[int]
    n_print: Optional[int]
    temperature_damping_timescale: Optional[float]
    pressure_damping_timescale: Optional[float]
    time_step: Optional[int]


@dataclass
class Units:
    length: str
    mass: str


@dataclass
class Cell:
    cell: np.ndarray  # 3 * 3   [Angstrom]
    pbc: np.ndarray  # 3


@dataclass
class Structure:
    dimension: int
    indices: np.array
    info: dict
    positions: np.ndarray  # N_atoms * 3  [Angstrom]
    species: List[str]
    cell: Cell
    units: Units
