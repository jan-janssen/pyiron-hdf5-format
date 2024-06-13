from pint import UnitRegistry

from pyiron_io.hdf import (
    convert_datacontainer_to_dictionary,
    convert_generic_parameters_to_dictionary,
    convert_generic_parameters_to_string,
)
from pyiron_io.dataclasses.generic import (
    Executable,
    Interactive,
    GenericInput,
    GenericOutput,
    OutputGenericDFT,
    GenericDict,
    Server,
    Structure,
    Units,
    Cell,
)
from pyiron_io.dataclasses.lammps import (
    LammpsJob,
    LammpsInput,
    LammpsOutput,
    LammpsPotential,
    LammpsInputFiles,
)
from pyiron_io.dataclasses.sphinx import (
    SphinxJob,
    SphinxInput,
    SphinxInputParameters,
    SphinxRho,
    SphinxAtom,
    SphinxMain,
    SphinxBasis,
    SphinxWaves,
    SphinxOutput,
    SphinxKpoint,
    SphinxElement,
    SphinxStructure,
    SphinxRicQN,
    SphinxInternalInput,
    Species,
    SphinxInitialGuess,
    SphinxPawHamiltonian,
    SphinxPreConditioner,
    SphinxElectronicStructure,
    SphinxElectrostaticPotential,
    SphinxDensityOfStates,
    ScfDiag,
    SphinxOutputChargeDensity,
    PawPot,
    BornOppenheimer,
)


def convert_sphinx_job_dict(job_dict: dict) -> SphinxJob:
    ureg = UnitRegistry()
    sphinx_input_parameter_dict = convert_datacontainer_to_dictionary(
        data_container_dict=job_dict["input"]["parameters"]
    )
    generic_input_dict = convert_generic_parameters_to_dictionary(
        generic_parameter_dict=job_dict["input"]["generic"],
    )
    output_dict = convert_datacontainer_to_dictionary(
        data_container_dict=job_dict["output"]["generic"]
    )
    return SphinxJob(
        executable=Executable(
            version=job_dict["executable"]["executable"]["version"],
            name=job_dict["executable"]["executable"]["name"],
            operation_system_nt=job_dict["executable"]["executable"][
                "operation_system_nt"
            ],
            executable=job_dict["executable"]["executable"]["executable"],
            mpi=job_dict["executable"]["executable"]["mpi"],
            accepted_return_codes=job_dict["executable"]["executable"][
                "accepted_return_codes"
            ],
        ),
        server=Server(
            user=job_dict["server"]["user"],
            host=job_dict["server"]["host"],
            run_mode=job_dict["server"]["run_mode"],
            queue=job_dict["server"]["queue"],
            qid=job_dict["server"]["qid"],
            cores=job_dict["server"]["cores"],
            threads=job_dict["server"]["threads"],
            new_h5=job_dict["server"]["new_h5"],
            structure_id=job_dict["server"]["structure_id"],
            run_time=job_dict["server"]["run_time"],
            memory_limit=job_dict["server"]["memory_limit"],
            accept_crash=job_dict["server"]["accept_crash"],
        ),
        calculation_input=SphinxInput(
            generic_dict=GenericDict(
                restart_file_list=job_dict["input"]["generic_dict"][
                    "restart_file_list"
                ],
                restart_file_dict=job_dict["input"]["generic_dict"][
                    "restart_file_dict"
                ],
                exclude_nodes_hdf=job_dict["input"]["generic_dict"][
                    "exclude_nodes_hdf"
                ],
                exclude_groups_hdf=job_dict["input"]["generic_dict"][
                    "exclude_groups_hdf"
                ],
            ),
            interactive=Interactive(
                interactive_flush_frequency=job_dict["input"]["interactive"][
                    "interactive_flush_frequency"
                ],
                interactive_write_frequency=job_dict["input"]["interactive"][
                    "interactive_flush_frequency"
                ],
            ),
            generic=GenericInput(
                calc_mode=generic_input_dict["calc_mode"],
                structure=generic_input_dict["structure"],
                fix_symmetry=generic_input_dict.get("fix_symmetry", None),
                k_mesh_spacing=generic_input_dict.get("k_mesh_spacing", None),
                k_mesh_center_shift=generic_input_dict.get("k_mesh_center_shift", None),
                reduce_kpoint_symmetry=generic_input_dict.get(
                    "reduce_kpoint_symmetry", None
                ),
                restart_for_band_structure=generic_input_dict.get(
                    "restart_for_band_structure", None
                ),
                path_name=generic_input_dict.get("path_name", None),
                n_path=generic_input_dict.get("n_path", None),
                fix_spin_constraint=generic_input_dict.get("fix_spin_constraint", None),
                max_iter=generic_input_dict.get("max_iter", None),
                temperature=generic_input_dict.get("temperature", None),
                n_ionic_steps=generic_input_dict.get("n_ionic_steps", None),
                n_print=generic_input_dict.get("n_print", None),
                temperature_damping_timescale=generic_input_dict.get(
                    "temperature_damping_timescale", None
                ),
                pressure_damping_timescale=generic_input_dict.get(
                    "pressure_damping_timescale", None
                ),
                time_step=generic_input_dict.get("time_step", None),
            ),
            parameters=SphinxInputParameters(
                sphinx=SphinxInternalInput(
                    paw_pot=PawPot(
                        species=[
                            Species(
                                name=s["name"],
                                pot_type=s["potType"],
                                element=s["element"],
                                potential=s["potential"],
                            )
                            for s in sphinx_input_parameter_dict["sphinx"]["pawPot"][
                                "species"
                            ]
                        ]
                    ),
                    structure=SphinxStructure(
                        cell=sphinx_input_parameter_dict["sphinx"]["structure"]["cell"],
                        species=[
                            SphinxElement(
                                element=s["element"],
                                atom=[
                                    SphinxAtom(
                                        label=a["label"],
                                        coords=a["coords"],
                                        movable=bool(a["movable"]),
                                    )
                                    for a in s["atom"]
                                ],
                            )
                            for s in sphinx_input_parameter_dict["sphinx"]["structure"][
                                "species"
                            ]
                        ],
                    ),
                    basis=SphinxBasis(
                        e_cut=sphinx_input_parameter_dict["sphinx"]["basis"]["eCut"],
                        k_point=SphinxKpoint(
                            coords=sphinx_input_parameter_dict["sphinx"]["basis"][
                                "kPoint"
                            ]["coords"],
                            weight=sphinx_input_parameter_dict["sphinx"]["basis"][
                                "kPoint"
                            ]["weight"],
                            relative=bool(
                                sphinx_input_parameter_dict["sphinx"]["basis"][
                                    "kPoint"
                                ]["relative"]
                            ),
                        ),
                        folding=sphinx_input_parameter_dict["sphinx"]["basis"][
                            "folding"
                        ],
                        save_memory=bool(
                            sphinx_input_parameter_dict["sphinx"]["basis"]["saveMemory"]
                        ),
                    ),
                    paw_hamilton=SphinxPawHamiltonian(
                        number_empty_states=sphinx_input_parameter_dict["sphinx"][
                            "PAWHamiltonian"
                        ]["nEmptyStates"],
                        ekt=sphinx_input_parameter_dict["sphinx"]["PAWHamiltonian"][
                            "ekt"
                        ],
                        methfessel_paxton=bool(
                            sphinx_input_parameter_dict["sphinx"]["PAWHamiltonian"][
                                "MethfesselPaxton"
                            ]
                        ),
                        xc=sphinx_input_parameter_dict["sphinx"]["PAWHamiltonian"][
                            "xc"
                        ],
                        spin_polarized=sphinx_input_parameter_dict["sphinx"][
                            "PAWHamiltonian"
                        ]["spinPolarized"],
                    ),
                    initial_guess=SphinxInitialGuess(
                        waves=SphinxWaves(
                            paw_basis=sphinx_input_parameter_dict["sphinx"][
                                "initialGuess"
                            ]["waves"]["pawBasis"],
                            lcao=sphinx_input_parameter_dict["sphinx"]["initialGuess"][
                                "waves"
                            ]["lcao"],
                        ),
                        rho=SphinxRho(
                            atomic_orbitals=sphinx_input_parameter_dict["sphinx"][
                                "initialGuess"
                            ]["rho"]["atomicOrbitals"]
                        ),
                        no_waves_storage=bool(
                            sphinx_input_parameter_dict["sphinx"]["initialGuess"][
                                "noWavesStorage"
                            ]
                        ),
                    ),
                    main=SphinxMain(
                        ric_qn=SphinxRicQN(
                            max_steps=int(
                                sphinx_input_parameter_dict["sphinx"]["main"]["ricQN"][
                                    "maxSteps"
                                ]
                            ),
                            max_step_length=float(
                                sphinx_input_parameter_dict["sphinx"]["main"]["ricQN"][
                                    "maxStepLength"
                                ]
                            ),
                            born_oppenheimer=BornOppenheimer(
                                scf_diag=ScfDiag(
                                    rho_mixing=float(
                                        sphinx_input_parameter_dict["sphinx"]["main"][
                                            "ricQN"
                                        ]["bornOppenheimer"]["scfDiag"]["rhoMixing"]
                                    ),
                                    spin_mixing=float(
                                        sphinx_input_parameter_dict["sphinx"]["main"][
                                            "ricQN"
                                        ]["bornOppenheimer"]["scfDiag"]["spinMixing"]
                                    ),
                                    delta_energy=sphinx_input_parameter_dict["sphinx"][
                                        "main"
                                    ]["ricQN"]["bornOppenheimer"]["scfDiag"]["dEnergy"],
                                    max_steps=sphinx_input_parameter_dict["sphinx"][
                                        "main"
                                    ]["ricQN"]["bornOppenheimer"]["scfDiag"][
                                        "maxSteps"
                                    ],
                                    preconditioner=SphinxPreConditioner(
                                        type=sphinx_input_parameter_dict["sphinx"][
                                            "main"
                                        ]["ricQN"]["bornOppenheimer"]["scfDiag"][
                                            "preconditioner"
                                        ][
                                            "type"
                                        ],
                                        scaling=sphinx_input_parameter_dict["sphinx"][
                                            "main"
                                        ]["ricQN"]["bornOppenheimer"]["scfDiag"][
                                            "preconditioner"
                                        ][
                                            "scaling"
                                        ],
                                        spin_scaling=sphinx_input_parameter_dict[
                                            "sphinx"
                                        ]["main"]["ricQN"]["bornOppenheimer"][
                                            "scfDiag"
                                        ][
                                            "preconditioner"
                                        ][
                                            "spinScaling"
                                        ],
                                    ),
                                    block_ccg=sphinx_input_parameter_dict["sphinx"][
                                        "main"
                                    ]["ricQN"]["bornOppenheimer"]["scfDiag"][
                                        "blockCCG"
                                    ],
                                ),
                            ),
                        ),
                    ),
                ),
                encut=float(sphinx_input_parameter_dict["EnCut"]),
                kpointcoords=sphinx_input_parameter_dict["KpointCoords"],
                kpointfolding=sphinx_input_parameter_dict["KpointFolding"],
                empty_states=int(sphinx_input_parameter_dict["EmptyStates"]),
                methfessel_paxton=bool(sphinx_input_parameter_dict["MethfesselPaxton"]),
                sigma=float(sphinx_input_parameter_dict["Sigma"]),
                xcorr=sphinx_input_parameter_dict["Xcorr"],
                vasppot=bool(sphinx_input_parameter_dict["VaspPot"]),
                e_step=int(sphinx_input_parameter_dict["Estep"]),
                ediff=float(sphinx_input_parameter_dict["Ediff"]),
                write_waves=bool(sphinx_input_parameter_dict["WriteWaves"]),
                kj_xc=bool(sphinx_input_parameter_dict["KJxc"]),
                save_memory=bool(sphinx_input_parameter_dict["SaveMemory"]),
                rho_mixing=float(sphinx_input_parameter_dict["rhoMixing"]),
                spin_mixing=float(sphinx_input_parameter_dict["spinMixing"]),
                rho_residual_scaling=float(
                    sphinx_input_parameter_dict["rhoResidualScaling"]
                ),
                spin_residual_scaling=float(
                    sphinx_input_parameter_dict["spinResidualScaling"]
                ),
                check_overlap=bool(sphinx_input_parameter_dict["CheckOverlap"]),
                threads=bool(sphinx_input_parameter_dict["THREADS"]),
                use_on_the_fly_cg_optimization=bool(
                    sphinx_input_parameter_dict["use_on_the_fly_cg_optimization"]
                ),
                ionic_step=int(sphinx_input_parameter_dict["Istep"]),
            ),
            structure=Structure(
                dimension=job_dict["input"]["structure"]["dimension"],
                indices=job_dict["input"]["structure"]["indices"],
                info=job_dict["input"]["structure"]["info"],
                positions=job_dict["input"]["structure"]["positions"],
                species=job_dict["input"]["structure"]["species"],
                cell=Cell(
                    cell=job_dict["input"]["structure"]["cell"]["cell"],
                    pbc=job_dict["input"]["structure"]["cell"]["pbc"],
                ),
                units=Units(
                    length=job_dict["input"]["structure"]["units"]["length"],
                    mass=job_dict["input"]["structure"]["units"]["mass"],
                ),
            ),
        ),
        calculation_output=SphinxOutput(
            charge_density=SphinxOutputChargeDensity(
                total=job_dict["output"]["charge_density"]["total"]
            ),
            electronic_structure=SphinxElectronicStructure(
                efermi=job_dict["output"]["electronic_structure"]["efermi"],
                eig_matrix=job_dict["output"]["electronic_structure"]["eig_matrix"],
                k_points=job_dict["output"]["electronic_structure"]["k_points"],
                k_weights=job_dict["output"]["electronic_structure"]["k_weights"],
                occ_matrix=job_dict["output"]["electronic_structure"]["occ_matrix"],
                dos=SphinxDensityOfStates(
                    energies=job_dict["output"]["electronic_structure"]["dos"][
                        "energies"
                    ],
                    int_densities=job_dict["output"]["electronic_structure"]["dos"][
                        "int_densities"
                    ],
                    tot_densities=job_dict["output"]["electronic_structure"]["dos"][
                        "tot_densities"
                    ],
                ),
            ),
            electrostatic_potential=SphinxElectrostaticPotential(
                total=job_dict["output"]["electrostatic_potential"]["total"]
            ),
            generic=GenericOutput(
                cells=output_dict["cells"] * ureg.angstrom,
                energy_pot=output_dict["energy_pot"] * ureg.eV,
                energy_tot=output_dict["energy_tot"] * ureg.eV,
                forces=output_dict["forces"] * ureg.eV / ureg.angstrom,
                indices=output_dict.get("indices", None),
                natoms=output_dict.get("natoms", None),
                positions=output_dict["positions"] * ureg.angstrom,
                pressures=output_dict.get("pressures", None),
                steps=output_dict.get("steps", None),
                temperature=output_dict.get("temperature", None),
                unwrapped_positions=output_dict.get("unwarapped_positions", None),
                velocities=output_dict.get("velocities", None),
                volume=output_dict["volume"]
                * ureg.angstrom
                * ureg.angstrom
                * ureg.angstrom,
                dft=OutputGenericDFT(
                    energy_free=output_dict["dft"]["energy_free"] * ureg.eV,
                    n_valence=output_dict["dft"]["n_valence"],
                    bands_k_weights=output_dict["dft"]["bands_k_weights"],
                    kpoints_cartesian=output_dict["dft"]["kpoints_cartesian"],
                    bands_e_fermi=output_dict["dft"]["bands_e_fermi"],
                    bands_occ=output_dict["dft"]["bands_occ"],
                    bands_eigen_values=output_dict["dft"]["bands_eigen_values"],
                    scf_convergence=output_dict["dft"]["scf_convergence"],
                    scf_energy_int=output_dict["dft"]["scf_energy_int"] * ureg.eV,
                    scf_energy_free=output_dict["dft"]["scf_energy_free"] * ureg.eV,
                    scf_computation_time=output_dict["dft"]["scf_computation_time"],
                    scf_energy_zero=output_dict["dft"]["scf_energy_zero"],
                    scf_energy_band=output_dict["dft"]["scf_energy_band"],
                    scf_electronic_entropy=output_dict["dft"]["scf_electronic_entropy"],
                    scf_residue=output_dict["dft"]["scf_residue"],
                    energy_int=output_dict["dft"]["energy_int"],
                    computation_time=output_dict["dft"]["computation_time"],
                    energy_zero=output_dict["dft"]["energy_zero"] * ureg.eV,
                    energy_band=output_dict["dft"]["energy_band"] * ureg.eV,
                    electronic_entropy=output_dict["dft"]["electronic_entropy"],
                    residue=output_dict["dft"]["residue"],
                ),
            ),
        ),
        job_id=job_dict["job_id"],
        status=job_dict["status"],
    )


def convert_lammps_job_dict(job_dict: dict) -> LammpsJob:
    ureg = UnitRegistry()
    generic_input_dict = convert_generic_parameters_to_dictionary(
        generic_parameter_dict=job_dict["input"]["generic"],
    )
    return LammpsJob(
        calculation_input=LammpsInput(
            generic_dict=GenericDict(
                restart_file_list=job_dict["input"]["generic_dict"][
                    "restart_file_list"
                ],
                restart_file_dict=job_dict["input"]["generic_dict"][
                    "restart_file_dict"
                ],
                exclude_nodes_hdf=job_dict["input"]["generic_dict"][
                    "exclude_nodes_hdf"
                ],
                exclude_groups_hdf=job_dict["input"]["generic_dict"][
                    "exclude_groups_hdf"
                ],
            ),
            interactive=Interactive(
                interactive_flush_frequency=job_dict["input"]["interactive"][
                    "interactive_flush_frequency"
                ],
                interactive_write_frequency=job_dict["input"]["interactive"][
                    "interactive_flush_frequency"
                ],
            ),
            generic=GenericInput(
                calc_mode=generic_input_dict["calc_mode"],
                structure=generic_input_dict["structure"],
                temperature=generic_input_dict["temperature"],
                n_ionic_steps=generic_input_dict["n_ionic_steps"],
                n_print=generic_input_dict["n_print"],
                temperature_damping_timescale=generic_input_dict[
                    "temperature_damping_timescale"
                ],
                pressure_damping_timescale=generic_input_dict[
                    "pressure_damping_timescale"
                ],
                time_step=generic_input_dict["time_step"],
                fix_symmetry=generic_input_dict.get("fix_symmetry", None),
                k_mesh_spacing=generic_input_dict.get("k_mesh_spacing", None),
                k_mesh_center_shift=generic_input_dict.get("k_mesh_center_shift", None),
                reduce_kpoint_symmetry=generic_input_dict.get(
                    "reduce_kpoint_symmetry", None
                ),
                restart_for_band_structure=generic_input_dict.get(
                    "restart_for_band_structure", None
                ),
                path_name=generic_input_dict.get("path_name", None),
                n_path=generic_input_dict.get("n_path", None),
                fix_spin_constraint=generic_input_dict.get("fix_spin_constraint", None),
                max_iter=generic_input_dict.get("max_iter", None),
            ),
            structure=Structure(
                dimension=job_dict["input"]["structure"]["dimension"],
                indices=job_dict["input"]["structure"]["indices"],
                info=job_dict["input"]["structure"]["info"],
                positions=job_dict["input"]["structure"]["positions"],
                species=job_dict["input"]["structure"]["species"],
                cell=Cell(
                    cell=job_dict["input"]["structure"]["cell"]["cell"],
                    pbc=job_dict["input"]["structure"]["cell"]["pbc"],
                ),
                units=Units(
                    length=job_dict["input"]["structure"]["units"]["length"],
                    mass=job_dict["input"]["structure"]["units"]["mass"],
                ),
            ),
            potential=LammpsPotential(
                citation=job_dict["input"]["potential_inp"]["potential"][
                    "Citations"
                ],
                config=job_dict["input"]["potential_inp"]["potential"]["Config"],
                filename=job_dict["input"]["potential_inp"]["potential"][
                    "Filename"
                ],
                model=job_dict["input"]["potential_inp"]["potential"]["Model"],
                name=job_dict["input"]["potential_inp"]["potential"]["Name"],
                species=job_dict["input"]["potential_inp"]["potential"][
                    "Species"
                ],
            ),
            input_files=LammpsInputFiles(
                control_inp=convert_generic_parameters_to_string(
                    generic_parameter_dict=job_dict["input"]["control_inp"]
                ),
                potential_inp=convert_generic_parameters_to_string(
                    generic_parameter_dict=job_dict["input"]["potential_inp"]
                ),
            ),
        ),
        executable=Executable(
            version=job_dict["executable"]["executable"]["version"],
            name=job_dict["executable"]["executable"]["name"],
            operation_system_nt=job_dict["executable"]["executable"][
                "operation_system_nt"
            ],
            executable=job_dict["executable"]["executable"]["executable"],
            mpi=job_dict["executable"]["executable"]["mpi"],
            accepted_return_codes=job_dict["executable"]["executable"][
                "accepted_return_codes"
            ],
        ),
        server=Server(
            user=job_dict["server"]["user"],
            host=job_dict["server"]["host"],
            run_mode=job_dict["server"]["run_mode"],
            queue=job_dict["server"]["queue"],
            qid=job_dict["server"]["qid"],
            cores=job_dict["server"]["cores"],
            threads=job_dict["server"]["threads"],
            new_h5=job_dict["server"]["new_h5"],
            structure_id=job_dict["server"]["structure_id"],
            run_time=job_dict["server"]["run_time"],
            memory_limit=job_dict["server"]["memory_limit"],
            accept_crash=job_dict["server"]["accept_crash"],
        ),
        calculation_output=LammpsOutput(
            generic=GenericOutput(
                cells=job_dict["output"]["generic"]["cells"] * ureg.angstrom,
                energy_pot=job_dict["output"]["generic"]["energy_pot"] * ureg.eV,
                energy_tot=job_dict["output"]["generic"]["energy_tot"] * ureg.eV,
                forces=job_dict["output"]["generic"]["forces"]
                * ureg.eV
                / ureg.angstrom,
                indices=job_dict["output"]["generic"]["indices"],
                natoms=job_dict["output"]["generic"]["natoms"],
                positions=job_dict["output"]["generic"]["positions"]
                * ureg.angstrom,
                pressures=job_dict["output"]["generic"]["pressures"] * ureg.GPa,
                steps=job_dict["output"]["generic"]["steps"],
                temperature=job_dict["output"]["generic"]["temperature"]
                * ureg.kelvin,
                unwrapped_positions=job_dict["output"]["generic"][
                    "unwrapped_positions"
                ]
                * ureg.angstrom,
                velocities=job_dict["output"]["generic"]["velocities"]
                * ureg.angstrom
                / ureg.fs,
                volume=job_dict["output"]["generic"]["volume"]
                * ureg.angstrom
                * ureg.angstrom
                * ureg.angstrom,
                dft=None,
            ),
        ),
        job_id=job_dict["job_id"],
        status=job_dict["status"],
    )
