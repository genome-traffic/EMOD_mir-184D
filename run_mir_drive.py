#!/usr/bin/env python3
import copy
import pathlib  # for a join

# idmtools ...
from idmtools.builders import SweepArm, ArmType, ArmSimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import add_report_vector_genetics, add_malaria_summary_report, \
    add_report_vector_stats

from helpers import *
import params as params
import manifest as manifest

from functools import partial
from itertools import product
import pandas as pd


def get_serialization_paths(platform, serialization_exp_id):
    exp = Experiment.from_id(serialization_exp_id, children=False)
    exp.simulations = platform.get_children(exp.id, exp.item_type,
                                            children=["tags", "configuration", "files", "hpc_jobs"])

    sim_dict = {'Larval_Capacity': [], 'Outpath': []}
    for simulation in exp.simulations:
        if simulation.tags['Run_Number'] == '0':
            string = simulation.get_platform_object().hpc_jobs[0].working_directory.replace('internal.idm.ctr', 'mnt')
            string = string.replace('\\', '/')
            string = string.replace('IDM2', 'idm2')

            sim_dict['Larval_Capacity'] += [float(simulation.tags['Larval_capacity'])]
            sim_dict['Outpath'] += [string]

    df = pd.DataFrame(sim_dict)
    return df


def general_sim(serialization=0, serialized_exp_id=None):
    """
    This function is designed to be a parameterized version of the sequence of things we do
    every time we run an emod experiment.
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM")

    # create EMODTask
    print("Creating EMODTask (from files)...")

    task = EMODTask.from_default2(
        config_path="my_config.json",
        eradication_path=manifest.eradication_path,
        ep4_custom_cb=None,
        campaign_builder=None,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        demog_builder=build_demographics
    )
    task.set_sif(manifest.sif_path)

    # Create simulation sweep with builder
    builder = ArmSimulationBuilder()

    # Add asset
    task.common_assets.add_asset("~/EMOD_mir-184D/download/schema.json")

    if serialized_exp_id:
        serialized_population_path_df = get_serialization_paths(platform=platform,
                                                                serialization_exp_id=serialized_exp_id)

        # define our first Sweep Arm
        arm_baseline = SweepArm(type=ArmType.cross)
        arm_drive = SweepArm(type=ArmType.cross)

        exp_name = params.exp_name + '_high_func_resistance_high_trans'

        # Sweep larval capacit
        gambiae_lh = [9.2, 9.5, 9.8]
        species = ['gambiae']
        func = partial(update_serialize, species=species, serialization=serialization, sim_duration=6 * 365,
                       serialized_population_path_df=serialized_population_path_df)
        arm_drive.add_sweep_definition(func, gambiae_lh)
        arm_baseline.add_sweep_definition(func, gambiae_lh)

        # Sweep run number
        arm_drive.add_sweep_definition(update_sim_random_seed, range(params.nSims))
        arm_baseline.add_sweep_definition(update_sim_random_seed, range(params.nSims))

        # Create baseline
        func = partial(update_camp_type, serialize=serialization, sim_duration=6 * 365)
        arm_baseline.add_sweep_definition(func, [(True, 'baseline', 0, 0, 0, 0, 0)])

        # Create full drive release
        # Sweep mortality
        arm_drive.add_sweep_definition(sweep_mortality, [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                         1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5])

        # Sweep blood feed mortality
        bloodmeal_mortalities = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        bloodmeal_species = ['gambiae'] * len(bloodmeal_mortalities)
        arm_drive.add_sweep_definition(sweep_bloodmeal_mortality,
                                       [(bloodmeal_species[i], bloodmeal_mortalities[i])
                                        for i in range(len(bloodmeal_mortalities))])
        arm_baseline.add_sweep_definition(sweep_bloodmeal_mortality, [('gambiae', 0)])

        baseline = [False]
        start_times = [200]
        release_number = [100]
        release_frequency = [1]  # how often to release in days
        release_duration = [1]
        repetitions = [1]
        rel_params = list(product(baseline, start_times, release_number, release_frequency,
                                  release_duration, repetitions))
        func = partial(update_camp_type, serialize=serialization, sim_duration=6 * 365)
        arm_drive.add_sweep_definition(func, rel_params)
        builder.add_arm(arm_drive)
        builder.add_arm(arm_baseline)

        add_report_vector_genetics(task, manifest, species='gambiae', stratify_by='ALLELE_FREQ',
                                   include_vector_state=False,
                                   include_death_state=False, combine_similar_genomes=True)

        # Add vector stats report
        add_report_vector_stats(task, manifest, include_death_state=True)

        # Add malaria summary report
        add_malaria_summary_report(task, manifest, age_bins=[5, 15, 125], reporting_interval=365,
                                   filename_suffix='Annual', start_day=200, max_number_reports=100)


    else:
        arm = SweepArm(type=ArmType.cross)

        gambiae = [9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8]
        species = ['gambiae']
        func = partial(update_serialize, species=species, serialization=serialization, sim_duration=40 * 365,
                       serialized_population_path_df=None)
        arm.add_sweep_definition(func, gambiae)

        arm.add_sweep_definition(update_sim_random_seed, range(3))

        func = partial(update_camp_type, serialize=serialization, sim_duration=40 * 365)
        arm.add_sweep_definition(func, [(1, 0, 0, 0, 0, 0)])

        builder.add_arm(arm)

        exp_name = params.exp_name + '_serialization'
        # exp_name = 'vector_aging_serialization'

    # create experiment from builder
    print(f"Prompting for COMPS creds if necessary...")
    experiment = Experiment.from_builder(builder, task, name=exp_name)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")

    # Save experiment id to file
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)


if __name__ == "__main__":

    serialization = 0
    serialization_experiment_id = '98d7be21-b6b5-ef11-aa1a-b88303911bc1'
    # serialization = 1
    # serialization_experiment_id = None
    general_sim(serialization=serialization, serialized_exp_id=serialization_experiment_id)
