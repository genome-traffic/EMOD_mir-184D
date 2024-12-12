import os
import pandas as pd
import numpy as np
from functools import \
    partial  # for setting Run_Number. In Jonathan Future World, Run_Number is set by dtk_pre_proc based on generic param_sweep_value...
import manifest as manifest

import emod_api.demographics.Demographics as Demographics

from emodpy_malaria import malaria_config as malconf
from emodpy_malaria import vector_config as vecconf
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.mosquitorelease import add_scheduled_mosquito_release

import emod_api.campaign as camp


def build_demographics():

    demog = Demographics.from_file(
        "~/EMOD_mir184D/input_files/single_node_demographics.json")

    return demog


def find_genome_index_in_trait_modifiers(vsp, genome, trait_of_interest):
    flat_gen = [i for j in genome for i in j]
    flat_gen.sort()

    trait_modifiers = vsp['Gene_To_Trait_Modifiers']

    for i, t in enumerate(trait_modifiers):

        config_gen = [m for j in t['Allele_Combinations'] for m in j]
        config_gen.sort()
        traits = [j['Trait'] for j in t['Trait_Modifiers']]
        t_index = None
        if trait_of_interest in traits:
            t_index = traits.index(trait_of_interest)

        if (t_index is not None) and (config_gen == flat_gen):
            return i, t_index

    return None


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def sweep_sim_larval_capacity(simulation, species, value):

    for j, s in enumerate(species):
        for i, vsp in enumerate(simulation.task.config.parameters.Vector_Species_Params):
            if vsp['Name'] == s:
                simulation.task.config.parameters.Vector_Species_Params[i].Habitats[0]['Max_Larval_Capacity'] \
                    = pow(10, value)

    return None


def sweep_mortality(simulation, value):
    vsp = simulation.task.config.parameters.Vector_Species_Params[0]

    trait = "MORTALITY"

    # Both males and females
    genome = [["a1", "a1"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    genome = [["a3", "a3"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    genome = [["a1", "a3"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = value

    genome = [["a0", "a1"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 + (value - 1) / 2

    genome = [["a0", "a3"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 + (value - 1) / 2

    genome = [["a1", "a2"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 + (value - 1) / 2

    genome = [["a2", "a3"]]
    idx1, idx2 = find_genome_index_in_trait_modifiers(vsp=vsp, genome=genome, trait_of_interest=trait)
    vsp.Gene_To_Trait_Modifiers[idx1].Trait_Modifiers[idx2].Modifier = 1 + (value - 1) / 2

    return {"Mortality": value}


def sweep_bloodmeal_mortality(simulation, value):
    species = value[0]
    bloodmeal_mortality = value[1]

    for i, vsp in enumerate(simulation.task.config.parameters.Vector_Species_Params):
        if vsp['Name'] == species:
            bmm_gp = simulation.task.config.parameters.Vector_Species_Params[
                i].Blood_Meal_Mortality.Genetic_Probabilities
            for j, gp in enumerate(bmm_gp):
                bmm_gp[j].Probability = bloodmeal_mortality

    return {"Bloodmeal_Mortality": value[1]}


def update_serialize(simulation, larval_multiplier, species, serialization=0, sim_duration=40 * 365,
                     serialized_population_path_df=None):
    if serialization:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Time_Steps = [sim_duration]
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)

    else:
        serialized_population_path = serialized_population_path_df[
            (serialized_population_path_df['Larval_Capacity'] == larval_multiplier)]['Outpath'].values[0]
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Mask_Node_Read = 0
        simulation.task.config.parameters.Serialization_Mask_Node_Write = 0
        simulation.task.config.parameters.Serialized_Population_Path = os.path.join(serialized_population_path,
                                                                                    'output')
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'READ'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Filenames = ['state-14600.dtk']
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)

    return {"Serialization": serialization, 'Larval_capacity': larval_multiplier}


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.
    """
    config = malconf.set_team_defaults(config, manifest)

    config.parameters.Age_Initialization_Distribution_Type = "DISTRIBUTION_OFF"
    config.parameters.Base_Rainfall = 150
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    config.parameters.Enable_Demographics_Risk = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Enable_Vector_Aging = 1
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Demographics_Filenames = ['demographics.json']
    config.parameters.Enable_Demographics_Birth = 1

    config.parameters.Enable_Spatial_Output = 1
    config.parameters.Spatial_Output_Channels = ['Prevalence', 'Population']

    config.parameters.Serialization_Mask_Node_Read = 0

    # Vector species params
    df = get_ento_values()

    # gambiae
    malconf.add_species(config, manifest, ['gambiae'])
    values = list(df['gambiae'])
    habitats = vecconf.configure_linear_spline(manifest, max_larval_capacity=pow(10, 8),
                                               capacity_distribution_number_of_years=1,
                                               capacity_distribution_over_time={
                                                   "Times": [0.0, 30.417, 60.833, 91.25,
                                                             121.667, 152.083,
                                                             182.5, 212.917, 243.333,
                                                             273.75, 304.167, 334.583],
                                                   "Values": values})
                                                    # "Values": [1]*12})
    vecconf.set_species_param(config, 'gambiae', 'Habitats', [habitats], overwrite=True)
    vecconf.set_species_param(config, 'gambiae', 'Adult_Life_Expectancy', 20, overwrite=True)
    vecconf.set_species_param(config, 'gambiae', 'Anthropophily', 1, overwrite=True)

    # Vector Genetics
    # --------- Define classic gene drive functions
    # - Initial resistance frequency (defined in set_classic_genes)
    rr0 = 0  # initial frequency of functional resistance
    # - Gene drive params (defined in add_classic_gene_drives)
    cut_rate = 1 # overall cut_rate (in most cases 100%)
    d = 0.99  # transmission rate of drive (cut rate * HDR rate)
    u12 = 0.1 # Functional:nonFunctional NHEJ repair ratio (most repair events are non-functional)

    # - Fitness params (defined in add_classic_fitness_costs)

    sd = 1  # cost of target site disruption
    hd = 0.5  # dominance coefficient for target site disruption

    # DEFINE ALLELES
    vecconf.add_genes_and_alleles(config, manifest, 'gambiae', [('a0', 1.0 - rr0), ('a1', 0), ('a2', rr0), ('a3', 0)])

    # ADD DRIVERS
    a1_ctl = [("a0", (1 - cut_rate)),  # Wild type
              ("a1", cut_rate * d),  # Driver
              ("a2", cut_rate * (1 - d) * u12),  # Functional Resistance
              ("a3", cut_rate * (1 - d) * (1 - u12))]  # nonFunctional Resistance

    vecconf.add_species_drivers(config, manifest, species='gambiae', driving_allele='a1',
                                driver_type="CLASSIC", to_copy='a1', to_replace='a0',
                                likelihood_list=a1_ctl)

    # MORTALITY
    # drive/LOF heterozygous
    mortality = vecconf.create_trait(manifest, trait='MORTALITY', modifier=1 + hd * sd)
    vecconf.add_trait(config, manifest, 'gambiae', [["a0", "a1"]],
                      [mortality])
    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a2"]],
                      [mortality])
    vecconf.add_trait(config, manifest, 'gambiae', [["a0", "a3"]],
                      [mortality])
    vecconf.add_trait(config, manifest, 'gambiae', [["a2", "a3"]],
                      [mortality])

    # drive/LOF homozygous
    mortality = vecconf.create_trait(manifest, trait='MORTALITY', modifier=1 + sd)
    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a1"]],
                      [mortality])

    vecconf.add_trait(config, manifest, 'gambiae', [["a1", "a3"]],
                      [mortality])

    vecconf.add_trait(config, manifest, 'gambiae', [["a3", "a3"]],
                      [mortality])

    # BLOODMEAL MORTALITY
    vecconf.add_blood_meal_mortality(config, manifest,
                                     default_probability_of_death=0,
                                     species="gambiae",
                                     allele_combo=[["a1", "a1"]],
                                     probability_of_death_for_allele_combo=0.5)

    vecconf.add_blood_meal_mortality(config, manifest,
                                     default_probability_of_death=0,
                                     species="gambiae",
                                     allele_combo=[["a1", "a3"]],
                                     probability_of_death_for_allele_combo=0.5)

    vecconf.add_blood_meal_mortality(config, manifest,
                                     default_probability_of_death=0,
                                     species="gambiae",
                                     allele_combo=[["a3", "a3"]],
                                     probability_of_death_for_allele_combo=0.5)

    return config


def update_camp_type(simulation, release_params, serialize=0, sim_duration=40 * 365):
    # simulation.task.config.parameters.Run_Number = value
    build_camp_partial = partial(build_camp, serialize=serialize,
                                 sim_duration=sim_duration, release_params=release_params)
    simulation.task.create_campaign_from_callback(build_camp_partial)

    return {"Baseline": release_params[0], "Start_time": release_params[1],
            "Release_number": release_params[2],
            "Release_frequency": release_params[3],
            "Release_duration": release_params[4],
            "Repetitions": release_params[5]}


def build_camp(release_params, serialize=0, sim_duration=40 * 365):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    camp.set_schema(manifest.schema_file)

    if not serialize:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                              "agemax": 70, "rate": 0.4},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.9}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=0,
                              broadcast_event_name='Received_Treatment'
                              )

        if not release_params[0]:

            rounds_per_year = 1
            if release_params[4] > 1:
                rounds_per_year = np.round(release_params[5]/release_params[4])

            if release_params[2] > 0:

                # Number of years to do this for
                for i in range(release_params[5]):

                    add_scheduled_mosquito_release(camp, start_day=i*365 + release_params[1],
                                                   released_number=release_params[2],
                                                   repetitions=rounds_per_year,
                                                   timesteps_between_repetitions=release_params[3],
                                                   released_species='gambiae',
                                                   released_genome=[["X", "Y"], ["a1", "a1"]])

    else:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                              "agemax": 70, "rate": 0.4},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.9}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=sim_duration - 10 * 365,
                              broadcast_event_name='Received_Treatment'
                              )
    return camp


def get_ento_values():
    file = '~/Input to mosquito data file'

    df = pd.read_csv(file)
    df = df[['Month', 'Fune_collected', 'Gamb_collected']]
    df.rename(columns={'Fune_collected': 'funestus', 'Gamb_collected': 'gambiae'}, inplace=True)
    df['funestus'] = df['funestus'] / sum(df['funestus'])
    df['gambiae'] = df['gambiae'] / sum(df['gambiae'])

    return df