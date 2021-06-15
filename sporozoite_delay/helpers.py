import sporozoite_delay.set_config as set_config
from emodpy_malaria import config as malconf
from emodpy_malaria import vector_config as vecconf
from functools import \
    partial  # for setting Run_Number. In Jonathan Future World, Run_Number is set by dtk_pre_proc based on generic param_sweep_value...
import sporozoite_delay.manifest as manifest

import emod_api.config.default_from_schema_no_validation as dfs
import emod_api.campaign as camp
import emodpy_malaria.interventions.bednet as bednet
from emodpy_malaria.interventions.treatment_seeking import add_healthseeking
from emodpy_malaria.interventions.mosquitorelease import mosquitorelease


def update_sim_bic(simulation, value):
    simulation.task.config.parameters.Base_Infectivity_Constant = value * 0.1
    return {"Base_Infectivity": value}


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def update_camp_start_day(simulation, value):
    # simulation.task.config.parameters.Run_Number = value
    build_camp_partial = partial(build_camp, start_day_in=80 + value * 10)
    simulation.task.create_campaign_from_callback(build_camp_partial)
    return {"Start_Day": 80 + value * 10}


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.
    """
    config = set_config.set_config(config)

    config.parameters.Base_Rainfall = 150
    config.parameters.Simulation_Duration = 365
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    # config.parameters.Serialization_Times = [ 365 ]
    config.parameters.Enable_Vector_Species_Report = 1
    # config["parameters"]["Insecticides"] = [] # emod_api gives a dict right now.
    config.parameters.pop("Serialized_Population_Filenames")

    # # Create new MalariaDrugParams
    # config = set_mdp(config, manifest)

    # Gene drive - using full dictionary here but hope to change to something more elegant in the future
    drivers = {"gambiae": [
        {
            "Alleles_Driven": [
                {
                    "Allele_To_Copy": "a1",
                    "Allele_To_Replace": "a0",
                    "Copy_To_Likelihood": {
                        "a0": 0.1,
                        "a1": 0.9
                    }
                },
                {
                    "Allele_To_Copy": "b1",
                    "Allele_To_Replace": "b0",
                    "Copy_To_Likelihood": {
                        "b0": 0.1,
                        "b1": 0.9
                    }
                }
            ],
            "Driver_Type": "INTEGRAL_AUTONOMOUS",
            "Driving_Allele": "a1"
        },
    ]}

    # Vector Genetics
    malconf.set_species(config, ['gambiae'])

    vecconf.add_alleles(['a0', 'a1', 'a2'], [1.0, 0.0, 0.0])
    vecconf.add_alleles(['b0', 'b1', 'b2'], [1.0, 0.0, 0.0])

    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b1"]], "INFECTED_PROGRESS", 0.5)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b0"]], "INFECTED_PROGRESS", 0.75)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b2"]], "INFECTED_PROGRESS", 0.75)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b1"]], "TRANSMISSION_TO_HUMAN", 0.6)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b0"]], "TRANSMISSION_TO_HUMAN", 0.8)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b2"]], "TRANSMISSION_TO_HUMAN", 0.8)

    vecconf.set_species_drivers(config, drivers)

    # # Crete new Vector Species Params
    # config = set_vsp(config, manifest)
    return config


def build_camp(start_day_in=100):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    camp.schema_path = manifest.schema_file

    # print( f"Telling emod-api to use {manifest.schema_file} as schema." )
    camp.add(bednet.Bednet(camp, start_day=start_day_in, coverage=0.5, killing_eff=0.5, blocking_eff=0.5, usage_eff=0.5,
                           insecticide="pyrethroid"), first=True)
    camp.add(add_healthseeking(camp, targets=[{ "trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                                "agemax": 70, "seek": 0.4, "rate": 0.3},
                                              {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                               drug=['Artemether', 'Lumefantrine'],
                               start_day=0,
                               broadcast_event_name='Received_Treatment'
                               ))
    camp.add(add_healthseeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                               "agemax": 70, "seek": 0.4, "rate": 0.3},
                                              {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                               drug=['Artemether', 'Lumefantrine'],
                               start_day=0,
                               broadcast_event_name='Received_Treatment'
                               ))
    return camp


def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    Also right now this function takes care of the config updates that are required as a result of specific demog settings. We do NOT want the emodpy-disease developers to have to know that. It needs to be done automatically in emod-api as much as possible.
    TBD: Pass the config (or a 'pointer' thereto) to the demog functions or to the demog class/module.

    """
    import emodpy_malaria.demographics.MalariaDemographics as Demographics  # OK to call into emod-api
    import emod_api.demographics.DemographicsTemplates as DT

    demog = Demographics.from_template_node(lat=0, lon=0, pop=10000, name=1, forced_id=1)
    return demog


def rvg_config_builder(params):
    params.Include_Vector_State_Columns = False
    params.Allele_Combinations_For_Stratification = [
        ["tom"],
        ["dick"]
    ]
    # params.Alleles_For_Stratification = [ "tom", "dick" ] # this works
    # params.Specific_Genome_Combinations_For_Stratification =
    """
    E.g.,
    [
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a1", "*"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a0"]
            ]
        }
    ]
    """
    return params
