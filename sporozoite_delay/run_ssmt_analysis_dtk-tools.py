from simtools.Managers.WorkItemManager import WorkItemManager
from simtools.SetupParser import SetupParser
from simtools.AssetManager.FileList import FileList

# ------------------------------------
# - Use w/ user_files = FileList(root='analyzers')

# wi_name = 'spatialinside_integral2l4a_GM_only_aEIR10_sweep_rc_d1_rr20_se2'
# wi_name = "spatial_vector_genetics_classicdrive3allele_VC_and_GM_sweep_sd_sne_genomes_4nodesonly"
# command = "python run_analysis_nash_reproduction_spatial.py"

# wi_name = "single_node_whole_construct_classic_drive_1_locus_sweep_hrc_rc"
# command = "python run_analysis_nash_reproduction_single_node.py"

# ------------------------------------
# - Use w/ user_files = FileList(root='analyzers_large_data')
# and docker_image = "docker-production.packages.idmod.org/dtk-tools/dtk-tools-feather-tqdm-hdf5:1.0.1"

wi_name = "sporozoite_delay_VG"

command = "python3 run_large_data_analysis.py"

# ------------------------------------
# user_files = FileList(root='analyzers')
user_files = FileList(root='C:\\Users\\pselvaraj\\Github\\emodpy-vector_genetics\\dtk-tools_analyzers')

if __name__ == "__main__":
    SetupParser.default_block = "HPC"
    SetupParser.init()

    wim = WorkItemManager(item_name=wi_name, command=command, user_files=user_files,
                          docker_image="docker-production.packages.idmod.org/dtk-tools/dtk-tools-feather-tqdm-hdf5:1.0.1",  # use w/ analyzers_large_data
                          related_experiments=["8bbe4c48-6326-ec11-9ecd-9440c9bee941"])
    wim.comps_env = 'Calculon'
    wim.execute(True)