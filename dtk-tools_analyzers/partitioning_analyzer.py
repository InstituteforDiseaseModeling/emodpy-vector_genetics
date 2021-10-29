import os
import sys

sys.path.append(os.path.dirname(__file__))
import itertools
import time
import pandas as pd
#import humanfriendly as humanfriendly
import numpy as np
import psutil
from math import floor
from concurrent.futures import as_completed, ProcessPoolExecutor
from functools import partial
from pandas.core.dtypes.dtypes import CategoricalDtype
from tqdm import tqdm
from simtools.Utilities.COMPSCache import COMPSCache
from simtools.Utilities.COMPSUtilities import COMPS_login, get_asset_files_for_simulation_id
from simtools.Utilities.RetryDecorator import retry
from collections import defaultdict
from simtools.Analysis.DataRetrievalProcess import retrieve_COMPS_AM_files, retrieve_SSMT_files, set_exception
from simtools.Utilities import pluralize, on_off
from simtools.Analysis.OutputParser import SimulationOutputParser
from simtools.Analysis.AnalyzeManager import AnalyzeManager


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(itertools.islice(iterable, n))


def retrieve_data_for_simulation(simulation, analyzers, path_mapping):
    # Filter first and get the filenames from filtered analysis
    filtered_analysis = [a for a in analyzers if a.filter(simulation)]
    filenames = set(itertools.chain(*(a.filenames for a in filtered_analysis)))

    # We dont have anything to do :)
    if not filtered_analysis:
        return

    # The byte_arrays will associate filename with content
    byte_arrays = {}

    # Retrieval for SSMT
    if path_mapping:
        byte_arrays = retrieve_SSMT_files(simulation, filenames, path_mapping)

    # Retrieval for normal HPC Asset Management
    elif simulation.experiment.location == "HPC":
        byte_arrays = retrieve_COMPS_AM_files(simulation, filenames)

    # Retrieval for local file
    else:
        for filename in filenames:
            path = os.path.join(simulation.get_path(), filename)
            with open(path, 'rb') as output_file:
                byte_arrays[filename] = output_file.read()

    # Selected data will be a dict with analyzer.uid => data
    selected_data = {}
    for analyzer in filtered_analysis:
        # If the analyzer needs the parsed data, parse
        if analyzer.parse:
            data = {filename: SimulationOutputParser.parse(filename, content)
                    for filename, content in byte_arrays.items()}
        else:
            # If the analyzer doesnt wish to parse, give the raw data
            data = byte_arrays

        # Retrieve the selected data for the given analyzer
        selected_data[analyzer.uid] = analyzer.select_simulation_data(data, simulation)
    return selected_data


class PartitioningDataAnalyzeManager(AnalyzeManager):
    def __init__(self, partitionable_columns, optimize_types=False, discover_partitions=False, ram_limit=2**30*8,
                 **kwargs):
        """
        :param partitionable_columns: List of columns that we can partition by
        :param optimize_types:  Optimize columns to categorical types which are more efficient in memory and matrix
            operations
        :param discover_partitions: Should we evaulute how many sims are in each partition. useful for more accurate
            progress but can be slow to load.
        :param ram_limit: str on how much ram to use(default 8gb limit). Effects the required amount of partitions needed
        :param kwargs:
        """
        super().__init__(**kwargs)
        self.partionable_columns = partitionable_columns
        self.optimize_types = optimize_types
        self.discover_partitions = discover_partitions
        if ram_limit == 'auto':
            self.ram_limit = floor(psutil.virtual_memory() / 1.5)
        else:
            self.ram_limit = ram_limit
            #self.ram_limit = humanfriendly.parse_size(ram_limit)

    def analyze_partition(self, pool, ssmt_path_mapping, max_threads, simulations, tag_types, pbar,
                          partition_vars, partition_vars_vals):
        final = defaultdict(list)
        get_data = partial(retrieve_data_for_simulation, analyzers=self.analyzers, path_mapping=ssmt_path_mapping)
        futures = {pool.submit(get_data, s): s for s in simulations}
        for future in as_completed(list(futures.keys())):
            if future.result() is not None:  # Since retrieve_data_for_simulation could filter and return none
                for k, v in future.result().items():
                    if self.optimize_types:
                        # optimize the types
                        for tag, ty in tag_types.items():
                            v[tag] = v[tag].astype(ty)
                        final[k].append(v)
                    else:
                        final[k].append(v)
            pbar.update()
            # delete reference to free memory
            # del futures[future]
        del futures
        final = {k: pd.concat(v) for k, v in final.items()}
        for a in tqdm(self.analyzers, desc="Writing partition analysis"):
            # - NEW: could update self.partition_vars_vals, self.partition_vars here instead
            final[a.uid] = a.finalize(final[a.uid], partition_vars, partition_vars_vals)

    def analyze(self):
        start_time = time.time()
        # Clear the cache
        self.cache = self.initialize_cache(shards=self.max_threads)

        # Check if we are on SSMT
        ssmt_path_mapping = os.environ.get("COMPS_DATA_MAPPING", None)
        if ssmt_path_mapping:
            # - for analyzing runs from Calculon
            ssmt_path_mapping = ssmt_path_mapping.split(';')
            # - for analyzing runs from Belegost
            # ssmt_path_mapping = ssmt_path_mapping.lower().split(';')

        # If any of the analyzer needs the dir map, create it
        # Or if we are on SSMT
        if ssmt_path_mapping or any(a.need_dir_map for a in self.analyzers):
            # preload the global dir map
            from simtools.Utilities.SimulationDirectoryMap import SimulationDirectoryMap
            for experiment in self.experiments:
                SimulationDirectoryMap.preload_experiment(experiment)

        # Run the per experiment on the analyzers
        # for exp in self.experiments:
        #     for a in self.analyzers:
        #         a.per_experiment(exp)

        # Display some info

        # - TEST 2
        # subset = take(100, self.simulations.keys())
        # self.simulations = {k:v for k, v in self.simulations.items() if k in subset}

        scount = len(self.simulations)

        # determine how to partition data
        tags = defaultdict(set)
        for sim_id, siminfo in self.simulations.items():
            for tag_name, value in siminfo.tags.items():
                if tag_name in self.partionable_columns:
                    tags[tag_name].add(value)

        print('Possible Partitions: ' + str(tags))
        tag_types = dict()
        for tag, values in tags.items():
            tag_types[tag] = CategoricalDtype(categories=list(sorted(values)), ordered=True)
        # now sort the tags
        tag_lookup = tags
        tags = sorted(tags.items(), key=lambda x: len(x[1]))
        possible_partitions = np.prod([len(tag[1]) for tag in tags])
        print(f"Possible Partitions based on colums: {possible_partitions}")
        max_threads = min(self.max_threads, scount if scount != 0 else 1)

        if self.verbose:
            print("Analyze Manager")
            print(" | {} simulation{} - {} experiment{}"
                  .format(scount, pluralize(scount), len(self.experiments), pluralize(self.experiments)))
            print(" | force_analyze is {} and {} simulation{} ignored"
                  .format(on_off(self.force_analyze), len(self.ignored_simulations),
                          pluralize(self.ignored_simulations)))
            print(" | Analyzer{}: ".format(pluralize(self.analyzers)))
            for a in self.analyzers:
                print(" |  - {} (Directory map: {} / File parsing: {} / Use cache: {})"
                      .format(a.uid, on_off(a.need_dir_map), on_off(a.parse), on_off(hasattr(a, "cache"))))
            print(" | Pool of {} analyzing processes".format(max_threads))

        if scount == 0 and self.verbose:
            print("No experiments/simulations for analysis.")
            return False

        # determine memory required for each analyzer and what columns are available in the result
        sample_data = retrieve_data_for_simulation(next(iter(self.simulations.values())), self.analyzers, ssmt_path_mapping)
        sample_data = {k: df.memory_usage(index=True).sum() for k, df in sample_data.items()}
        mem_per_sim = sum(sample_data.values())
        total_mem = mem_per_sim * len(self.simulations)
        print(f'Memory needed for analysis required: {round(total_mem/2**30, 2)} GB')
        # print(f'Memory needed for analysis required: {humanfriendly.format_size(total_mem)}')
        partitions_required = total_mem / self.ram_limit * 4

        # - TEST 4
        # partitions_required = 1.1

        print(f'Partition Required: {partitions_required}')

        # - TEST 1
        # exit()

        partitions = tags
        sets = self.get_sim_sets(partitions, partitions_required, tags)

        # - TEST 3
        # print(sets)
        # exit()

        total_sims = 0
        pbar = tqdm(desc="Analyzing partitions", total=len(self.simulations), smoothing=0.01)
        pool = ProcessPoolExecutor(max_workers=max_threads)
        for sims in sets:

            # - NEW
            partition_vars = list(sims.keys())
            partition_vars_vals = list(sims.values())

            # assume it is still generated combo group so we need to filter sims first
            if type(sims) is dict:
                sims = list(filter(lambda x: all([p in x.tags and x.tags[p] == sims[p] for p in sims.keys()]),
                                   self.simulations.values()))

            if len(sims) > 0:
                total_sims += len(sims)
                self.analyze_partition(pool, ssmt_path_mapping, max_threads, sims, tag_types, pbar,
                                       partition_vars, partition_vars_vals)  # - NEW
        print(f"Process {total_sims}")

    def get_sim_sets(self, partitions, partitions_required, tags):
        last_i = None
        for i in range(len(tags)):
            if np.prod([len(tag[1]) for tag in tags][:i]) >= partitions_required:
                last_i = i
                partitions = tags[:last_i]
                break
        partitions = [tag[0] for tag in partitions]
        vlist = [tags[i][1] for i in range(len(partitions))]
        sets = self.get_sets_from_partitions(partitions, vlist, self.discover_partitions)
        if len(sets) < partitions_required:
            if last_i < len(tags):
                print(f'Only found {len(sets)} using {str(partitions)}. Trying more fields for partitioning')
                partitions = tags[:last_i + 1]
                partitions = [tag[0] for tag in partitions]
                vlist = [tags[i][1] for i in range(len(partitions))]
                sets = self.get_sets_from_partitions(partitions, vlist, self.discover_partitions)
            else:  # not enough partitions to segment data
                return sets
        print(f"Partitioned by {str(partitions)}")
        return sets

    def get_sets_from_partitions(self, partitions, vlist, filter_sims=True):
        combos = [dict(zip(partitions, v)) for v in itertools.product(*vlist)]
        if filter_sims:
            sets = []
            for partition in tqdm(combos, desc="Gathering partition info"):
                # filter sims for ones that match this partition
                sim_set = list(filter(lambda x: all([p in x.tags and x.tags[p] == partition[p] for p in partitions]),
                                      self.simulations.values()))
                if sim_set:
                    sets.append(sim_set)
            print(f'Found {len(sets)} partitions in the data')
            return sets
        return combos  # return combos if we don't filter