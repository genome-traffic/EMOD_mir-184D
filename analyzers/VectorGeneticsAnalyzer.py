import json
import os
import pandas as pd
import numpy as np
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import matplotlib as mpl
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation

from logging import getLogger

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.analysis.csv_analyzer import CSVAnalyzer
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

mpl.use('Agg')


class VectorGeneticsAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm'):
        super().__init__(filenames=[
                                    "output\\ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv"])
        self.tags = ['Baseline', 'Run_Number', 'Mortality', 'Bloodmeal_Mortality', 'Larval_capacity']
        # self.tags = ['Baseline', 'Run_Number', 'Bloodmeal_Mortality', 'Fecundity', 'Mortality']
        print(title)

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:
        """
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        """
        Extracts the Statistical Population, Data channel from InsetChart.
        Called for Each WorkItem/Simulation.
        Args:
            data: Data mapping str to content of file
            item: Item to Extract Data from(Usually a Simulation)
        Returns:
        """
        df = data[self.filenames[0]]

        df = df[df['Alleles'].isin(['a0', 'a1', 'a2', 'a3'])]
        df = df[['Time', 'Alleles', 'VectorPopulation']]
        df["VectorPopulation"] = df["VectorPopulation"] / df.groupby(["Time"])['VectorPopulation'].transform("sum")

        pivoted_df = df.pivot(index="Time", columns="Alleles", values="VectorPopulation").reset_index()

        return pivoted_df

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot
        Args:
            all_data: Populate data_simulation from all the Simulations
        Returns:
            None
        """
        output_dir = os.path.join(self.working_dir, "output")
        print('Starting append\n')

        df_final = pd.DataFrame()
        for s, v in all_data.items():
            dftemp = v.copy()
            for t in self.tags:
                dftemp[t] = [s.tags[t]]*len(v)
            dftemp.set_index(self.tags)
            df_final = pd.concat([df_final, dftemp])
        df_final.to_pickle(os.path.join(output_dir, 'MIR_drive_allele_freq_full.csv'))

        groupby_tags = self.tags
        groupby_tags.remove('Run_Number')

        df_allele_final = df_final.groupby(groupby_tags + ['Time'])[['a0', 'a1', 'a2', 'a3']].mean().reset_index()

        df_allele_final_std = df_final.groupby(groupby_tags + ['Time'])[['a0', 'a1', 'a2', 'a3']].apply(np.std)
        for c in ['a0', 'a1', 'a2', 'a3']:
            df_allele_final[c + '_std'] = list(df_allele_final_std[c])

        df_allele_final.to_csv(os.path.join(output_dir, "MIR_drive_allele_freq.csv"))


if __name__ == '__main__':

    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    logger = getLogger()
    with Platform('CALCULON') as platform:

        # Initialize the analyser class with the path of the output csv file
        analyzers = [VectorGeneticsAnalyzer()]

        # Set the experiment id you want to analyze
        experiment_id = 'f36d11e6-8cdf-eb11-a9ec-b88303911bc1'

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()
