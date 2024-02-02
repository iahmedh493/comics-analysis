import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from sanbomics.tools import id_map
import pandas as pd

from .load_data import assemble
from .utils1 import open_file

class ExplorAnalysis:

    def __init__(self, param = []) -> None:
        self.param = param

        if len(self.param[2]) == 1:
            self.df = open_file(self.param[2][0], self.param[1])
        elif len(self.param[2]) >1 :
            path_ = os.path.split(self.param[2][0])[0]
            files = []
            for item in self.param[2]:
                files.append(os.path.split(item)[1])
            self.df = assemble(path_, files, self.param[1])
        self.df.index.name = "GeneId"
        
        self.meta = self.meta_data()

    
    def filter_zero(self):
        self.df = self.df[self.df.sum(axis = 1) > 0]
        

    def meta_data(self):
        self.df = self.df.T
        meta = pd.DataFrame(zip(self.df.index, self.param[3]), columns=["Sample", "Condition"])
        meta = meta.set_index('Sample')
        return meta
    
    def deseq_obj(self):
        dds = DeseqDataSet(
            counts=self.df,
            metadata=self.meta,
            design_factors="Condition"
            )
        dds.deseq2()
        stats_res = DeseqStats(dds, contrast = ('Condition', 'LP', 'ML'))
        stats_res.summary()
        res = stats_res.results_df
        return res , dds, stats_res
    
    #def id_mapping(self):
       # res, dds = self.deseq_obj()
        #mapper = id_map(species = 'human')
        #res['Symbol'] = res.index.map(mapper.mapper)
        

