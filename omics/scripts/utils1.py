#
import pandas as pd

def open_file(file_name, format):

    if format == "tsv":
        df = pd.read_csv(file_name, delimiter = '\t', quoting = 3, header=0)
        return df
    
    elif format == "txt":
        df = pd.read_csv(file_name, sep='\t', header=0)
        return df
    elif format == "csv":
        df = pd.read_csv(file_name, sep='\t', header=0)
        return df

