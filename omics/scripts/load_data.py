
# Data pre-process
# RNASeq workflow
import numpy as np
import pandas as pd
from rnanorm import CPM
import csv
import os
import pathlib
from .utils1 import open_file


#workDir = "/Users/ibrahimahmed/October23/Data"
'''
files = ("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
   "GSM1545545_JMS9-P8c.txt") 

type = ""
'''

def samples_dict(dir_, files_, type):
    main_dict = {}
    for file_name in files_:
        path_ = os.path.join(dir_, file_name)
        temp, sample_name = read_data(path_, file_name, type)
        main_dict[sample_name] = temp

    return main_dict

def read_data(file_path, file_, type):

    sample_name = file_[11:-4]
    #sample_name = file_
    temp = {}
    #dataf = pd.read_csv(file_path, sep="\t")
    dataf = open_file(file_path, type)
    dataf1 = dataf[["EntrezID", "Count"]]
    temp = {}
    for index, row in dataf1.iterrows():
        temp[row['EntrezID']] = row['Count']
    return temp, sample_name


def unique_ids(main_dict):
    ids_list = []
    for key, values in main_dict.items():
        for item in values.keys():
            if item not in ids_list:
                ids_list.append(item)
    return ids_list


def create_matrix(main_dict, ids_data):

    samples = main_dict.keys()
    values_ = list(main_dict.values())
    matrix_dict = {}
    for item in ids_data:
        temp = []
        for list_ in values_:
            if item in list_.keys():
                temp.append(list_[item])
            else:
                temp.append(0)
        matrix_dict[item] = temp
    df2 = pd.DataFrame.from_dict(matrix_dict, orient='index', columns=samples)
    return df2


def assemble(work_dir_, files_, type):
    dict_main = samples_dict(work_dir_, files_, type)
    ids = unique_ids(dict_main)
    df = create_matrix(dict_main, ids)
    # df_normalized = CPM().set_output(transform="pandas").fit_transform(df1)
    ##print(df.head(4))
    # df.set_index('EntrezID')
    # df['EntrezID'] = df.index
    dir_path = pathlib.Path().resolve()
    file_name = "temporary.csv"
    full_file_path = os.path.join(dir_path, file_name)
    print(full_file_path)
    return df
    # df.to_csv(file_name, index=True)

# print(df.head(4))
# df1 = df.T


#df = assemble(workDir, files, "txt")
### print(list(df_normalized.index))
#print(df.head(4))

#print(df.columns)
#print(df.head().index)

