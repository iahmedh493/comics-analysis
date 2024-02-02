
# this module perform data exploration

from typing import Union, List, Any, Literal, Iterable, Tuple
from pathlib import Path
import pandas as pd
import numpy as np

def isiterable(obj):
    """
    Returns True if obj is Iterable (list, str, tuple, set, dict, etc'), and False otherwise. \
    This function does not consume iterators/generators.

    :param obj: the object to be tested for iterability
    :type obj: Any
    :return: True if obj is Iterable, False otherwise.
    :rtype: bool
    """
    try:
        _ = iter(obj)
    except TypeError:
        return False
    else:
        return True

def isinstanceiter(iterable: Iterable, object_class: Union[type, Tuple[type, ...]]):
    """
    Returns True if all members of an Iterable object are instances of a class or of a subclass thereof. \
    This function consumes iterators/generators. Always returns True when 'iterable' is empty.

    :param iterable: the Iterable object whose members' types should be checked.
    :type iterable: Iterable (list, tuple, str, dict, set, etc')
    :param object_class: the class/classes to check 'isinstance' against
    :type object_class: type (e.g. list, tuple, int, bool) or tuple of types
    :return: True if all members of 'iterable' are of type 'object_class', and False otherwise.
    :rtype: bool
    """
    assert isiterable(iterable), f"Object of type {type(iterable)} is not iterable."
    return all([isinstance(i, object_class) for i in iterable])

def data_to_list(data: Any, sort: bool = False) -> list:
    if isinstance(data, list):
        lst = data
    elif isinstance(data, (set, tuple, np.ndarray)):
        lst = list(data)
    elif isinstance(data, (dict, int, float, bool, str, pd.DataFrame, pd.Series)):
        lst = [data]
    elif data is None:
        lst = [None]

    elif callable(data):
        lst = [data]
    else:
        try:
            lst = list(data)
        except TypeError:
            lst = [data]

    if sort:
        lst.sort()
    return lst

def load_table(filename: Union[str, Path], index_col: int = None, drop_columns: Union[str, List[str]] = False,
             comment: str = None, engine: Literal['pyarrow', 'auto'] = 'auto'):
    """
    loads a csv/parquet table into a pandas dataframe.

    :type filename: str or pathlib.Path
    :param filename: name of the csv file to be loaded
    :type index_col: int, default None
    :param index_col: number of column to be used as index. default is None, meaning no column will be used as index.
    :type drop_columns: str, list of str, or False (default False)
    :param drop_columns: if a string or list of strings are specified, \
    the columns of the same name/s will be dropped from the loaded DataFrame.
    :type squeeze: bool, default False
    :param squeeze: If the parsed data only contains one column then return a Series.
    :type comment: str (optional)
    :param comment: Indicates remainder of line should not be parsed. \
    If found at the beginning of a line, the line will be ignored altogether. This parameter must be a single character.
    :return: a pandas dataframe of the csv file
    """

    assert isinstance(filename,
                      (str, Path)), f"Filename must be of type str or pathlib.Path, is instead {type(filename)}."
    filename = Path(filename)
    assert filename.exists() and filename.is_file(), f"File '{filename.as_posix()}' does not exist!"
    assert filename.suffix.lower() in {'.csv', '.tsv', '.txt'}, \
        f"RNAlysis cannot load files of type '{filename.suffix}'. " \
        f"Please convert your file to a .csv, .tsv, .txt"

   
    kwargs = dict(sep=None, engine='python' if engine == 'auto' else engine, encoding='ISO-8859-1', comment=comment,
        skipinitialspace=True)
    if index_col is not None:
        kwargs['index_col'] = index_col
    df = pd.read_csv(filename, **kwargs)

    if index_col is not None:
        df.index = df.index.astype('str')
    df.index = [ind.strip() if isinstance(ind, str) else ind for ind in df.index]
    if isinstance(df, pd.DataFrame):
        df.columns = [col.strip() if isinstance(col, str) else str(col).strip() for col in df.columns]

        for col in df.columns:
            # check if the columns contains string data
            if pd.api.types.is_string_dtype(df[col]):
                df[col] = df[col].str.strip()
    else:
        if pd.api.types.is_string_dtype(df):
            df = df.str.strip()
    # if there remained only empty string "", change to Nan
    df = df.replace({"": np.nan})
    if drop_columns:
        drop_columns_lst = data_to_list(drop_columns)
        assert isinstanceiter(drop_columns_lst,
                                         str), f"'drop_columns' must be str, list of str, or False; " \
                                               f"is instead {type(drop_columns)}."
        for col in drop_columns_lst:
            col_stripped = col.strip()
            if col_stripped in df:
                df.drop(col_stripped, axis=1, inplace=True)
            else:
                raise IndexError(f"The argument {col} in 'drop_columns' is not a column in the loaded csv file!")
    return df

class ExploratoryAnalysis(object):

    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, list(str)] = None, is_normalized: bool = False):

        self.df = load_table(fname, drop_columns)
        self.is_normalized = is_normalized


    
if __name__ == "__main__":

    path = Path("/Users/ibrahimahmed/projects/GUI/counts_data.csv")
    df = load_table(path)
    print(df.head())