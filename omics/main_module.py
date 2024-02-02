
"""
This module can filter, normalize, intersect and visualize tabular data such as read counts and differential expression data.
Any tabular data saved in a csv format can be imported. \
Use this module to perform various filtering operations on your data,  normalize your data, \
perform set operations (union, intersection, etc), run basic exploratory analyses and plots \
(such as PCA, clustergram, violin plots, scatter, etc), \
save the filtered data to your computer, and more.
When you save filtered/modified data, its new file name will include by default \
 all of the operations performed on it, in the order they were performed, to allow easy traceback of your analyses.

"""
import copy
import os
import re
import types
from pathlib import Path
from typing import Any, Iterable, List, Tuple, Union, Callable, Sequence, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from grid_strategy import strategies
from scipy.stats import spearmanr
from scipy.stats.mstats import gmean
from sklearn.decomposition import PCA
from sklearn.preprocessing import PowerTransformer, StandardScaler
from tqdm.auto import tqdm

from .utils import generic, settings, validation, differential_expression, \
    param_typing
from .utils.generic import readable_name
from .utils.param_typing import *
import warnings
from pydantic import BaseModel, PositiveInt, NegativeInt, NonNegativeInt
from .utils import parsing, io


class Filter:
    """
    An all-purpose Filter object.


    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the DESeq output file contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.
    """
   
    def __init__(self, fname: Union[str, Path], drop_columns: Union[str, List[str]] = None):

        """
        Load a table.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)

        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/counted.csv")

        """
        # init from a file (load the csv into a DataFrame/Series)
        assert isinstance(fname, (str, Path))
        self.fname = Path(fname)
        self.df = io.load_table(fname, 0, squeeze=True, drop_columns=drop_columns)
        if isinstance(self.df, pd.Series):
            self.df = self.df.to_frame()
        # check for duplicate indices
        if self.df.index.has_duplicates:
            warnings.warn("This Filter object contains multiple rows with the same name/index.")

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, name: Union[str, Path]) -> 'Filter':
        obj = cls.__new__(cls)
        obj.df = df.copy(deep=True)
        obj.fname = Path(name)
        if obj.df.index.has_duplicates:
            warnings.warn("This Filter object contains multiple rows with the same name/index.")
        return obj

    def _init_warnings(self):
        pass

    def __repr__(self):
        return f"{type(self).__name__}('{self.fname.as_posix()}')"

    def __str__(self):
        return f"{self.readable_name} named '{self.fname.stem}{self.fname.suffix}'"

    def __len__(self):
        """
        Returns the number of rows in the Filter object

        """
        return self.shape[0]

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        if self.df.equals(other.df) and self.shape == other.shape:
            return True
        return False

    def __contains__(self, item):
        return True if item in self.df.index else False

    def __iter__(self):
        yield from self.df.index

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname)

    @property
    def columns(self) -> list:
        """
        The columns of df.

        :return: a list of the columns in the Filter object.
        :rtype: list
        """
        return list(self.df.columns)

    @property
    def shape(self) -> Tuple[int, int]:
        return self.df.shape

    def _update(self, **kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self, key, val)
            except AttributeError:
                raise AttributeError(f"Cannot update attribute {key} for {type(self)} object: attribute does not exist")

    def _inplace(self, new_df: pd.DataFrame, opposite: bool, inplace: bool, suffix: str,
                 printout_operation: str = 'filter', **filter_update_kwargs):

        """
        Executes the user's choice whether to filter in-place or create a new instance of the Filter object.

        :param new_df: the post-filtering DataFrame
        :type new_df: pd.DataFrame
        :param opposite: Determines whether to return the filtration ,or its opposite.
        :type opposite: bool
        :param inplace: Determines whether to filter in-place or not.
        :type inplace: bool
        :param suffix: The suffix to be added to the filename
        :type suffix: str
        :return: If inplace is False, returns a new instance of the Filter object.

        """
        legal_operations = {'filter': 'Filtering', 'normalize': 'Normalization', 'sort': 'Sorting',
                            'transform': 'Transformation', 'translate': 'Translation'}
        assert isinstance(inplace, bool), "'inplace' must be True or False!"
        assert isinstance(opposite, bool), "'opposite' must be True or False!"
        assert printout_operation.lower() in legal_operations, \
            f"Invalid input for variable 'printout_operation': {printout_operation}"
        # when user requests the opposite of a filter, return the Set Difference between the filtering result and self
        if opposite:
            new_df = self.df.loc[self.df.index.difference(new_df.index)]
            suffix += 'opposite'

        # update filename with the suffix of the operation that was just performed
        new_fname = Path(os.path.join(str(self.fname.parent), f"{self.fname.stem}{suffix}{self.fname.suffix}"))

        # generate printout for user ("Filtered X features, leaving Y... filtered inplace/not inplace")
        printout = ''
        if printout_operation.lower() == 'filter':
            printout += f"Filtered {self.shape[0] - new_df.shape[0]} features, leaving {new_df.shape[0]} " \
                        f"of the original {self.shape[0]} features. "
        elif printout_operation.lower() == 'translate':
            printout += f"Translated the gene IDs of {new_df.shape[0]} features. "
            if self.shape[0] != new_df.shape[0]:
                printout += f"Filtered {self.shape[0] - new_df.shape[0]} unmapped features, " \
                            f"leaving {new_df.shape[0]} of the original {self.shape[0]} features. "

        else:
            printout += f'{printout_operation.capitalize().rstrip("e")}ed {new_df.shape[0]} features. '
        # if inplace, modify the df, fname and other attributes of self
        if inplace:
            printout += f'{printout_operation.capitalize().rstrip("e")}ed inplace.'
            print(printout)
            self._update(df=new_df, fname=new_fname, **filter_update_kwargs)
        # if not inplace, copy self, modify the df/fname properties of the copy, and return it
        else:
            printout += f'{legal_operations[printout_operation]} result saved to new object.'
            print(printout)
            new_obj = self.__copy__()
            new_obj._update(df=new_df, fname=new_fname, **filter_update_kwargs)
            return new_obj

    def save_table(self, suffix: Literal['.csv', '.tsv', '.parquet'] = '.csv',
                   alt_filename: Union[None, str, Path] = None):

        """
        Save the current filtered data table.

        :param suffix: the file suffix
        :type suffix: '.csv', '.tsv', or '.parquet' (default='.csv')
        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        # save with the default filename if no alternative filename was given
        if alt_filename is None:
            alt_filename = self.fname.with_suffix(suffix)
        else:
            assert isinstance(alt_filename, (str, Path)), \
                f"'alt_filename' must be a string or Path object. Instead got {type(alt_filename)}."
            alt_filename = self.fname.parent.joinpath(alt_filename).with_suffix(suffix)
        io.save_table(self.df, alt_filename)

    def save_csv(self, alt_filename: Union[None, str, Path] = None):

        """
        Saves the current filtered data to a .csv file.

        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        self.save_table('.csv', alt_filename)

    def save_parquet(self, alt_filename: Union[None, str, Path] = None):

        """
        Saves the current filtered data to a .parquet file.

        :param alt_filename: If None, file name will be generated automatically \
        according to the filtering methods used. \
        If it's a string, it will be used as the name of the saved file. Example input: 'myfilename'
        :type alt_filename: str, pathlib.Path, or None (default)

        """
        self.save_table('.parquet', alt_filename)

    @staticmethod
    def _from_string(msg: str = '', delimiter: str = '\n'):

        """
        Takes a manual string input from the user, and then splits it using a delimiter into a list of values.

        :param msg: a promprt to be printed to the user
        :param delimiter: the delimiter used to separate the values. Default is '\n'
        :return: A list of the comma-seperated values the user inserted.
        """
        string = input(msg)
        split = string.split(sep=delimiter)
        if split[-1] == '':
            split = split[:-1]
        return split

    @readable_name('Table head')
    def head(self, n: PositiveInt = 5) -> pd.DataFrame:

        """
        Return the first n rows of the Filter object. See pandas.DataFrame.head documentation.

        :type n: positive int, default 5
        :param n: Number of rows to show.
        :return: returns the first n rows of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> d.head()
                               baseMean  log2FoldChange  ...         pvalue           padj
            WBGene00000002  6820.755327        7.567762  ...   0.000000e+00   0.000000e+00
            WBGene00000003  3049.625670        9.138071  ...  4.660000e-302  4.280000e-298
            WBGene00000004  1432.911791        8.111737  ...  6.400000e-237  3.920000e-233
            WBGene00000005  4028.154186        6.534112  ...  1.700000e-228  7.800000e-225
            WBGene00000006  1230.585240        7.157428  ...  2.070000e-216  7.590000e-213
            <BLANKLINE>
            [5 rows x 6 columns]

            >>> d.head(3) # return only the first 3 rows
                               baseMean  log2FoldChange  ...         pvalue           padj
            WBGene00000002  6820.755327        7.567762  ...   0.000000e+00   0.000000e+00
            WBGene00000003  3049.625670        9.138071  ...  4.660000e-302  4.280000e-298
            WBGene00000004  1432.911791        8.111737  ...  6.400000e-237  3.920000e-233
            <BLANKLINE>
            [3 rows x 6 columns]

        """
        return self.df.head(n)

    @readable_name('Table tail')
    def tail(self, n: PositiveInt = 5) -> pd.DataFrame:

        """
        Return the last n rows of the Filter object. See pandas.DataFrame.tail documentation.

        :type n: positive int, default 5
        :param n: Number of rows to show.
        :rtype: pandas.DataFrame
        :return: returns the last n rows of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> d.tail()
                               baseMean  log2FoldChange  ...        pvalue          padj
            WBGene00000025  2236.185837        2.477374  ...  1.910000e-81  1.460000e-78
            WBGene00000026   343.648987       -4.037191  ...  2.320000e-75  1.700000e-72
            WBGene00000027   175.142856        6.352044  ...  1.580000e-74  1.120000e-71
            WBGene00000028   219.163200        3.913657  ...  3.420000e-72  2.320000e-69
            WBGene00000029  1066.242402       -2.811281  ...  1.420000e-70  9.290000e-68
            <BLANKLINE>
            [5 rows x 6 columns]


            >>> d.tail(8) # returns the last 8 rows
                               baseMean  log2FoldChange  ...        pvalue          padj
            WBGene00000022   365.813048        6.101303  ...  2.740000e-97  2.400000e-94
            WBGene00000023  3168.566714        3.906719  ...  1.600000e-93  1.340000e-90
            WBGene00000024   221.925724        4.801676  ...  1.230000e-84  9.820000e-82
            WBGene00000025  2236.185837        2.477374  ...  1.910000e-81  1.460000e-78
            WBGene00000026   343.648987       -4.037191  ...  2.320000e-75  1.700000e-72
            WBGene00000027   175.142856        6.352044  ...  1.580000e-74  1.120000e-71
            WBGene00000028   219.163200        3.913657  ...  3.420000e-72  2.320000e-69
            WBGene00000029  1066.242402       -2.811281  ...  1.420000e-70  9.290000e-68
            <BLANKLINE>
            [8 rows x 6 columns]

        """
        return self.df.tail(n)

 
    def filter_duplicate_ids(self, keep: Literal['first', 'last', 'neither'] = 'first', opposite: bool = False,
                             inplace: bool = True):
        """
        Filter out rows with duplicate names/IDs (index).

        :param keep: determines which of the duplicates to keep for each group of duplicates. 'first' will keep the \
        first duplicate found for each group; 'last' will keep the last; \
        and 'neither' will remove *all* of the values in the group.
        :type keep: 'first', 'last', or 'neither' (default='first')
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = f'_dropduplicateskeep{keep}'
        if keep == 'neither':
            keep = False
        new_df = self.df[~self.df.index.duplicated(keep=keep)]
        return self._inplace(new_df, opposite, inplace, suffix)

    
    def filter_by_row_name(self, row_names: Union[str, List[str]], opposite: bool = False, inplace: bool = True):
        """
        Filter out specific rows from the table by their name (index).

        :param row_names: list of row names to be removed from the table.
        :type row_names: str or list of str
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = '_filterbyrowname'
        row_names = parsing.data_to_list(row_names)
        for name in row_names:
            assert name in self.df.index, f"'{name}' is now a row name in the table!"
        new_df = self.df.drop(row_names)
        return self._inplace(new_df, opposite, inplace, suffix)

    
    def drop_columns(self, columns: param_typing.ColumnNames, inplace: bool = True):
        """
        Drop specific columns from the table.

        :param columns: The names of the column/columns to be dropped fro mthe table.
        :type columns: str or list of str
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.
        """
        suffix = '_dropcolumns'
        columns = parsing.data_to_list(columns)
        for col in columns:
            assert col in self.columns, f"column '{col}' does not exist!"
        new_df = self.df.drop(columns, axis=1)
        return self._inplace(new_df, False, inplace, suffix, 'transform')


    
    def filter_percentile(self, percentile: param_typing.Fraction, column: param_typing.ColumnName,
                          opposite: bool = False, inplace: bool = True):

        """
        Removes all entries above the specified percentile in the specified column. \
        For example, if the column were 'pvalue' and the percentile was 0.5, then all features whose pvalue is above \
        the median pvalue will be filtered out.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new and filtered instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # keep only the rows whose value in the column 'log2FoldChange' is below the 75th percentile
            >>> d.filter_percentile(0.75,'log2FoldChange')
            Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.

            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> # keep only the rows vulse value in the column 'log2FoldChange' is above the 25th percentile
            >>> d.filter_percentile(0.25,'log2FoldChange',opposite=True)
            Filtered 7 features, leaving 21 of the original 28 features. Filtered inplace.


        """
        assert isinstance(percentile, (float, int)) and 0 <= percentile <= 1, \
            "percentile must be a float between 0 and 1!"
        assert isinstance(column, str) and column in self.columns, "Invalid column name!"
        suffix = f'_below{percentile}percentile'
        new_df = self.df[self.df[column] <= self.df[column].quantile(percentile)]
        return self._inplace(new_df, opposite, inplace, suffix)

   
    def split_by_percentile(self, percentile: param_typing.Fraction, column: param_typing.ColumnName) -> tuple:

        """
        Splits the features in the Filter object into two non-overlapping Filter objects: \
        one containing features below the specified percentile in the specfieid column, \
        and the other containing features about the specified percentile in the specified column.

        :type percentile: float between 0 and 1
        :param percentile: The percentile that all features above it will be filtered out.
        :type column: str
        :param column: Name of the DataFrame column according to which the filtering will be performed.
        :rtype: Tuple[filtering.Filter, filtering.Filter]
        :return: a tuple of two Filter objects: the first contains all of the features below the specified percentile, \
        and the second contains all of the features above and equal to the specified percentile.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/test_deseq.csv")
            >>> below, above = d.split_by_percentile(0.75,'log2FoldChange')
            Filtered 7 features, leaving 21 of the original 28 features. Filtering result saved to new object.
            Filtered 21 features, leaving 7 of the original 28 features. Filtering result saved to new object.

        """
        return self.filter_percentile(percentile=percentile, column=column, opposite=False,
                                      inplace=False), self.filter_percentile(percentile=percentile, column=column,
                                                                             opposite=True, inplace=False)

   



    @readable_name('Filter with a number filter')
    def number_filters(self, column: param_typing.ColumnName,
                       operator: Literal['greater than', 'equals', 'lesser than'], value: float,
                       opposite: bool = False, inplace: bool = True):

        """
        Applay a number filter (greater than, equal, lesser than) on a particular column in the Filter object.

        :type column: str
        :param column: name of the column to filter by
        :type operator: str: 'gt' / 'greater than' / '>', 'eq' / 'equals' / '=', 'lt' / 'lesser than' / '<'
        :param operator: the operator to filter the column by (greater than, equal or lesser than)
        :type value: float
        :param value: the value to filter by
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current Filter object. If False, \
        the function will return a new Filter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','gt',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','greater than',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

            >>> filt = filtering.Filter('tests/test_files/test_deseq.csv')
            >>> #keep only rows that have a value greater than 5900 in the column 'baseMean'.
            >>> filt.number_filters('baseMean','>',5900)
            Filtered 26 features, leaving 2 of the original 28 features. Filtered inplace.

        """
        # determine whether operator is valid
        operator_dict = {'gt': 'gt', 'greater than': 'gt', '>': 'gt', 'eq': 'eq', 'equals': 'eq', '=': 'eq', 'lt': 'lt',
                         'lesser than': 'lt', '<': 'lt', 'equal': 'eq'}
        operator = operator.lower()
        assert operator in operator_dict, f"Invalid operator {operator}"
        op = operator_dict[operator]
        # determine that 'value' is a number
        assert isinstance(value, (int, float)), "'value' must be a number!"
        # determine that the column is legal
        assert column in self.columns, f"column {column} not in DataFrame!"

        suffix = f"_{column}{op}{value}"
        # perform operation according to operator
        if op == 'eq':
            new_df = self.df[self.df[column] == value]
        elif op == 'gt':
            new_df = self.df[self.df[column] > value]
        elif op == 'lt':
            new_df = self.df[self.df[column] < value]

        # noinspection PyUnboundLocalVariable
        return self._inplace(new_df, opposite, inplace, suffix)

   
#*********************************
@readable_name('Differential expression table')
class DESeqFilter(Filter):
    """
    A class that receives a DESeq output file and can filter it according to various characteristics.

    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the DESeq output file contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.

    """
    __slots__ = {'log2fc_col': 'name of the log2 fold change column', 'padj_col': 'name of the adjusted p-value column'}

    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, List[str]] = None,
                 log2fc_col: str = 'log2FoldChange', padj_col: str = 'padj', suppress_warnings: bool = False):
        """
        Load a differential expression table. A valid differential expression table should have \
        a column containing log2(fold change) values for each gene, and another column containing \
        adjusted p-values for each gene.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)
        :param log2fc_col: name of the table column containing log2(fold change) values.
        :type log2fc_col: str (default='Log2FoldChange')
        :param padj_col: name of the table column containing adjusted p-values.
        :type padj_col: str (default='padj')
        :param suppress_warnings: if True, RNAlysis will not issue warnings about the loaded table's \
        structure or content.
        :type suppress_warnings: bool (default=False)
        """
        super().__init__(fname, drop_columns)
        self.log2fc_col = log2fc_col
        self.padj_col = padj_col
        if not suppress_warnings:
            self._init_warnings()

    def _init_warnings(self):
        if self.log2fc_col not in self.columns:
            warnings.warn(f"The specified log2fc_col '{self.log2fc_col}' does not appear in the DESeqFilter's columns: "
                          f"{self.columns}. DESeqFilter-specific functions that depend on "
                          f"log2(fold change) may fail to run. ")
        if self.padj_col not in self.columns:
            warnings.warn(f"The specified padj_col '{self.padj_col}' does not appear in the DESeqFilter's columns: "
                          f"{self.columns}. DESeqFilter-specific functions that depend on p-values may fail to run. ")

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, name: Union[str, Path],
                       log2fc_col: str = 'log2FoldChange', padj_col: str = 'padj',
                       suppress_warnings: bool = False) -> 'DESeqFilter':
        obj = cls.__new__(cls)
        obj.df = df.copy(deep=True)
        obj.fname = Path(name)
        obj.log2fc_col = log2fc_col
        obj.padj_col = padj_col
        if not suppress_warnings:
            obj._init_warnings()
        return obj

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname, self.log2fc_col, self.padj_col, suppress_warnings=True)

    def _assert_padj_col(self):
        if self.padj_col not in self.columns:
            raise KeyError(f"A column with adjusted p-values under the name padj_col='{self.padj_col}' "
                           f"could not be found. Try setting a different value for the parameter 'padj_col' "
                           f"when creating the DESeqFilter object.")

    def _assert_log2fc_col(self):
        if self.log2fc_col not in self.columns:
            raise KeyError(f"A column with log2 fold change values under the name log2fc_col='{self.log2fc_col}' "
                           f"could not be found. Try setting a different value for the parameter 'log2fc_col' "
                           f"when creating the DESeqFilter object.")

    @readable_name('Filter by statistical significance')
    def filter_significant(self, alpha: param_typing.Fraction = 0.1, opposite: bool = False, inplace: bool = True):

        """
        Removes all features which did not change significantly, according to the provided alpha.

        :param alpha: the significance threshold to determine which genes will be filtered. between 0 and 1.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_significant(0.1) # keep only rows whose adjusted p-value is <=0.1
            Filtered 4 features, leaving 25 of the original 29 features. Filtered inplace.

             >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_significant(0.1, opposite=True) # keep only rows whose adjusted p-value is >0.1
            Filtered 25 features, leaving 4 of the original 29 features. Filtered inplace.

        """
        assert isinstance(alpha, float), "alpha must be a float!"
        self._assert_padj_col()

        new_df = self.df[self.df[self.padj_col] <= alpha]
        suffix = f"_sig{alpha}"
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by absolute log2 fold-change magnitude')
    def filter_abs_log2_fold_change(self, abslog2fc: float = 1, opposite: bool = False, inplace: bool = True):

        """
        Filters out all features whose absolute log2 fold change is below the indicated threshold. \
        For example: if log2fc is 2.0, all features whose log2 fold change is between 1 and -1 (went up less than \
        two-fold or went down less than two-fold) will be filtered out.

        :param abslog2fc: The threshold absolute log2 fold change for filtering out a feature. Float or int. \
        All features whose absolute log2 fold change is lower than log2fc will be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_abs_log2_fold_change(2) # keep only rows whose log2(fold change) is >=2 or <=-2
            Filtered 1 features, leaving 28 of the original 29 features. Filtered inplace.

        """
        assert isinstance(abslog2fc, (float, int)), "abslog2fc must be a number!"
        assert abslog2fc >= 0, "abslog2fc must be non-negative!"
        self._assert_log2fc_col()

        suffix = f"_{abslog2fc}abslog2foldchange"
        new_df = self.df[np.abs(self.df[self.log2fc_col]) >= abslog2fc]
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Filter by log2 fold-change direction')
    def filter_fold_change_direction(self, direction: Literal['pos', 'neg'] = 'pos', opposite: bool = False,
                                     inplace: bool = True):

        """
        Filters out features according to the direction in which they changed between the two conditions.

        :param direction: 'pos' or 'neg'. If 'pos', will keep only features that have positive log2foldchange. \
        If 'neg', will keep only features that have negative log2foldchange.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current DESeqFilter object. If False, \
        the function will return a new DESeqFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of DESeqFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('pos') # keep only rows with a positive log2(fold change) value
            Filtered 3 features, leaving 26 of the original 29 features. Filtered inplace.

            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('neg') # keep only rows with a negative log2(fold change) value
            Filtered 27 features, leaving 2 of the original 29 features. Filtered inplace.

            >>> d = filtering.DESeqFilter('tests/test_files/sample_deseq.csv')
            >>> d.filter_fold_change_direction('pos', opposite=True) # keep only rows with a non-positive log2(fold change) value
            Filtered 26 features, leaving 3 of the original 29 features. Filtered inplace.

        """
        assert isinstance(direction, str), \
            "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. "
        self._assert_log2fc_col()

        if direction == 'pos':
            new_df = self.df[self.df[self.log2fc_col] > 0]
            suffix = '_PositiveLog2FC'
        elif direction == 'neg':
            new_df = self.df[self.df[self.log2fc_col] < 0]
            suffix = '_NegativeLog2FC'
        else:
            raise ValueError(
                "'direction' must be either 'pos' for positive fold-change, or 'neg' for negative fold-change. ")
        return self._inplace(new_df, opposite, inplace, suffix)

    @readable_name('Split by log2 fold-change direction')
    def split_fold_change_direction(self) -> tuple:

        """
        Splits the features in the DESeqFilter object into two non-overlapping DESeqFilter \
        objects, based on the direction of their log2foldchange. The first object will contain only features with a \
        positive log2foldchange, the second object will contain only features with a negative log2foldchange.

        :rtype: Tuple[filtering.DESeqFilter, filteirng.DESeqFilter]
        :return: a tuple containing two DESeqFilter objects: the first has only features with positive log2 fold change, \
        and the other has only features with negative log2 fold change.


        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.DESeqFilter('tests/test_files/test_deseq.csv')
            >>> pos, neg = d.split_fold_change_direction()
            Filtered 2 features, leaving 26 of the original 28 features. Filtering result saved to new object.
            Filtered 26 features, leaving 2 of the original 28 features. Filtering result saved to new object.

        """
        return self.filter_fold_change_direction(direction='pos', inplace=False), self.filter_fold_change_direction(
            direction='neg', inplace=False)

    @readable_name('Volcano plot')
    def volcano_plot(self, alpha: param_typing.Fraction = 0.1, log2fc_threshold: Union[float, None] = 1,
                     title: Union[str, Literal['auto']] = 'auto', title_fontsize: float = 16,
                     label_fontsize: float = 16, tick_fontsize: float = 12,
                     annotation_fontsize: float = 10, point_size: float = 10, opacity: param_typing.Fraction = 0.65,
                     interactive: bool = True, show_cursor: bool = False) -> plt.Figure:

        """
        Plots a volcano plot (log2(fold change) vs -log10(adj. p-value)) of the DESeqFilter object. \
        Significantly upregulated features are colored in red, \
        and significantly downregulated features are colored in blue. \
        If the plot is generated in interactive mode, data points can be labeled by clicking on them.

        :type alpha: float between 0 and 1
        :param alpha: the significance threshold to paint data points as significantly up/down-regulated.
        :param log2fc_threshold: the absolute log2(fold change) threshold to paint data as \
        significantly up/down-regulated. if log2fc_threshold is None, no threshold will be used.
        :type log2fc_threshold: non-negative float or None (default=1)
        :param title: The title of the plot. If 'auto', a title will be generated automatically.
        :type title: str or 'auto' (default='auto')
        :param title_fontsize: determines the font size of the graph title.
        :type title_fontsize: float (default=30)
        :param label_fontsize: determines the font size of the X and Y axis labels.
        :type label_fontsize: float (default=15)
         :param tick_fontsize: determines the font size of the X and Y tick labels.
        :type tick_fontsize: float (default=10)
        :param annotation_fontsize: determines the font size of the point annotations created in interactive mode.
        :type annotation_fontsize: float (default=10)
        :param opacity: float between 0 and 1 (default=0.65)
        :type opacity: determines the opacity of the points in the scatter plot. 0 indicates completely transparent, \
        while 1 indicates completely opaque.
        :param point_size: determines the size of the points in the scatter plot
        :type point_size: float (default=10)
        :param interactive: if True, turns on interactive mode. While in interactive mode, you can click on a data \
        point to label it with its gene name/ID, or click on a labeled data point to unlabel it.
        :type interactive: bool (default=True)
        :param show_cursor: if True, show the cursor position on the plot during interactive mode
        :type show_cursor: bool (default=False)
        :rtype: A matplotlib Figure

        .. figure:: /figures/volcano.png
           :align:   center
           :scale: 70 %

           Example plot of volcano_plot()

        """
        self._assert_padj_col()
        self._assert_log2fc_col()

        if log2fc_threshold is None:
            log2fc_threshold = 0
        else:
            assert isinstance(log2fc_threshold, (int, float)) and log2fc_threshold >= 0, \
                "'log2fc_threshold' must be a non-negative number!"

        if interactive:
            fig = plt.figure(constrained_layout=True, FigureClass=generic.InteractiveScatterFigure,
                             labels=parsing.data_to_tuple(self.df.index), annotation_fontsize=annotation_fontsize,
                             show_cursor=show_cursor)
            ax = fig.axes[0]
            kwargs = {'picker': 5}
        else:
            fig = plt.figure(constrained_layout=True)
            ax = fig.add_subplot(111)
            kwargs = {}
        colors = pd.Series(index=self.df.index, dtype='float64')
        colors.loc[(self.df[self.padj_col] <= alpha) & (self.df[self.log2fc_col] > log2fc_threshold)] = 'tab:red'
        colors.loc[(self.df[self.padj_col] <= alpha) & (self.df[self.log2fc_col] < -log2fc_threshold)] = 'tab:blue'
        colors.fillna('grey', inplace=True)
        ax.scatter(self.df[self.log2fc_col], -np.log10(self.df[self.padj_col]), c=colors, s=point_size, alpha=opacity,
                   **kwargs)
        if title == 'auto':
            title = f"Volcano plot of {self.fname.stem}"
        ax.set_title(title, fontsize=title_fontsize)
        ax.set_xlabel(r"$\log_2$(Fold Change)", fontsize=label_fontsize)
        ax.set_ylabel(r'-$\log_{10}$(Adj. P-Value)', fontsize=label_fontsize)
        ax.tick_params(axis='both', which='both', labelsize=tick_fontsize)

        if log2fc_threshold > 0:
            ax.axvline(log2fc_threshold, linestyle='--', color='black')
            ax.axvline(-log2fc_threshold, linestyle='--', color='black')
        generic.despine(ax)

        plt.show()
        return fig


#*********************************

class CountFilter(Filter):
    """
    A class that receives a count matrix and can filter it according to various characteristics.

    **Attributes**

    df: pandas DataFrame
        A DataFrame that contains the count matrix contents. \
        The DataFrame is modified upon usage of filter operations.
    shape: tuple (rows, columns)
        The dimensions of df.
    columns: list
        The columns of df.
    fname: pathlib.Path
        The path and filename for the purpose of saving df as a csv file. \
        Updates automatically when filter operations are applied.
    index_set: set
        All of the indices in the current DataFrame (which were not removed by previously used filter methods) \
        as a set.
    index_string: string
        A string of all feature indices in the current DataFrame separated by newline.
    triplicates: list
        Returns a nested list of the column names in the CountFilter, grouped by alphabetical order into triplicates. \
        For example, if counts.columns is ['A_rep1','A_rep2','A_rep3','B_rep1','B_rep2',_B_rep3'], then \
        counts.triplicates will be  [['A_rep1','A_rep2','A_rep3'],['B_rep1','B_rep2',_B_rep3']]

    """
    def __init__(self, fname: Union[str, Path, tuple], drop_columns: Union[str, List[str]] = None,
                 is_normalized: bool = False):
        """
        Load a count matrix. A valid count matrix should have one row per gene/genomic feature \
        and one column per condition/RNA library. The contents of the count matrix can be raw or pre-normalized.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)
        :param is_normalized: indicates whether this count table is pre-normalized. \
        RNAlysis issues a warning when a function meant for normalized tables is applied to a \
        table that was not already normalized.
        :type is_normalized: bool (default=False)
        """
        super().__init__(fname, drop_columns)
        self._is_normalized = is_normalized

    def _init_warnings(self):
        if len(self._numeric_columns) < len(self.columns):
            warnings.warn(f"The following columns in the CountFilter are not numeric, and will therefore be ignored "
                          f"when running some CountFilter-specific functions: "
                          f"{set(self.columns).difference(self._numeric_columns)}")

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, name: Union[str, Path], is_normalized: bool = False,
                       suppress_warnings: bool = False) -> 'CountFilter':
        obj = cls.__new__(cls)
        obj.df = df.copy(deep=True)
        obj.fname = Path(name)
        obj._is_normalized = is_normalized
        if not suppress_warnings:
            obj._init_warnings()
        return obj

    def __copy__(self):
        return type(self).from_dataframe(self.df, self.fname, self.is_normalized, suppress_warnings=True)

  
    def is_normalized(self) -> bool:
        return self._is_normalized

   
    def _numeric_columns(self) -> list:
        """
        Returns a list of the numeric (int/float) columns in the DataFrame.
        """
        return list(self.df.columns[[dtype in self._numeric_dtypes for dtype in self.df.dtypes]])

   
    def triplicates(self):

        """
        Returns a nested list of the column names in the CountFilter, grouped by alphabetical order into triplicates. \
        For example, if counts.columns is ['A_rep1','A_rep2','A_rep3','B_rep1','B_rep2',_B_rep3'], then \
        counts.triplicates will be  [['A_rep1','A_rep2','A_rep3'],['B_rep1','B_rep2',_B_rep3']]

        """

        multiplier = 3
        numeric_cols = sorted(self._numeric_columns)
        n_cols = len(numeric_cols)
        triplicate = [numeric_cols[i * multiplier:(1 + i) * multiplier] for i in range(n_cols // multiplier)]
        if len(numeric_cols[(n_cols // multiplier) * multiplier::]) > 0:
            triplicate.append(numeric_cols[(n_cols // multiplier) * multiplier::])
            warnings.warn(f'Number of samples {n_cols} is not divisible by 3. '
                          f'Appending the remaining {n_cols % multiplier} as an inncomplete triplicate.')
        return triplicate

    def _diff_exp_assertions(self, design_mat_df: pd.DataFrame):
        assert design_mat_df.shape[0] == self.shape[1], f"The number of items in the design matrix " \
                                                        f"({design_mat_df.shape[0]}) does not match the number of " \
                                                        f"columns in the count matrix ({self.shape[1]})."
        assert sorted(design_mat_df.index) == sorted(self.columns), f"The sample names in the design matrix do not " \
                                                                    f"match the sample names in the count matrix: " \
                                                                    f"{sorted(design_mat_df.index)} != " \
                                                                    f"{sorted(self.columns)}"

        for factor in design_mat_df.columns:
            assert generic.sanitize_variable_name(
                factor) == factor, f"Invalid factor name '{factor}': contains invalid characters." \
                                   f" \nSuggested alternative name: '{generic.sanitize_variable_name(factor)}'. "

   
   
    def differential_expression_edger2(self, design_matrix: Union[str, Path],
                                       comparisons: Iterable[Tuple[str, str, str]],
                                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                       output_folder: Union[str, Path, None] = None
                                       ) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_todays_cache_dir().joinpath(self.fname.name)
        design_mat_path = None
        i = 0
        while design_mat_path is None or design_mat_path.exists():
            design_mat_path = io.get_todays_cache_dir().joinpath(f'design_mat_{i}.csv')
            i += 1

        io.save_table(self.df.round(), data_path)
        # use Pandas to automatically detect file delimiter type, then export it as a CSV file.
        design_mat_df = io.load_table(design_matrix, index_col=0)
        self._diff_exp_assertions(design_mat_df)
        io.save_table(design_mat_df, design_mat_path)
        differential_expression.run_edger2_analysis(data_path, design_mat_path, comparisons,
                                                                   r_installation_folder)
        '''
        r_output_dir = differential_expression.run_edger2_analysis(data_path, design_mat_path, comparisons,
                                                                   r_installation_folder)
        outputs = []
        for item in r_output_dir.iterdir():
            if not item.is_file():
                continue
            if item.suffix == '.csv':
                outputs.append(DESeqFilter(item))

            if output_folder is not None:
                with open(item) as infile, open(output_folder.joinpath(item.name), 'w') as outfile:
                    outfile.write(infile.read())

        return parsing.data_to_tuple(outputs)
        '''

    def differential_expression_edger3(self, design_matrix: Union[str, Path],
                                       reference: str,
                                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                       output_folder: Union[str, Path, None] = None
                                       ) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_todays_cache_dir().joinpath(self.fname.name)
        design_mat_path = None
        i = 0
        while design_mat_path is None or design_mat_path.exists():
            design_mat_path = io.get_todays_cache_dir().joinpath(f'design_mat_{i}.csv')
            i += 1

        io.save_table(self.df.round(), data_path)
        # use Pandas to automatically detect file delimiter type, then export it as a CSV file.
        design_mat_df = io.load_table(design_matrix, index_col=0)
        self._diff_exp_assertions(design_mat_df)
        io.save_table(design_mat_df, design_mat_path)
        differential_expression.run_edger3_analysis(data_path, design_mat_path, reference,
                                                                   r_installation_folder)
      

    def differential_expression_deseq2(self, design_matrix: Union[str, Path],
                                       comparisons: Iterable[Tuple[str, str, str]],
                                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                       output_folder: Union[str, Path, None] = None
                                       ) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_todays_cache_dir().joinpath(self.fname.name)
        design_mat_path = None
        i = 0
        while design_mat_path is None or design_mat_path.exists():
            design_mat_path = io.get_todays_cache_dir().joinpath(f'design_mat_{i}.csv')
            i += 1

        io.save_table(self.df.round(), data_path)
        # use Pandas to automatically detect file delimiter type, then export it as a CSV file.
        design_mat_df = io.load_table(design_matrix, index_col=0)
        self._diff_exp_assertions(design_mat_df)
        io.save_table(design_mat_df, design_mat_path)

        r_output_dir = differential_expression.run_deseq2_analysis(data_path, design_mat_path, comparisons,
                                                                   r_installation_folder)
        outputs = []
        for item in r_output_dir.iterdir():
            if not item.is_file():
                continue
            if item.suffix == '.csv':
                outputs.append(DESeqFilter(item))

            if output_folder is not None:
                with open(item) as infile, open(output_folder.joinpath(item.name), 'w') as outfile:
                    outfile.write(infile.read())

        return parsing.data_to_tuple(outputs)
 
    def differential_expression_pathway(self,
                                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                                       output_folder: Union[str, Path, None] = None
                                       ) -> Tuple['DESeqFilter', ...]:
        """
        Run differential expression analysis on the count matrix using the \
        `DESeq2 <https://doi.org/10.1186/s13059-014-0550-8>`_ algorithm. \
        The count matrix you are analyzing should be unnormalized (meaning, raw read counts). \
        The analysis will be based on a design matrix supplied by the user. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix. \
        To run this function, a version of R must be installed.

        :param design_matrix: path to a csv file containing the experiment's design matrix. \
        The design matrix should contain at least two columns: the first column contains all the sample names, \
        and each of the following columns contains an experimental design factor (e.g. 'condition', 'replicate', etc). \
        (see the User Guide and Tutorial for a complete example). \
        The analysis formula will contain all the factors in the design matrix.
        :type design_matrix: str or Path
        :param comparisons: specifies what comparisons to build results tables out of. \
        each individual comparison should be a tuple with exactly three elements: \
        the name of a factor in the design formula, the name of the numerator level for the fold change, \
        and the name of the denominator level for the fold change.
        :type comparisons: Iterable of tuple(factor, numerator_value, denominator_value)
        :param r_installation_folder: Path to the installation folder of R. For example: \
        'C:/Program Files/R/R-4.2.1'
        :type r_installation_folder: str, Path, or 'auto' (default='auto')
        :param output_folder: Path to a folder in which the analysis results, \
        as well as the log files and R script used to generate them, will be saved. \
        if output_folder is None, the results will not be saved to a specified directory.
        :type output_folder: str, Path, or None
        :return: a tuple of DESeqFilter objects, one for each comparison
        """
        if output_folder is not None:
            output_folder = Path(output_folder)
            assert output_folder.exists(), 'Output folder does not exist!'

        self._validate_is_normalized(expect_normalized=False)
        data_path = io.get_data_dir().joinpath(self.fname.name)
    
        io.save_table(self.df.round(), data_path)
    
        r_output_dir = differential_expression.run_pathway_analysis(data_path, output_folder,
                                                                   r_installation_folder)
        outputs = []
        for item in r_output_dir.iterdir():
            if not item.is_file():
                continue
            if item.suffix == '.csv':
                outputs.append(DESeqFilter(item))

            if output_folder is not None:
                with open(item) as infile, open(output_folder.joinpath(item.name), 'w') as outfile:
                    outfile.write(infile.read())

        return parsing.data_to_tuple(outputs)

 
    

    def _validate_is_normalized(self, expect_normalized: bool = True):
        if not self.is_normalized and expect_normalized:
            warnings.warn("This function is meant for normalized values, and your count matrix may not be normalized. ")
        elif self.is_normalized and not expect_normalized:
            warnings.warn(
                "This function is meant for raw, unnormalize counts, and your count matrix appears to be normalized. ")

    def _norm_scaling_factors(self, scaling_factors: Union[pd.DataFrame, pd.Series]):
        numeric_cols = self._numeric_columns
        scaling_factors = scaling_factors.squeeze()
        new_df = self.df.copy()

        if isinstance(scaling_factors, pd.Series):
            assert scaling_factors.shape[0] == len(numeric_cols), \
                f"Number of scaling factors ({scaling_factors.shape[0]}) does not match " \
                f"number of numeric columns in your data table ({len(numeric_cols)})!"

            for column in new_df.columns:
                if column in numeric_cols:
                    norm_factor = scaling_factors[column]
                    new_df[column] /= norm_factor
        elif isinstance(scaling_factors, pd.DataFrame):
            assert scaling_factors.shape[0] >= self.shape[0] and scaling_factors.shape[1] == len(numeric_cols), \
                f"Dimensions of scaling factors table ({scaling_factors.shape}) does not match the " \
                f"dimensions of your data table ({(self.shape[0], len(numeric_cols))} - numeric columns only)!"
            new_df[numeric_cols] = new_df[numeric_cols].div(scaling_factors, axis='rows')

        return new_df

   
    def normalize_to_rpm_htseqcount(self, special_counter_fname: Union[str, Path], inplace: bool = True,
                                    return_scaling_factors: bool = False):

        """
        Normalizes the count matrix to Reads Per Million (RPM). \
        Uses a table of feature counts (ambiguous, no feature, not aligned, etc) from HTSeq-count's output. \
        Divides each column in the CountFilter object by (total reads + ambiguous + no feature)*10^-6 .

        :param special_counter_fname: the .csv file which contains feature information about the RNA library \
        (ambiguous, no feature, not aligned, etc).
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_rpm_htseqcount("tests/test_files/uncounted.csv")
           Normalized 22 features. Normalized inplace.

        """
        suffix = '_normtoRPMhtseqcount'
        scaling_factors = {}

        if isinstance(special_counter_fname, (str, Path)):
            features = io.load_table(special_counter_fname, 0)
        elif isinstance(special_counter_fname, pd.DataFrame):
            features = special_counter_fname
        else:
            raise TypeError("Invalid type for 'special_counter_fname'!")
        numeric_cols = self._numeric_columns
        for column in self.df.columns:
            if column in numeric_cols:
                norm_factor = (self.df[column].sum() + features.loc[r'__ambiguous', column] + features.loc[
                    r'__no_feature', column] + features.loc[r'__alignment_not_unique', column]) / (10 ** 6)
                scaling_factors[column] = norm_factor

        scaling_factors = pd.Series(scaling_factors)
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)


    def normalize_to_rpm(self, inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Reads Per Million (RPM). \
        Divides each column in the count matrix by (total reads)*10^-6 .

        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_rpm()
           Normalized 22 features. Normalized inplace.

        """
        suffix = '_normtoRPM'
        scaling_factors = {}

        numeric_cols = self._numeric_columns
        for column in self.df.columns:
            if column in numeric_cols:
                norm_factor = self.df[column].sum() / (10 ** 6)
                scaling_factors[column] = norm_factor

        scaling_factors = pd.Series(scaling_factors)
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    def normalize_to_rpkm(self, gtf_file: Union[str, Path], feature_type: Literal['gene', 'transcript'] = 'gene',
                          method: Literal[LEGAL_GENE_LENGTH_METHODS] = 'mean', inplace: bool = True,
                          return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Reads Per Kilobase Million (RPKM). \
        Divides each column in the count matrix by (total reads)*(gene length / 1000)*10^-6. \


       :param gtf_file: Path to a GTF/GFF3 annotation file. This file will be used to determine the length of each \
        gene/transcript. The gene/transcript names in this annotation file should match the ones in count matrix.
        :type gtf_file: str or Path
        :param feature_type: the type of features in your count matrix. if feature_type is 'transcript', \
        lengths will be calculated per-transcript, and the 'method' parameter is ignored. \
        Otherwise, lengths will be aggregated per gene according to the method specified in the 'method' parameter.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :param method: if feature_type='gene', this determines the aggregation method to calculate gene lengths. \
        'mean', 'median', 'min', and 'max' will calculate the \
        mean/median/min/max of all transcripts' lengths of the given gene. \
        'geometric_mean' will calculate the goemetric mean of all transcripts' lengths of the given gene. \
        'merged_exons' will calculate the total lengths of all exons of a gene across all of its transcripts, \
        while counting overlapping exons/regions exactly once.
        :type method: 'mean', 'median', 'min', 'max', 'geometric_mean', or 'merged_exons' (deafult='mean')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.
        """
        suffix = f"_normtoRPKM{method.replace('_', '')}"
        gene_lengths_kbp = pd.Series(
            genome_annotation.get_genomic_feature_lengths(gtf_file, feature_type, method)) / 1000
        numeric_cols = self._numeric_columns
        scaling_factors = pd.concat([gene_lengths_kbp] * len(numeric_cols), axis=1)
        scaling_factors.columns = numeric_cols
        scaling_factors = scaling_factors.mul(self.df[numeric_cols].sum() / (10 ** 6), axis='columns')

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)


    def normalize_to_tpm(self, gtf_file: Union[str, Path], feature_type: Literal['gene', 'transcript'] = 'gene',
                         method: Literal[LEGAL_GENE_LENGTH_METHODS] = 'mean', inplace: bool = True,
                         return_scaling_factors: bool = False):
        """
        Normalizes the count matrix to Transcripts Per Million (TPM). \
        First, normalizes each gene to Reads Per Kilobase (RPK) \
        by dividing each gene in the count matrix by its length in Kbp (gene length / 1000). \
        Then, divides each column in the RPK matrix by (total RPK in column)*10^-6. \
        This calculation is similar to that of Reads Per Kilobase Million (RPKM), but in the opposite order: \
        the "per million" normalization factors are calculated **after** normalizing to gene lengths, not before.

       :param gtf_file: Path to a GTF/GFF3 annotation file. This file will be used to determine the length of each \
        gene/transcript. The gene/transcript names in this annotation file should match the ones in count matrix.
        :type gtf_file: str or Path
        :param feature_type: the type of features in your count matrix. if feature_type is 'transcript', \
        lengths will be calculated per-transcript, and the 'method' parameter is ignored. \
        Otherwise, lengths will be aggregated per gene according to the method specified in the 'method' parameter.
        :type feature_type: 'gene' or 'transcript' (default='gene')
        :param method: if feature_type='gene', this determines the aggregation method to calculate gene lengths. \
        'mean', 'median', 'min', and 'max' will calculate the \
        mean/median/min/max of all transcripts' lengths of the given gene. \
        'geometric_mean' will calculate the goemetric mean of all transcripts' lengths of the given gene. \
        'merged_exons' will calculate the total lengths of all exons of a gene across all of its transcripts, \
        while counting overlapping exons/regions exactly once.
        :type method: 'mean', 'median', 'min', 'max', 'geometric_mean', or 'merged_exons' (deafult='mean')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.
        """
        suffix = f"_normtoTPM{method.replace('_', '')}"
        gene_lengths_kbp = pd.Series(
            genome_annotation.get_genomic_feature_lengths(gtf_file, feature_type, method)) / 1000
        tmp_df = self.df.copy()
        numeric_cols = self._numeric_columns
        tmp_df[numeric_cols] = tmp_df[numeric_cols].div(gene_lengths_kbp, axis='rows')

        scaling_factors = pd.concat([gene_lengths_kbp] * len(numeric_cols), axis=1)
        scaling_factors.columns = numeric_cols
        scaling_factors = scaling_factors.mul(tmp_df[numeric_cols].sum() / (10 ** 6), axis='columns')

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    def normalize_to_quantile(self, quantile: param_typing.Fraction = 0.75, inplace: bool = True,
                              return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the quantile method, generalized from \
        `Bullard et al 2010 <https://doi.org/10.1186/1471-2105-11-94>`_. \
        This is the default normalization method used by R's Limma. \
        To calculate the Quantile Method scaling factors, you first calculate the given quantile of gene expression \
        within each sample, excluding genes that have 0 reads in all samples. \
        You then divide those quantile values by the total number of reads in each sample, \
        which yields the scaling factors for each sample.

        :param quantile: the quantile from which scaling factors will be calculated.
        :type quantile: float between 0 and 1 (default=0.75)
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_to_quantile(0.75)
           Normalized 22 features. Normalized inplace.
        """
        suffix = f'_normto{quantile}quantile'
        data = self.df[self._numeric_columns]
        expressed_genes = data[data.sum(axis=1) != 0]
        quantiles = expressed_genes.quantile(quantile, axis=0)
        if quantiles.min() == 0:
            warnings.warn("One or more quantiles are zero")

        scaling_factors = quantiles / quantiles.mean()
        new_df = self._norm_scaling_factors(scaling_factors)

        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    def normalize_tmm(self, log_ratio_trim: float = 0.3, sum_trim: float = 0.05,
                      a_cutoff: Union[float, None] = -1 * 10 ** 10,
                      ref_column: Union[Literal['auto'], param_typing.ColumnName] = 'auto',
                      inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'trimmed mean of M values' (TMM) method \
        `(Robinson and Oshlack 2010) <https://doi.org/10.1186/gb-2010-11-3-r25>`_. \
        This is the default normalization method used by R's edgeR. \
        To calculate the Trimmed Mean of M Values scaling factors, you first calculate the M-values of each gene \
        between each sample and the reference sample (log2 of each sample Minus log2 of the reference sample), \
        and the A-values of each gene between each sample and the reference sample \
        (log2 of each sample Added to log2 of the reference sample). \
        You then trim out genes with extreme values that are likely to be differentially expressed or non-indicative, \
        by trimming the top and bottom X% of M-values, the top and bottom Y% of A-values, all A-values which are \
        smaller than the specified cutuff, and all genes with 0 reads (to avoid log2 values of inf or -inf). \
        Next, a weighted mean is calculated on the filtered M-values, with the weights being an inverse of \
        an approximation of variance of each gene, which gives out the scaling factors for each sample. \
        Finally, the scaling factors are adjusted, for symmetry, so that they multiply to 1.

        :param log_ratio_trim: the fraction of M-values that should be trimmed from each direction (top and bottom X%).
        :type log_ratio_trim:  float between 0 and 0.5 (default=0.3)
        :param sum_trim: the fraction of A-values that should be trimmed from each direction (top and bottom Y%).
        :type sum_trim:  float between 0 and 0.5 (default=0.05)
        :param a_cutoff: a lower bound on the A-values that should be included in the trimmed mean. \
        If set to None, no lower bound will be used.
        :type a_cutoff: float or None (default = -1e10)
        :param ref_column: the column to be used as reference for normalization. If 'auto', \
        then the reference column will be chosen automatically to be \
        the column whose upper quartile is closest to the mean upper quartile.
        :type ref_column: name of a column or 'auto' (default='auto')
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_tmm()
           Normalized 22 features. Normalized inplace.
        """
        assert 0 <= log_ratio_trim < 0.5, "'log_ratio_trim' must be a value in the range 0 <= log_ratio_trim < 5"
        assert 0 <= sum_trim < 0.5, "'sum_trim' must be a value in the range 0 <= sum_trim < 5"

        suffix = '_normTMM'

        if isinstance(ref_column, str) and ref_column.lower() == 'auto':
            upper_quartiles = self.df[self._numeric_columns].quantile(0.75, axis=0)
            ref_index = np.argmin(abs(upper_quartiles - upper_quartiles.mean()))
            ref_column = self.columns[ref_index]

        columns = self._numeric_columns
        m_data, a_data = self._calculate_ma(ref_column, columns)
        scaling_factors = {}
        norm_data = self.df / self.df.sum(axis=0)
        weights = (1 - norm_data).divide(self.df)
        weights = 1 / (weights.add(weights[ref_column], axis=0))
        for i, col in enumerate(columns):
            a_post_cutoff = a_data[i][a_data[i] > a_cutoff] if a_cutoff is not None else a_data[i]
            a_post_cutoff = a_post_cutoff.dropna()
            this_m = m_data[i].dropna()

            m_trim_number = int(np.ceil(log_ratio_trim * m_data[i].notna().sum()))
            a_trim_number = int(np.ceil(sum_trim * a_post_cutoff.notna().sum()))
            trimmed_m = this_m.sort_values().iloc[m_trim_number:this_m.shape[0] - m_trim_number]
            trimmed_a = a_post_cutoff.sort_values().iloc[a_trim_number:a_post_cutoff.shape[0] - a_trim_number]

            zero_genes = self.df[col][self.df[col] == 0].index
            genes_post_trimming = trimmed_m.index.intersection(trimmed_a.index).difference(zero_genes)
            tmm = np.average(trimmed_m.loc[genes_post_trimming], weights=weights.loc[genes_post_trimming, col])
            scaling_factors[col] = tmm

        scaling_factors = pd.Series(scaling_factors)
        # adjust scaling factors to multiply, for symmetry, to 1
        scaling_factors -= scaling_factors.mean()
        scaling_factors = 2 ** scaling_factors
        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

    
    def normalize_rle(self, inplace: bool = True, return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'Relative Log Expression' (RLE) method \
        `(Anders and Huber 2010) <https://doi.org/10.1186/gb-2010-11-10-r106>`_. \
        This is the default normalization method used by R's DESeq2. \
        To calculate the Relative Log Expression scaling factors, you first generate a pseudo-sample by calculating \
        the geometric mean expression of each gene across samples. You then calculate the gene-wise ratio \
        of expression between each sample and the pseudo-sample. You then pick the median ratio within each \
        sample as the scaling factor for that sample.

        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_rle()
           Normalized 22 features. Normalized inplace.
        """
        suffix = '_normRLE'
        data = self.df[self._numeric_columns].dropna(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'):
            pseudo_sample = pd.Series(gmean(data, axis=1), index=data.index)
            ratios = data.divide(pseudo_sample, axis=0)
            scaling_factors = ratios.dropna(axis=0).median(axis=0)

        # TODO: implement a 'control genes' parameter that calculates ratios only for the given control genes

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

   
    def normalize_median_of_ratios(self, sample_grouping: param_typing.GroupedColumns,
                                   reference_group: NonNegativeInt = 0, inplace: bool = True,
                                   return_scaling_factors: bool = False):
        """
        Normalizes the count matrix using the 'Median of Ratios Normalization' (MRN) method \
        `(Maza et al 2013) <https://doi.org/10.4161%2Fcib.25849>`_. \
        This normalization method uses information about the experimental condition of each sample. \
        To calculate the Median of Ratios scaling factors, you first calculate the weighted mean expression of \
        each gene within the replicates of each experimental condition. You then calculate per gene the ratio between \
        each weighted mean in the experimental condition and those of the reference condition. \
        You then pick the median ratio for each experimental condition, and calculate the scaling factor for each \
        sample by multiplying it with the sample's total number of reads. \
        Finally, the scaling factors are adjusted, for symmetry, so that they multiply to 1.

        :type sample_grouping: nested list of column names
        :param sample_grouping: grouping of the samples into conditions. \
        Each grouping should containg all replicates of the same condition.
        :type reference_group: int (default=0)
        :param reference_group: the index of the sample group to be used as the reference condition. \
        Must be an integer between 0 and the number of sample groups -1.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :param return_scaling_factors: if True, return a DataFrame containing the calculated scaling factors.
        :type return_scaling_factors: bool (default=False)
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_median_of_ratios([['cond1','cond2'],['cond3','cond4']])
           Normalized 22 features. Normalized inplace.
        """
        flat_grouping = parsing.flatten(sample_grouping)
        assert len(flat_grouping) >= len(self._numeric_columns), f"'sample_grouping' must include all columns. " \
                                                                 f"Only {len(flat_grouping)} out of " \
                                                                 f"{len(self._numeric_columns)} " \
                                                                 f"numeric columns were included. "
        assert isinstance(reference_group, int) and reference_group >= 0, \
            f"Invalid value for 'reference_group': {reference_group}"
        assert reference_group < len(sample_grouping), f"'reference_group' value {reference_group} " \
                                                       f"is larger than the number of sample groups!"

        suffix = '_normMRN'
        data = self.df[self._numeric_columns]
        weighed_expression = data / data.sum(axis=0)
        weighted_means = []
        for grp in sample_grouping:
            weighted_means.append(weighed_expression.loc[:, grp].mean(axis=1))
        scaling_factors = {}
        for i, grp in enumerate(sample_grouping):
            ratios = weighted_means[i] / weighted_means[reference_group]
            ratios[ratios == np.inf] = np.nan
            median_of_ratios = ratios.median()
            for cond in grp:
                scaling_factors[cond] = median_of_ratios * self.df[cond].sum()
        scaling_factors = pd.Series(scaling_factors)

        # adjust scaling factors to multiply, for symmetry, to 1
        scaling_factors = scaling_factors / gmean(scaling_factors)

        new_df = self._norm_scaling_factors(scaling_factors)
        if return_scaling_factors:
            return [self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True), scaling_factors]
        return self._inplace(new_df, False, inplace, suffix, 'normalize', _is_normalized=True)

 
    def normalize_with_scaling_factors(self, scaling_factor_fname: Union[str, Path], inplace: bool = True):

        """
        Normalizes the reads in the CountFilter using pre-calculated scaling factors. \
        Receives a table of sample names and their corresponding scaling factors, \
        and divides each column in the CountFilter by the corresponding scaling factor.

        :type scaling_factor_fname: str or pathlib.Path
        :param scaling_factor_fname: the .csv file which contains scaling factors for the different libraries.
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If inplace is False, returns a new instance of the Filter object.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter("tests/test_files/counted.csv")
            >>> c.normalize_with_scaling_factors("tests/test_files/scaling_factors.csv")
           Normalized 22 features. Normalized inplace.

        """
        suffix = '_normwithscalingfactors'
        if isinstance(scaling_factor_fname, (str, Path)):
            scaling_factors = io.load_table(scaling_factor_fname).squeeze()
        elif isinstance(scaling_factor_fname, pd.DataFrame):
            scaling_factors = scaling_factor_fname.squeeze()
        elif isinstance(scaling_factor_fname, pd.Series):
            scaling_factors = scaling_factor_fname
        else:
            raise TypeError("Invalid type for 'scaling_factor_fname'!")
        new_df = self._norm_scaling_factors(scaling_factors)
        return self._inplace(new_df, opposite=False, inplace=inplace, suffix=suffix, printout_operation='normalize',
                             _is_normalized=True)

   
   
    def filter_low_reads(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):

        """
        Filter out features which are lowly-expressed in all columns, keeping only features with at least 'threshold' \
        reads in at least one column.

        :type threshold: float
        :param threshold: The minimal number of reads (counts, rpm, rpkm, tpm, etc) a feature should have \
        in at least one sample in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_low_reads(5) # remove all rows whose values in all columns are all <5
            Filtered 6 features, leaving 16 of the original 22 features. Filtered inplace.

        """
        validation.validate_threshold(threshold)
        self._validate_is_normalized()

        new_df = self.df.loc[[True if max(vals) > threshold else False for gene, vals in
                              self.df.loc[:, self._numeric_columns].iterrows()]]
        suffix = f"_filt{threshold}reads"
        return self._inplace(new_df, opposite, inplace, suffix)

    

  
    def filter_by_row_sum(self, threshold: float = 5, opposite: bool = False, inplace: bool = True):

        """
        Removes features/rows whose sum is belove 'threshold'.

        :type threshold: float
        :param threshold: The minimal sum a row should have in order not to be filtered out.
        :type opposite: bool
        :param opposite: If True, the output of the filtering will be the OPPOSITE of the specified \
        (instead of filtering out X, the function will filter out anything BUT X). \
        If False (default), the function will filter as expected.
        :type inplace: bool (default=True)
        :param inplace: If True (default), filtering will be applied to the current CountFilter object. If False, \
        the function will return a new CountFilter instance and the current instance will not be affected.
        :return: If 'inplace' is False, returns a new instance of CountFilter.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter('tests/test_files/counted.csv')
            >>> c.filter_by_row_sum(5) # remove all rows whose sum is <5
            Filtered 4 features, leaving 18 of the original 22 features. Filtered inplace.

        """
        validation.validate_threshold(threshold)
        self._validate_is_normalized()

        new_df = self.df.loc[self.df.sum(axis=1) >= threshold]
        suffix = f"_filt{threshold}sum"
        return self._inplace(new_df, opposite, inplace, suffix)


    
    
    @classmethod
    def from_folder(cls, folder_path: str, save_csv: bool = False, fname: str = None, input_format: str = '.txt'
                    ) -> 'CountFilter':
        """
        Iterates over count .txt files in a given folder and combines them into a single CountFilter table. \
        Can also save the count data table and the uncounted data table to .csv files.

        :param folder_path: str or pathlib.Path. Full path of the folder that contains individual htcount .txt files.
        :param save_csv: bool. If True, the joint DataFrame of count data and uncounted data will be saved \
        to two separate .csv files. The files will be saved in 'folder_path', and named according to the parameters \
        'counted_fname' for the count data, and 'uncounted_fname' for the uncounted data (unaligned, \
        alignment not unique, etc).
        :param fname: str. Name under which to save the combined count data table. Does not need to include \
        the '.csv' suffix.
        :param input_format: the file format of the input files. Default is '.txt'.
        :return: an CountFilter object containing the combined count data from all individual htcount .txt files in the \
        specified folder.


        :Examples:
            >>> from rnalysis import filtering
            >>> c = filtering.CountFilter.from_folder('tests/test_files/test_count_from_folder')
            """
        file_suffix = '.csv'
        if save_csv:
            assert isinstance(fname, str)

            if not fname.endswith(file_suffix):
                fname += file_suffix

            counted_fname = os.path.join(folder_path, fname)

        folder = Path(folder_path)
        counts = pd.DataFrame()
        for item in sorted(folder.iterdir()):
            if item.is_file() and item.suffix == input_format:
                counts = pd.concat([counts, pd.read_csv(item, sep='\t', index_col=0, names=[item.stem])], axis=1)
        assert not counts.empty, f"No valid files with the suffix '{input_format}' were found in '{folder_path}'."

        if save_csv:
            io.save_table(df=counts, filename=counted_fname)

        fname = counted_fname if save_csv else os.path.join(folder.absolute(), folder.name + file_suffix)
        count_filter_obj = cls.from_dataframe(counts, Path(fname))
        return count_filter_obj

    @classmethod
    def from_folder_htseqcount(cls, folder_path: str, norm_to_rpm: bool = False, save_csv: bool = False,
                               counted_fname: str = None, uncounted_fname: str = None,
                               input_format: str = '.txt') -> 'CountFilter':

        """
            Iterates over HTSeq count .txt files in a given folder and combines them into a single CountFilter table. \
            Can also save the count data table and the uncounted data table to .csv files, and normalize the CountFilter \
            table to reads per million (RPM). Note that the saved data will always be count data, and not normalized data, \
            regardless if the CountFilter table was normalized or not.

            :param folder_path: str or pathlib.Path. Full path of the folder that contains individual htcount .txt files.
            :param norm_to_rpm: bool. If True, the CountFilter table will be automatically normalized to reads per \
            million (RPM). If False (defualt), the CountFilter object will not be normalized, and will instead contain \
            absolute count data (as in the original htcount .txt files). \
            Note that if save_csv is True, the saved .csv fill will contain ABSOLUTE COUNT DATA, as in the original \
            htcount .txt files, and NOT normalized data.
            :param save_csv: bool. If True, the joint DataFrame of count data and uncounted data will be saved \
            to two separate .csv files. The files will be saved in 'folder_path', and named according to the parameters \
            'counted_fname' for the count data, and 'uncounted_fname' for the uncounted data (unaligned, \
            alignment not unique, etc).
            :param counted_fname: str. Name under which to save the combined count data table. Does not need to include \
            the '.csv' suffix.
            :param uncounted_fname: counted_fname: str. Name under which to save the combined uncounted data. \
            Does not need to include the '.csv' suffix.
            :param input_format: the file format of the input files. Default is '.txt'.
            :return: an CountFilter object containing the combined count data from all individual htcount .txt files in the \
            specified folder.


            :Examples:
                >>> from rnalysis import filtering
                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder')

                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', norm_to_rpm=True) # This will also normalize the CountFilter to reads-per-million (RPM).
               Normalized 10 features. Normalized inplace.

                >>> c = filtering.CountFilter.from_folder_htseqcount('tests/test_files/test_count_from_folder', save_csv=True, counted_fname='name_for_reads_csv_file', uncounted_fname='name_for_uncounted_reads_csv_file') # This will also save the counted reads and uncounted reads as separate .csv files

            """
        file_suffix = '.csv'
        if save_csv:
            assert isinstance(counted_fname, str)
            assert isinstance(uncounted_fname, str)

            if not counted_fname.endswith(file_suffix):
                counted_fname += file_suffix
            if not uncounted_fname.endswith(file_suffix):
                uncounted_fname += file_suffix

            counted_fname = os.path.join(folder_path, counted_fname)
            uncounted_fname = os.path.join(folder_path, uncounted_fname)

        folder = Path(folder_path)
        df = pd.DataFrame()
        for item in sorted(folder.iterdir()):
            if item.is_file() and item.suffix == input_format:
                df = pd.concat([df, pd.read_csv(item, sep='\t', index_col=0, names=[item.stem])], axis=1)
        assert not df.empty, f"Error: no valid files with the suffix '{input_format}' were found in '{folder_path}'."

        uncounted = df.loc[
            ['__no_feature', '__ambiguous', '__alignment_not_unique', '__too_low_aQual', '__not_aligned']]
        counts = df.drop(uncounted.index, inplace=False)

        if save_csv:
            io.save_table(df=counts, filename=counted_fname)
            io.save_table(df=uncounted, filename=uncounted_fname)

        fname = counted_fname if save_csv else os.path.join(folder.absolute(), folder.name + file_suffix)
        count_filter_obj = cls.from_dataframe(counts, Path(fname))
        if norm_to_rpm:
            count_filter_obj.normalize_to_rpm_htseqcount(uncounted)
        return count_filter_obj

