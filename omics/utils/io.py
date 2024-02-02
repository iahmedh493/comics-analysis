import asyncio
import concurrent.futures
import contextlib
import ftplib
import functools
import gzip
import hashlib
import inspect
import json
import os
import queue
import random
import re
import shlex
import shutil
import subprocess
import threading
import time
import typing
import warnings
from datetime import date, datetime
from functools import lru_cache
from io import StringIO
from itertools import chain
from pathlib import Path
from sys import executable
from typing import List, Set, Union, Iterable, Tuple, Dict, Any, Callable, Literal
from urllib.parse import urlparse, parse_qs, urlencode

import aiohttp
import aiolimiter
import appdirs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import tenacity
import yaml
from defusedxml import ElementTree
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

__version__ = "3.10.1"
from . import parsing, validation, __path__


def get_gui_cache_dir() -> Path:
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    return cache_dir.joinpath('rnalysis_gui')


def get_data_dir() -> Path:
    data_dir = Path("/Users/ibrahimahmed/projects/GUI/result_dir")
    #data_dir = Path(appdirs.user_data_dir('RNAlysis', roaming=True))
    return data_dir


def get_tutorial_videos_dir() -> Path:
    data_dir = get_data_dir()
    return data_dir.joinpath('videos')


def get_todays_cache_dir() -> Path: 
    today = date.today().strftime('%Y_%m_%d')
    cache_dir = Path("/Users/ibrahimahmed/projects/GUI/result_dir")
    #cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    todays_dir = cache_dir.joinpath(today)
    todays_dir.mkdir(parents=True, exist_ok=True)
    return todays_dir


def load_cached_file(filename: str):
    directory = get_todays_cache_dir()
    file_path = directory.joinpath(filename)
    if file_path.exists():
        with open(file_path) as f:
            return f.read()
    else:
        return None


def cache_file(content: str, filename: str):
    directory = get_todays_cache_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    file_path = directory.joinpath(filename)
    with open(file_path, 'w') as f:
        f.write(content)


def clear_directory(directory: Union[str, Path]):
    directory = Path(directory)
    if not directory.exists():
        return

    for item in directory.iterdir():
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item, ignore_errors=True)


def clear_cache():
    cache_dir = Path(appdirs.user_cache_dir('RNAlysis'))
    clear_directory(cache_dir)


def clear_gui_cache():
    directory = get_gui_cache_dir()
    clear_directory(directory)


def load_cached_gui_file(filename: Union[str, Path], load_as_obj: bool = True) -> Union[
    str, set, pd.DataFrame, bytes, None]:
    """
    Load a cached file from the GUI cache directory.

    :param filename: The name of the file to load.
    :type filename: str
    :param load_as_obj: Whether to load the file as an object or raw content as string. Defaults to True.
    :type load_as_obj: bool
    :return: The contents of the file, loaded as either a string, a Pandas DataFrame or a set.
    :rtype: str or pandas.DataFrame or set or Bytes or None
    """
    directory = get_gui_cache_dir()
    file_path = directory.joinpath(filename)
    if file_path.exists():
        if file_path.suffix in {'.csv', '.tsv', '.parquet'} and load_as_obj:
            return load_table(file_path, index_col=0)
        elif file_path.suffix in {'.txt'} and load_as_obj:
            with open(file_path) as f:
                return {item.strip() for item in f.readlines()}

        else:
            with open(file_path, 'rb') as f:
                return f.read()
    else:
        return None


def cache_gui_file(item: Union[pd.DataFrame, set, str], filename: str):
    directory = get_gui_cache_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    file_path = directory.joinpath(filename)
    if isinstance(item, (pd.DataFrame, pd.Series)):
        save_table(item, file_path, index=True)
    elif isinstance(item, set):
        save_gene_set(item, file_path)
    elif isinstance(item, str):
        with open(file_path, 'w') as f:
            f.write(item)
    elif isinstance(item, plt.Figure):
        item.savefig(file_path)
    else:
        raise TypeError(type(item))


def check_changed_version():  # pragma: no cover
    data_dir = get_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    filename = data_dir.joinpath('latest_version.txt')
    if not filename.exists():
        ver = ''
    else:
        with open(filename) as f:
            ver = f.read()
    current_ver = __version__
    # update latest version to current version
    with open(filename, 'w') as f:
        f.write(current_ver)

    return ver != current_ver


def save_gui_session(session_filename: Union[str, Path], file_names: List[str], item_names: List[str], item_types: list,
                     item_properties: list, pipeline_names: List[str], pipeline_files: List[str]):
    session_filename = Path(session_filename)
    session_folder = session_filename
    if session_folder.exists():
        if session_folder.is_dir():
            shutil.rmtree(session_folder)
        else:
            session_folder.unlink()
    session_folder.mkdir(parents=True)

    session_data = dict(files=dict(), pipelines=dict(), metadata=dict())
    for file_name, item_name, item_type, item_property in zip(file_names, item_names, item_types, item_properties):
        shutil.move(Path(get_gui_cache_dir().joinpath(file_name)), session_folder.joinpath(file_name))
        session_data['files'][file_name] = (item_name, item_type.__name__, item_property)

    for i, (pipeline_name, pipeline_file) in enumerate(zip(pipeline_names, pipeline_files)):
        pipeline_filename = session_folder.joinpath(f"pipeline_{i}.yaml")
        Path(pipeline_filename).write_text(pipeline_file)
        session_data['pipelines'][pipeline_filename.name] = pipeline_name

    session_data['metadata']['creation_time'] = get_datetime()
    session_data['metadata']['name'] = Path(session_filename).stem
    session_data['metadata']['n_tabs'] = len(session_data['files'])
    session_data['metadata']['n_pipelines'] = len(session_data['pipelines'])
    session_data['metadata']['tab_order'] = file_names

    with open(session_folder.joinpath('session_data.yaml'), 'w') as f:
        yaml.safe_dump(session_data, f)
    shutil.make_archive(session_folder.with_suffix(''), 'zip', session_folder)
    shutil.rmtree(session_folder)
    session_filename.with_suffix('.zip').replace(session_filename.with_suffix('.rnal'))


def load_gui_session(session_filename: Union[str, Path]):
    session_filename = Path(session_filename)
    try:
        session_filename.with_suffix('.rnal').rename(session_filename.with_suffix('.rnal.zip'))
        shutil.unpack_archive(session_filename.with_suffix('.rnal.zip'),
                              get_gui_cache_dir().joinpath(session_filename.name))
    finally:
        session_filename.with_suffix('.rnal.zip').rename(session_filename.with_suffix('.rnal'))

    session_dir = get_gui_cache_dir().joinpath(session_filename.name)
    assert session_dir.exists()

    items = []
    item_names = []
    item_types = []
    item_properties = []
    pipeline_files = []
    pipeline_names = []
    with open(session_dir.joinpath('session_data.yaml')) as f:
        session_data = yaml.safe_load(f)
    if 'tab_order' in session_data['metadata'] and len(session_data['metadata']['tab_order']) == len(
        session_data['files']):
        filenames = session_data['metadata']['tab_order']
    else:
        filenames = session_data['files'].keys()

    for file_name in filenames:
        file_path = session_dir.joinpath(file_name)
        assert file_path.exists() and file_path.is_file()
        item = load_cached_gui_file(Path(session_filename.name).joinpath(file_name))
        item_name, item_type, item_property = session_data['files'][file_name]
        items.append(item)
        item_names.append(item_name)
        item_types.append(item_type)
        item_properties.append(item_property)

    for pipeline_filename in session_data['pipelines'].keys():
        pipeline_path = session_dir.joinpath(pipeline_filename)
        assert pipeline_path.exists() and pipeline_path.is_file()
        pipeline_files.append(pipeline_path.read_text())
        pipeline_names.append(session_data['pipelines'][pipeline_filename])

    shutil.rmtree(session_dir)
    return items, item_names, item_types, item_properties, pipeline_names, pipeline_files


def load_table(filename: Union[str, Path], index_col: int = None, drop_columns: Union[str, List[str]] = False,
               squeeze=False, comment: str = None, engine: Literal['pyarrow', 'auto'] = 'auto'):
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
    assert filename.suffix.lower() in {'.csv', '.tsv', '.txt', '.parquet'}, \
        f"RNAlysis cannot load files of type '{filename.suffix}'. " \
        f"Please convert your file to a .csv, .tsv, .txt, or .parquet file and try again."

    if filename.suffix.lower() == '.parquet':
        df = pd.read_parquet(filename, engine=engine)
    else:
        kwargs = dict(sep=None, engine='python' if engine == 'auto' else engine, encoding='ISO-8859-1', comment=comment,
                      skipinitialspace=True)
        if index_col is not None:
            kwargs['index_col'] = index_col
        df = pd.read_csv(filename, **kwargs)

    if squeeze:
        df = df.squeeze("columns")

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
        drop_columns_lst = parsing.data_to_list(drop_columns)
        assert validation.isinstanceiter(drop_columns_lst,
                                         str), f"'drop_columns' must be str, list of str, or False; " \
                                               f"is instead {type(drop_columns)}."
        for col in drop_columns_lst:
            col_stripped = col.strip()
            if col_stripped in df:
                df.drop(col_stripped, axis=1, inplace=True)
            else:
                raise IndexError(f"The argument {col} in 'drop_columns' is not a column in the loaded csv file!")
    return df


def save_table(df: pd.DataFrame, filename: Union[str, Path], postfix: str = None, index: bool = True):
    """
    save a pandas DataFrame to csv/parquet file.

    :param df: pandas DataFrame to be saved
    :param filename: a string or pathlib.Path object stating the original name of the file
    :type postfix: str, default None
    :param postfix: A postfix to be added to the original name of the file. If None, no postfix will be added.
    :param index: if True, saves the DataFrame with the indices. If false, ignores the index.
    """
    fname = Path(filename)
    if postfix is None:
        postfix = ''
    else:
        assert isinstance(postfix, str), "'postfix' must be either str or None!"
    new_fname = os.path.join(fname.parent.absolute(), f"{fname.stem}{postfix}{fname.suffix}")
    if fname.suffix.lower() == '.parquet':
        if isinstance(df, pd.Series):
            df = df.to_frame()
        df.to_parquet(new_fname, index=index)
    else:
        df.to_csv(new_fname, header=True, index=index)



def run_r_script(script_path: Union[str, Path], r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    #if r_installation_folder == 'auto':
        #prefix = "Rscript"
    #else:
        #prefix = f'{Path(r_installation_folder).as_posix()}/bin/Rscript'
    prefix = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
    script_path = Path(script_path).as_posix()
    assert Path(script_path).exists() and Path(
        script_path).is_file(), f"Could not find the requested R script: {script_path}"

    try:
        return_code, _ = run_subprocess([prefix, "--help"], False, False)
        if return_code:
            raise FileNotFoundError

    except FileNotFoundError:
        raise FileNotFoundError("Failed to find R executable. "
                                "Please make sure your R installation folder is correct. \n"
                                "(For example: 'C:/Program Files/R/R-4.2.3')")

    return_code, stderr = run_subprocess([prefix, script_path])
    if return_code:
        full_err = 'See R log below. \n' + '\n'.join([s.rstrip() for s in stderr])
        short_err = stderr[-1].rstrip()
        for i, s in enumerate(stderr):
            if s.startswith('Error'):
                short_err = s.rstrip() + stderr[i + 1].rstrip()
                break
        raise ChildProcessError(f"R script failed to execute: '{short_err}'. See full error report below.") \
            from RuntimeError(full_err)


def stdout_reader(pipe, log_filename, lock, print_output: bool = True):
    with open(log_filename, 'a') if log_filename is not None else contextlib.nullcontext() as logfile:
        for line in (pipe if isinstance(pipe, list) else iter(pipe.readline, b'')):
            decoded_line = line.decode('utf8', errors="ignore")

            if print_output:
                print(decoded_line)

            if log_filename is not None:
                with lock:
                    logfile.write(decoded_line)
    if not isinstance(pipe, list):
        pipe.close()


def stderr_reader(pipe, stderr_record, log_filename, lock, print_output: bool = True):
    with open(log_filename, 'a') if log_filename is not None else contextlib.nullcontext() as logfile:
        for line in (pipe if isinstance(pipe, list) else iter(pipe.readline, b'')):
            decoded_line = line.decode('utf8', errors="ignore")
            stderr_record.append(decoded_line)

            if print_output:
                print(decoded_line)

            if log_filename is not None:
                with lock:
                    logfile.write(decoded_line)
    if not isinstance(pipe, list):
        pipe.close()


def run_subprocess(args: List[str], print_stdout: bool = True, print_stderr: bool = True,
                   log_filename: Union[str, None] = None, shell: bool = False) -> Tuple[int, List[str]]:
    # join List of args into a string of args when running in shell mode
    if shell:
        try:
            args = shlex.join(args)
        except AttributeError:
            args = ' '.join([shlex.quote(arg) for arg in args])

    stderr_record = []
    lock = threading.Lock()
    stdout = subprocess.PIPE
    stderr = subprocess.PIPE

    process = subprocess.Popen(args, stdout=stdout, stderr=stderr, shell=shell)

    stdout_thread = threading.Thread(target=stdout_reader, args=(process.stdout, log_filename, lock, print_stdout))
    stderr_thread = threading.Thread(target=stderr_reader,
                                     args=(process.stderr, stderr_record, log_filename, lock, print_stderr))

    stdout_thread.start()
    stderr_thread.start()
    stdout_thread.join()
    stderr_thread.join()
    process.wait()

    return process.returncode, stderr_record











