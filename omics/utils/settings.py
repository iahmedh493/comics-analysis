import os
from pathlib import Path
from typing import Union, Any

import yaml


from . import io


def get_settings_file_path():
    """
    Generates the full path of the 'settings.yaml' file.
    :returns: the path of the settings.yaml file.
    :rtype: pathlib.Path
    """
    directory = io.get_data_dir()
    if not directory.exists():
        directory.mkdir(parents=True)
    return directory.joinpath('settings.yaml')


def reset_settings():
    pth = get_settings_file_path()
    if pth.exists():
        pth.unlink()


def load_settings_file():
    """
    loads and parses the settings.yaml file into a dictionary.
    :rtype: dict
    """
    settings_pth = get_settings_file_path()
    if not settings_pth.exists():
        return dict()
    with settings_pth.open() as f:
        settings = yaml.safe_load(f)
        if settings is None:
            settings = dict()
        return settings


def update_settings_file(value: Any, key: str):
    """
    Receives a key and a value, and updates/adds the key and value to the settings.yaml file.

    :param value: the value to be added/updated (such as Reference Table path)
    :type value: str
    :param key: the key to be added/updated (such as __attr_file_key__)
    :type key: str
    """
    settings_pth = get_settings_file_path()
    out = load_settings_file()
    out[key] = value
    with settings_pth.open('w') as f:
        yaml.safe_dump(out, f)


def is_setting_in_file(key) -> bool:
    """
    Returns True if a settings value is defined for given key, and False otherwise.
    :type key: str
    :param key: the key in the settings file whose value to read.
    :rtype: bool
    """
    settings = load_settings_file()
    return key in settings



def make_temp_copy_of_settings_file():
    """
    Make a temporary copy ('temp_settings.yaml') of the default settings file ('settings.yaml').
    """
    pth = get_settings_file_path()
    try:
        remove_temp_copy_of_settings_file()
    except FileNotFoundError:
        print("no previous temporary settings file existed")
    if not pth.exists():
        print("no previous settings file exists")
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml'), 'w') as tempfile, pth.open() as originfile:
        tempfile.writelines(originfile.readlines())


def set_temp_copy_of_settings_file_as_default():
    """
    Copy the contents of the temporary settings file ('temp_settings.yaml') into the default settings file \
    ('settings.yaml'), if a temporary settings file exists.
    """
    pth = get_settings_file_path()
    if pth.exists():
        pth.unlink()
    if not Path(os.path.join(str(pth.parent), 'temp_settings.yaml')).exists():
        print("no temporary settings file exists")
        return
    with open(os.path.join(str(pth.parent), 'temp_settings.yaml')) as temp_file, pth.open('w') as original_file:
        original_file.writelines(temp_file.readlines())


def remove_temp_copy_of_settings_file():
    """
    Remove the temporary copy of the settings file ('temp_settings.yaml'), if such file exists.
    """
    pth = get_settings_file_path()
    tmp_pth = Path(os.path.join(str(pth.parent), 'temp_settings.yaml'))
    if tmp_pth.exists():
        tmp_pth.unlink()
