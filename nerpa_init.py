#!/usr/bin/env python3

import sys
from os.path import abspath, dirname, realpath, join, isdir

# developers configuration
nerpa_root_dir = abspath(dirname(realpath(__file__)))
bin_dir = join(nerpa_root_dir, "bin")
# common dirs: only basic values are set here, all are adjusted to the real dir paths in init()
python_modules_dir = join(nerpa_root_dir, 'src')
external_tools_dir = None
configs_dir = None


def init():
    global nerpa_root_dir
    global bin_dir
    global python_modules_dir
    global external_tools_dir
    global configs_dir

    # check whether we are in the "user" configuration (bin + share dirs) or in the "developer" one (default setting)
    install_prefix = dirname(nerpa_root_dir)
    possible_bin_dir = join(install_prefix, 'bin')
    possible_nerpa_root_dir = join(install_prefix, 'share', 'nerpa')
    if isdir(possible_bin_dir) and isdir(possible_nerpa_root_dir):  # "user" configuration (with ./bin/ & ./share/nerpa/)
        bin_dir = possible_bin_dir
        nerpa_root_dir = possible_nerpa_root_dir
        python_modules_dir = nerpa_root_dir

    python_modules_dir = join(python_modules_dir, 'nerpa_pipeline')
    external_tools_dir = join(nerpa_root_dir, 'external_tools')
    configs_dir = join(nerpa_root_dir, 'configs')


if __name__ == "__main__":
    nerpa_py_path = join(nerpa_root_dir, "nerpa.py")
    sys.stderr.write("Please use " + nerpa_py_path + " for running Nerpa\n")
