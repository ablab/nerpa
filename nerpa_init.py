#!/usr/bin/env python

import os
import sys
from os.path import abspath, dirname, realpath, join, isfile

source_dirs = [""]

# developers configuration
nerpa_home = abspath(dirname(realpath(__file__)))
bin_home = join(nerpa_home, "bin")
python_modules_home = join(nerpa_home, "src")

def init():
    global nerpa_home
    global bin_home
    global python_modules_home

    if isfile(os.path.join(nerpa_home, "NRPsMatcher")):
        install_prefix = dirname(nerpa_home)
        bin_home = join(install_prefix, "bin")
        nerpa_home = join(install_prefix, "share", "nerpa")
        python_modules_home = nerpa_home

    for dir in source_dirs:
        sys.path.append(join(python_modules_home, "nerpa_pipeline", dir))

if __name__ == "__main__":
    nerpa_py_path = join(dirname(realpath(__file__)), "nerpa.py")
    sys.stderr.write("Please use " + nerpa_py_path + " for running Nerpa\n")