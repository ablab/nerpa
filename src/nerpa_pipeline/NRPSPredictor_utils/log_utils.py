import sys


def error(msg, exit=True):
    print('ERROR! ' + msg)
    if exit:
        sys.exit(1)


def info(msg, verbose=True):
    if verbose:
        print(msg)
