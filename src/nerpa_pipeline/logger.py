import sys
import logging

log = None


def init(log_name):
    global log
    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    setup_log_handler(logging.StreamHandler(sys.stdout))


def add_file_handler(fpath):
    log_handler = logging.FileHandler(fpath, mode='w')
    setup_log_handler(log_handler)


def setup_log_handler(handler, level=logging.DEBUG):
    handler.setFormatter(logging.Formatter('%(asctime)s\t%(message)s', '%Y-%m-%d %H:%M:%S'))
    handler.setLevel(level)
    log.addHandler(handler)


def cleanup():
    if log is not None:
        for handler in list(log.handlers):
            log.removeHandler(handler)


def error(msg):
    # TODO: clean intermediate large files if needed
    # ...

    for handler in list(log.handlers):  # do not print to stdout
        if type(handler) == logging.StreamHandler:
            log.removeHandler(handler)
    setup_log_handler(logging.StreamHandler(sys.stderr), level=logging.ERROR)
    log.error('ERROR: ' + msg)
    log.error("\n\nIn case you have troubles running our tool, "
              "you can post an issue on https://github.com/ablab/nerpa/issues "
              "or write to aleksey.gurevich@spbu.ru\n")
    sys.exit(1)


def exception(e):
    if log is not None and log.handlers:
        log.error('')
        log.exception(e)
    else:
        sys.stderr.write(str(e) + '\n')


def info(msg, silent=False):
    if not silent:
        log.info('INFO: ' + msg)


def warning(msg):
    log.warning('WARNING: ' + msg)
