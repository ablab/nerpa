import sys
import logging
import datetime
import platform
import os


class NerpaLogger(object):
    _NAME = 'nerpa'
    _SINGLE_INDENT = '  '

    _logger = None
    _log_fpath = ''
    _start_time = None
    _num_warnings = 0

    def __init__(self):
        self._logger = logging.getLogger(self._NAME)
        self._logger.setLevel(logging.DEBUG)
        self._set_up_console_handler()

    def enable_debug_mode(self):  # by default: print all DEBUG-level to file but only INFO-level to stdout
        for handler in self._logger.handlers:
            handler.setLevel(logging.DEBUG)

    def set_up_file_handler(self, output_dir):
        self._log_fpath = os.path.join(output_dir, self._NAME + '.log')
        file_handler = logging.FileHandler(self._log_fpath, mode='w')
        file_handler.setLevel(logging.DEBUG)
        self._logger.addHandler(file_handler)

    def _set_up_console_handler(self):
        console_handler = logging.StreamHandler(sys.stdout)
        # console_handler.setFormatter(logging.Formatter(indent_val * '  ' + '%(message)s'))
        # handler.setFormatter(logging.Formatter('%(asctime)s\t%(message)s', '%Y-%m-%d %H:%M:%S'))
        console_handler.setLevel(logging.INFO)
        self._logger.addHandler(console_handler)

    def error(self, msg, to_stderr=False, is_exception=False):
        if msg:
            if is_exception:
                msg = 'EXCEPTION! ' + str(msg)
            else:
                msg = 'ERROR! ' + str(msg)
            msg += "\nIn case you have troubles running our tool, " \
                   "you can post an issue on https://github.com/ablab/nerpa/issues " \
                   "or write to aleksey.gurevich@spbu.ru"

        if to_stderr or self._logger is None or not self._logger.handlers:
            sys.stderr.write('\n' + msg + '\n')
        else:
            self._logger.error('')
            if is_exception:
                self._logger.exception(msg)
            else:
                self._logger.error(msg)

        sys.exit(1)

    def exception(self, e):
        self.error(e, is_exception=True)

    def info(self, msg, indent=0):
        # self._logger.info('INFO: ' + msg)
        self._logger.info(indent * self._SINGLE_INDENT + msg)

    def debug(self, msg):
        self._logger.debug(msg)

    def warning(self, msg):
        self._logger.warning('WARNING: ' + msg)
        self._num_warnings += 1

    def print_command_line(self):
        args = sys.argv[:]
        for i, arg in enumerate(args):
            if ' ' in arg or '\t' in arg:
                args[i] = "'" + arg + "'"
        self.info('')
        self.info('Started with command: ' + ' '.join(args))

    def print_system_info(self):
        self.info('')
        self.info("System information:")
        self.info("OS: " + platform.platform(), indent=1)
        self.info("Python version: " + '.'.join(map(str, sys.version_info[:2])), indent=1)
        try:
            import multiprocessing
            self.info("CPUs number: " + str(multiprocessing.cpu_count()), indent=1)
        except ImportError:
            self.info("Problem occurred when getting CPUs number information", indent=1)

    def print_timestamp(self, msg=''):
        now = datetime.datetime.now()
        current_time = now.strftime("%Y-%m-%d %H:%M:%S")
        self.info('')
        self.info(msg + current_time)
        return now

    def start(self):
        # self.print_tool_version()
        self.print_command_line()
        self.print_system_info()

        self._start_time = self.print_timestamp('Started: ')
        self.info('')
        self.info('Logging to ' + self._log_fpath)

    def finish(self):
        self.info('Log is saved to ' + self._log_fpath, indent=1)

        finish_time = self.print_timestamp('Finished: ')
        self.info('Elapsed time: ' + str(finish_time - self._start_time))
        if self._num_warnings:
            self.info("WARNINGs: %d" % self._num_warnings)

        self.info('\nThank you for using Nerpa!')

        for handler in self._logger.handlers:
            self._logger.removeHandler(handler)
