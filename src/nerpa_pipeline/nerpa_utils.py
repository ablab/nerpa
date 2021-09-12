import os
import re
import shutil
import datetime
import shlex
import subprocess

import nerpa_config


def set_up_output_dir(output_dirpath, crash_if_exists, log):
    make_latest_symlink = False

    if not output_dirpath:  # 'output dir was not specified with -o option'
        output_dirpath = os.path.join(os.path.abspath(
            nerpa_config.default_results_root_dirname),
            nerpa_config.default_results_dirname_prefix +
            datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S'))
        make_latest_symlink = True

        # in case of starting two jobs in the same second
        if os.path.isdir(output_dirpath):
            i = 2
            base_dirpath = output_dirpath
            while os.path.isdir(output_dirpath):
                output_dirpath = str(base_dirpath) + '__' + str(i)
                i += 1

    if os.path.isdir(output_dirpath) and crash_if_exists:
        log.error(f"output directory ({output_dirpath}) already exists! "
                  f"Rerun with --force-existing-outdir if you still want to use it as the output dir "
                  f"OR specify another (nonexistent) directory. Exiting now..", to_stderr=True)
    else:
        if os.path.isdir(output_dirpath):  # TODO: check whether we want to completely remove it always
            shutil.rmtree(output_dirpath)
        os.makedirs(output_dirpath)

    # 'latest' symlink
    if make_latest_symlink:
        prev_dirpath = os.getcwd()
        os.chdir(nerpa_config.default_results_root_dirname)

        latest_symlink = 'latest'
        if os.path.islink(latest_symlink):
            os.remove(latest_symlink)
        os.symlink(os.path.basename(output_dirpath), latest_symlink)

        os.chdir(prev_dirpath)

    return os.path.abspath(output_dirpath)


def sys_call(cmd, log, indent='  ', cwd=None, verbose=True):
    def _process_readline(line):
        return str(line, "utf-8").rstrip()

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    if verbose:
        log.info("\n== Running: %s\n" % (' '.join(cmd_list)))
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    # output = ""
    while not proc.poll():
        line = _process_readline(proc.stdout.readline())
        if line:
            log.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = _process_readline(line)
        if line:
            log.info(indent + line)
            # if log:
            #     log.info(line)
            # else:
            #     output += line + "\n"

    if proc.returncode:
        log.error("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))

    if verbose:
        log.info("\n== Done\n")
    # return output


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    cur_path = os.path.dirname(os.path.abspath(__file__))
    if is_exe(os.path.join(cur_path, fname)):
        return os.path.join(cur_path, fname)

    if fpath:
        if is_exe(program):
            return program
    elif "PATH" in os.environ:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def get_path_to_program(program, dirpath=None, min_version=None):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
            if not min_version or check_version(fpath, min_version):
                return True

    def check_version(fpath, min_version):
        p = subprocess.Popen([fpath, '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()

        version_pattern = re.compile('(?P<major_version>\d+)\.(?P<minor_version>\d+)')
        searchstring = stdout.decode('utf8').strip()

        # ad hoc workaround to AS 5.2.0 printing FutureWarning to stdout
        searchstring = searchstring.split('\n')[-1]

        v = version_pattern.search(searchstring)
        if not v.group('major_version') or not v.group('minor_version'):
            return False
        version, minor_version = map(int, min_version.split('.'))
        if int(v.group('major_version')) == version and int(v.group('minor_version')) >= minor_version:
            return True

    if dirpath:
        exe_file = os.path.join(dirpath, program)
        if is_exe(exe_file):
            return exe_file
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None