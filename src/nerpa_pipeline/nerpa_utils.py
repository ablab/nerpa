import os
import datetime

import nerpa_config


def set_up_output_dir(output_dirpath):
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

    if not os.path.isdir(output_dirpath):
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
