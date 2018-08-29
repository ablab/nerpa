# NRP Matcher

NRP Matcher is a tool which links gene cluster to known natural products.

## Dependencies

NRP Matcher requires a 64-bit Linux system or MAC OS and Python 3,
g++ (version 5.2 or higher) and cmake (version 3.5 or higher) to be
pre-installed on it.

Also https://github.com/ablab/dereplicator/ must be
added to PATH. Folders "Fragmentation_rule" and "configs" must be
located in one of the paths:
* \<DEREPLICATOR INSTALL DIR\>/
* \<DEREPLICATOR INSTALL DIR\>/../
* \<DEREPLICATOR INSTALL DIR\>/../../
* \<DEREPLICATOR INSTALL DIR\>/../share/
* \<DEREPLICATOR INSTALL DIR\>/../share/npdtools/


## Installation

To compile NRP Matcher you can download the NRP Matcher source code:

    git clone https://github.com/olga24912/NRPsMatcher.git
    cd NRPsMatcher

and build it with the following script:

    ./install.sh

NRP Matcher will be built in the directory ./bin. If you wish to install
NRP Matcher into another directory, you can specify full path of destination folder
by running the following command in bash or sh:

    PREFIX=<destination_dir> ./install.sh

for example:

    PREFIX=/usr/local ./install.sh

which will install NRP Matcher into /usr/local/bin.

Note: you should use absolute path for <destination_dir>.

After installation you will get NRPsMatcher and run_nrp_matcher.py
files in ./bin (or <destination_dir>/bin if you specified PREFIX)
directory.

We also suggest adding NRP Matcher installation directory to PATH variable.

## Running
### Input
### Command line options
./run_nrp_matcher.py --predictions \<path to file with paths to ctg1_nrpspredictor2_codes.txt files\> --lib_info \<path to file with paths to mol files\>
### Output
### Example