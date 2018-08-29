# NRP Matcher

NRP Matcher is a tool which links gene cluster to known natural products.

## Installation
### Requirements
   *  https://github.com/ablab/dereplicator/ must be added to PATH. Folders
   "Fragmentation_rule" and "configs" must be in one of the paths:
   "\<DEREPLICATOR INSTALL DIR\>/", \<DEREPLICATOR INSTALL DIR\>/../", \<DEREPLICATOR INSTALL DIR\>/../../",
   \<DEREPLICATOR INSTALL DIR\>/../share/", \<DEREPLICATOR INSTALL DIR\>/../share/npdtools/"


To install into CWD: ./install.sh

To install into specific dir: PREFIX=<destination_dir> ./install.sh

Note: you should use absolute path for <destination_dir>! Main executables will be placed under <destination_dir>/bin/

## Running
### Input
### Command line options
./run_nrp_matcher.py --predictions \<path to file with paths to ctg1_nrpspredictor2_codes.txt files\> --lib_info \<path to file with paths to mol files\>
### Output
### Example