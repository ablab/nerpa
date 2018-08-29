# NRP Matcher

NRP Matcher is a tool which links gene cluster to known natural products.

## Dependencies

NRP Matcher requires a 64-bit Linux system or MAC OS and Python 3,
g++ (version 5.2 or higher) and cmake (version 3.5 or higher) to be
pre-installed on it.

Also https://github.com/ablab/dereplicator/ must be
added to PATH. Folders "Fragmentation_rule" and "configs" must be
located in one of the paths:
- \<DEREPLICATOR INSTALL DIR\>/
- \<DEREPLICATOR INSTALL DIR\>/../
- \<DEREPLICATOR INSTALL DIR\>/../../
- \<DEREPLICATOR INSTALL DIR\>/../share/
- \<DEREPLICATOR INSTALL DIR\>/../share/npdtools/


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
NRP Matcher takes as input file with list of paths to gene cluster prediction files
and file with paths to files with NRP structure in MOL format.

#### Predictions

**TODO:**  нужен скрипт и бинарник antismash что бы получить нужное предсказание.
Или... долго описывать как нужный файлик получать. Или вообще передавать список геномов
и самостоятельно запускать antismash...

#### NRPs structures

By using command line interface you need specify info file,
where each line described one NRP in following format:

    <path to file with NRP structure in MOL format> <any extra information about NRP>

for example:

    streptomedb/streptomedb.1.mol geranylphenazinediol 348.184 1

You can read about MOL format [here](https://en.wikipedia.org/wiki/Chemical_table_file#Molfile).

If you have NRP structure in some other format we recommend use [molconvert](https://docs.chemaxon.com/display/docs/Molecule+file+conversion+with+Molconverter).
To convert smile string to required MOL file you can run:

    molconvert mol:V3+H --smiles <smile string> -o <nrp file>


Example of info and mols file you can find in

    <installing_dir>/share/library.info.streptomedb
    <installing_dir>/share/streptomedb/

### Command line options

To run NRP Matcher from the command line type

    python3 run_nrp_matcher.py [options]

#### Options

<p>
<code>-h</code> (or <code>--help</code>)    Print help
</p>

<p>
    <code>-p</code> (or <code>--predictions</code>) <code>&lt;file_name></code>  File with paths to prediction files. Required option.
</p>

<p>
    <code>--lib_info &lt;file_name></code>  File with paths to nrp structure description files in MOL format. Required option.
</p>

<p>
    <code>-o</code> (or <code>--local_output_dir</code>) <code>&lt;output_dir></code> Specify the output directory.
</p>

### Output



### Example