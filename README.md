# Nerpa

Nerpa is a tool for linking biosynthetic gene clusters (BGCs) to known nonribosomal peptides (NRPs).
BGCs are predicted in genome sequences (FASTA or GBK) with [antiSMASH](https://antismash.secondarymetabolites.org/). 
Known NRPs are accepted in the [SMILES format](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
 and processed with [rBAN](https://web.expasy.org/rban/).

## Dependencies

### Required
* 64-bit Linux system or macOS
* g++ v.5.2+ or clang v.3.2+
* cmake v.3.5+
* Python v.3.6+
* Python libraries:
    * [RDKit](https://www.rdkit.org/docs/Install.html)
    * [networkx](https://networkx.org/documentation/stable/install.html)

### Conditionally required
The following dependencies are required only for specific input data types.   

* java (only if NRPs are provided as SMILES)

## Installation

Get source code from GitHub:

    git clone https://github.com/ablab/nerpa.git
    cd nerpa

and build it with the following script:

    ./install.sh

Nerpa will be built in the directory ./bin. If you wish to install
the tool into another directory, you can specify destination folder
by running the following command in bash or sh:

    PREFIX=<destination_dir> ./install.sh

for example:

    PREFIX=/usr/local ./install.sh

which will install Nerpa into `/usr/local/bin/`.

Note: you should use absolute path for `<destination_dir>`.

After installation you should get `NRPsMatcher` and `nerpa.py`
files in `./bin` (or `<destination_dir>/bin` if you specified PREFIX)
directory.


## Running
### Minimal working example

```
python3 bin/nerpa.py -a test_data/NCBI_subset/genome_predictions/ --structures test_data/NCBI_subset/structure.info.monomers`
less nerpa_results/latest/report.csv
```


### Input
Nerpa takes as input BGC predictions generated with [antiSMASH](https://antismash.secondarymetabolites.org/) 
and NRP structures in the SMILES format.

#### BGC predictions

The most convenient way to get antiSMASH predictions of BGC in your genomic data with antiSMASH is to upload your 
FASTA or GBK file to their [webserver](https://antismash.secondarymetabolites.org/). 
When the server job is completed, you may download archive with results, unpack it and 
provide the path to the unpacked directory or just the main JSON file from it to Nerpa via option `-a`. 

You can also use [the command-line version](https://docs.antismash.secondarymetabolites.org/install/) of antiSMASH.
Nerpa supports outputs from v.3 and v.5. The recommended running parameters for v.5.1.1 are:  
`python run_antismash.py <your_genome.fasta> --output-dir <your_output_dir> --genefinding-tool prodigal --skip-zip-file --enable-nrps-pks --minimal`    

Note that you may specify an unlimited number of antiSMASH output files by using `-a` multiple times or by specifying a root directory with many inputs inside.
You may also write paths to all antiSMASH outputs in a single file and provide it via option `--antismash_output_list`.

#### NRP structures

NRP molecules should be specified in the [SMILES format](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system).
One option is to provide them as a space-separated list of SMILES strings via option `--smiles`.
Another way is to write all structures in a multi-column file and specify it via `--smiles-tsv`. 
Default column separator (`\t`), names of the SMILES column (`SMILES`) and the column with molecule IDs (*row index*) could be adjusted via option `--sep`, `--col_smiles`, and `--col_id`, respectively.


### Command line options

To get the full list of Nerpa options type

`python3 bin/nerpa.py --help`


### Output

Nerpa stores all output files in the user-defined directory if specified via `-o` (should not exist before the run) or 
in `nerpa_results/results_YYYY_MM_DD_HH_MM_SS` otherwise. In the latter case, a symlink `nerpa_results/latest` pointing to the output directory is also created.

The key files/directories inside the output directory are:
* `reports.csv` matched NRP-BGC pairs with scores
* `details_mols` directory with detailed descriptions and exact alignments for each match. 
Each file in the directory corresponds to an NRP and contains information on all matches with this NRP.

### More running examples

```
./nerpa.py -a test_data/MIBiG_subset/genome_predictions --smiles-tsv test_data/MIBiG_subset/structures_info.tsv
./nerpa.py -a test_data/NCBI_subset/genome_predictions/ --structures test_data/NCBI_subset/structure.info.monomers -o NCBI_v3_out_dir
./nerpa.py -a test_data/NCBI_subset/genome_predictions_v5/ --structures test_data/NCBI_subset/structure.info.monomers -o NCBI_v5_out_dir
```