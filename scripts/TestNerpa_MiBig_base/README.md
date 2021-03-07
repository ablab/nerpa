# Test Nerpa on MiBig dataset

````
https://github.com/ablab/nerpa/tree/master/scripts/TestNerpa_MiBig_base
````


* Create result directory
````
mkdir result
````

* Add print_structure and NERPA to PATH. Let's Nerpa is compiled to ``/Bmo/kolga/soft/Nerpa/bin``
````
PATH=$PATH:/Bmo/kolga/soft/Nerpa/bin
````

* Setup path to testing data and result dir
````
TEST_MIBIG_DATA_PATH=/Bmo/kolga/data/Nerpa/MiBigTest/NRP/
TEST_MIBIG_RESULT_PATH=/Bmo/kolga/run/Nerpa/MiBigResult/
````

##RUN
* Recompile Nerpa
````
PREFIX=/Bmo/kolga/soft/Nerpa ${PATH_TO_NERPA_GIT}/install.sh
````

* Run tests
````
./test_nerpa.sh
````

##Example of script for running tests
````
/Bmo/kolga/runs/Nerpa/MibigResult/run_test_mibig.sh
````
