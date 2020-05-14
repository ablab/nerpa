# Test Nerpa on MiBig dataset
## Preparation
* Install RDkit ``https://www.rdkit.org/docs/Install.html``
* Copy folder to local directory
````
https://github.com/ablab/nerpa/tree/master/scripts/TestNerpa_MiBig_base
````

For example: 
````
git clone https://github.com/ablab/nerpa.git
cp  -r nerpa/scripts/TestNerpa_MiBig_base /Bmo/kolga/scripts/Nerpa/
cd  /Bmo/kolga/scripts/Nerpa/TestNerpa_MiBig_base
````

* Extract data from data.tar.gz
````
tar -xvf data.tar.gz
````

* Create result directory
````
mkdir result
````

* Add print_structure and NERPA to PATH. Let's Nerpa is compiled to ``/Bmo/kolga/soft/Nerpa/bin``
````
PATH=$PATH:/Bmo/kolga/soft/Nerpa/bin
PATH=$PATH:/Bmo/kolga/soft/dereplicator_build/bin
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