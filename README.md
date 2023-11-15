# First Passage Time on Linear Framework Graphs

A C++ package with Python binding (Pybind11) for First Passage Time calculations on gene regulatory networks using the Linear Framework, partially based on [kmnam/markov-digraphs](https://github.com/kmnam/markov-digraphs.git) from Kee-Myoung Nam. The code has been tested on macOS Sonoma version 14.0, Ubuntu 20.4 and Centos7


# Compile the models

Regardless of the system you are using, the first step is to create a conda environment with the minimal packages requirments for compiling. To install minconda on your system follow these instructions: https://docs.conda.io/projects/miniconda/en/latest/. Once conda is installed, you can proceed with creating the environment:
```
cd ${path_to_fpt_repo}
conda env create -f environment.yml
conda activate fpt_env
```
## Compiling on macOS 
Within the environment double check the presence of pybind11 includes by typing 
```
python -m pybind11 --includes
```
The output should be an absolute path to a python3.11 include file and the pybind11 package include folder within the conda fpt_env folder in your home (double ckeck). If the pybind11 includes are in the conda env you can procede with compiling: 
```
make --file=Makefile_mac dependencies 
```
this command locally installes all the c/c++ dependencies required (boost, eigen, gmp, mpfr, mpc). Please consider doing this step even if you have previously installed these packages, the code is version specific. To test that everything went fine during this step you can compile a series of tests: 
```
make --file=Makefile_mac tests 
```
if this step is succesfull you will see binaries in the test folder. To compile the models run
```
make --file=Makefile_mac ladders 
```
the binaries are authomatically stored in the bin folder of the repository. You can make them available globally or just refer to them in a python script ysing sys.path.insert with the path to FPT/bin. Enjoy the computation :)

## Compiling on Linux
Within the environment double check the presence of pybind11 includes by typing 
```
python -m pybind11 --includes
```
The output should be an absolute path to a python3.11 include file and the pybind11 package include folder within the conda fpt_env folder in your home (double ckeck). If the pybind11 includes are in the conda env you can procede with compiling: 
```
make --file=Makefile_linux dependencies 
```
this command locally installes all the c/c++ dependencies required (boost, eigen, gmp, mpfr, mpc). Please consider doing this step even if you have previously installed these packages, the code is version specific. To test that everything went fine during this step you can compile a series of tests: 
```
make --file=Makefile_linux tests 
```
if this step is succesfull you will see binaries in the test folder. To compile the models run
```
make --file=Makefile_linux ladders 
```
the binaries are authomatically stored in the bin folder of the repository. You can make them available globally or just refer to them in a python script ysing sys.path.insert with the path to FPT/bin. Enjoy the computation :)


# Singularity Container

To build the singularity container run (soon will be uploaded to singularity-hub): 
```
sudo singularity build fpt_singularity_container.sif fpt_singularity_container.def &> fpt_singularity_build.log
```
This needs singularity >= 3.10.5. To install singularity follow these instructions https://apptainer.org/user-docs/master/quick_start.html. Once the container is build you can access the image as follows. Locate the .sif file and save the absolute path in the variable path_to_sif. To activate the container type: 

```
singularity shell ${path_to_sif}
```

You can also bind specific folders in your local system by using the B flag. 

```
singularity shell -B /path/to/bind:path/to/bind ${path_to_sif}
```

Within the container the directory is mounted at /path/to/bind. 
