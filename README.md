# First Passage Time on Linear Framework Graphs

A C++ package with Python binding (Pybind11) for First Passage Time calculations on gene regulatory networks using the Linear Framework, partially based on [kmnam/markov-digraphs](https://github.com/kmnam/markov-digraphs.git) from Kee-Myoung Nam. The code has been tested on macOS Sonoma version 14.0, Ubuntu 20.4 and Centos7

> :warning: **Compiler version:** to compile the code you need gcc/9.2.0 on your path!!!!


# Compile the models

Regardless of the system you are using, the first step is to create a conda environment with the minimal packages requirments for compiling. To install minconda on your system follow these instructions: https://docs.conda.io/projects/miniconda/en/latest/. Once conda is installed, you can proceed with creating the environment:
```
cd ${path_to_fpt_repo}
conda create -n "fpt_env" python=3.11.5
conda activate fpt_env
conda install pip
pip install -r minimal_requirments.txt
```
> :warning: **Pybind11 includes:** before proceeding any further double check the presence of pybind11 includes in the environment by typing ``` python -m pybind11 --includes ```. The output should be an absolute path to a ```python3.11``` include file and the ```pybind11``` package include folder within the conda ```fpt_env``` folder in your home (double ckeck).

If the pybind11 includes are in the conda env you can procede with compiling:
```
make dependencies
```
this installes locally all the c/c++ dependencies required (boost, eigen, gmp, mpfr, mpc). Please do this even if you have previously installed these packages, the code is version specific. To test that everything went fine during this step you can compile a series of tests: 
```
make tests
```
if this step is succesfull you will see binaries in the test folder. To compile the models run
```
make ladders 
```
The binaries are stored in the ```/bin``` folder. 

> :warning: **Importing the modules:** before importing the modules in python you need to source the ```config.sh``` file! This sets the dynamic library paths to the /FPT/lib folder. 

Enjoy :D

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
