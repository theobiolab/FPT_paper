# First Passage Time on Linear Framework Graphs

A C++ package with Python binding (Pybind11) for First Passage Time calculations on gene regulatory networks using the Linear Framework, partially based on [kmnam/markov-digraphs](https://github.com/kmnam/markov-digraphs.git) from Kee-Myoung Nam. The code has been tested on macOS Sonoma version 14.0, Ubuntu 20.4 and Centos7

> :warning: **Compiler version:** to compile the code you need gcc/9.2.0 or higher on your path!!!!

# Paper information
This repository contains all the data and code to reproduce the figures and results in the paper by Ravanelli, Nam, Gunawardena, Martinez-Corral 2025. Reference will be updated upon publication. 

The notebooks used to generate the figures are located in the Data folder after unzipping. The steps for manipulating the analytical equations presented in the text are outlined in the data/analytics notebooks. 

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


# How to run the Docker container locally

The docker-compose file allows to build the image and run the container with one command *from the root directory*:

```
docker compose up
```
The image will build only once. After the first time, the image will be cached and Docker won't build it again,
unless you force it by adding the ```--build``` flag. If you make some changes after the first build, make sure to rebuild the image.

```docker compose up``` will start the Jupyter Lab server and connect your terminal to the shell inside the container. 
When the container is started, you should see some urls printed in the terminal. To access Jupyter Lab from your browser,
use the url starting with ```http://127.0.0.1:8888/lab```. If you prefer, you can also substitute the ip ```127.0.0.1``` with the ```localhost``` domain.

To stop the container, you can hit ```CTRL``` + ```C``` in the terminal. This will only stop the container, but won't remove it.
To stop and remove the container, run the following command *from the root directory*:

```
docker compose down
```


# How to prevent files from being copied in the Docker container

If you want to exclude files from the container, add them in the ```.dockerignore``` file.


# How to share the image on Docker Hub

Once the image is ready to be shared, you can create an account on Docker Hub and create a repository.
The repository will be identified with the following name:
```
<account>/<repository>
```
This will also be the image name.

To build the image and name it with the correct repository name, run from the root directory:
```
docker build -t <account>/<repository> .
```
If you want to version your image, you could also specify a tag using the format ```<account>/<repository>:<tag>```.

Once the image is built, you can push it to Docker Hub with
```
docker push <image-name>
```
And you can download it with
```
docker pull <image-name>
```
