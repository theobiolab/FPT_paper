# First Passage Time on Linear Framwork Graphs

A C++ package with Python binding (Pybind11) for First Passage Time calculations on gene regulatory networks using the Linear Framework, partially based on [kmnam/markov-digraphs](https://github.com/kmnam/markov-digraphs.git) from Kee-Myoung Nam.

To build the singularity container run (soon will be uploaded to singularity-hub): 
```
sudo singularity build fpt_singularity_container.sif fpt_singularity_container.def &> fpt_singularity_build.log
```
This needs singularity >= 3.10.5. To install singularity follow these instructions https://apptainer.org/user-docs/master/quick_start.html

# Build outside the container

Rules for building outside the container are present in the makefile. These are based on ubuntu>=20.04 (for now). If you have a windows system we recomend installing Windows Sub Linux with ubuntu 22.04 (https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview). 

1) clone the repository

2) execute the following command to install the required dependencies (requires sudo permission): 

   ```
   make ubuntu_setup
   ```

3) To compile and install FPT binaries execute (warnings are ok): 

   ```
   make all_ubuntu
   make install 
   make clear
   ```

   the binaries are now installed in the ```bin``` folder
   
4) To access the modules from python we need to add ```bin``` to python path

   ```
   export PYTHONPATH=$PYTHONPATH:$(path_to_repo)/bin
   ```




# How to enter the container
Locate the .sif file and save the absolute path in the variable path_to_sif. To activate the container type: 

```
singularity shell ${path_to_sif}
```

You can also bind specific folders in your local system by using the B flag. 

```
singularity shell -B /path/to/bind:path/to/bind ${path_to_sif}
```

Within the container the directory is mounted at /path/to/bind. 
