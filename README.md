# First Passage Time on Linear Framwork Graphs

A C++ package for First Passage Time calculations on gene regulatory networks using the Linear Framework, partially based on kmnam/markov-digraphs from Kee-Myoung Nam.

To build the singularity container run: 
```
sudo singularity build fpt_singularity_container.sif fpt_singularity_container.def &> fpt_singularity_build.log
```
This needs singularity >= 3.10.5. To install singularity follow these instructions https://apptainer.org/user-docs/master/quick_start.html

To build fpt has outside the container follow these steps (assuming ubuntu). ${path_to_FPT_repo} represents the path to the FPT repository you cloned.  
1) Download libmpfr and libmpc headers:
   
    ```
    sudo apt-get install -y build-essential
    sudo apt-get install -y libmpfr-dev
    sudo apt-get install -y libmpc-dev
    ```
    
2) Download python3 and pybind11:

    ```
    apt-get install -y python3.8 python3-pip
    pip3 install pybind11
    ```
3) In ${path_to_FPT_repo}/lib download boost-1.81.0

    ```
    cd FPT/lib
    wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
    tar -xf boost_1_81_0.tar.gz
    rm boost_1_81_0.tar.gz
    ```

4) in ${path_to_FPT_repo}/lib download eigen-3.4.0

    ```
    cd ${path_to_FPT_repo}/lib
    wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    tar -xf eigen-3.4.0.tar.gz
    rm eigen-3.4.0.tar.gz
    ```

5) install fpt
    
    ```
    cd ${path_to_FPT_repo}
    make
    make install
    make clean
    export PYTHONPATH=$PYTHONPATH:/home/fpt/bin
    ```





