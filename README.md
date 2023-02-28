# First Passage Time on Linear Framwork Graphs

A C++ package for First Passage Time calculations on gene regulatory networks using the Linear Framework partially based on kmnam/markov-digraphs from Kee-Myoung Nam.

To build the singularity container run: 
```
sudo singularity build fpt_singularity_container.sif fpt_singularity_container.def &> fpt_singularity_build.log
```
This needs singularity >= 3.10.5. To install singularity follow these instructions https://apptainer.org/user-docs/master/quick_start.html

