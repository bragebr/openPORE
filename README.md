# openPORE
openPORE is an open repository aiming to solve advection-diffusion-reaction flow problems in complex media. 

# Usage (LINUX-BASED SYSTEMS) 
To run the examples provided in this repository, the openLB framework is needed. Download through:
~~~
git clone https://gitlab.com/openlb/release.git
~~~
and preferably install/update

~~~
brew install gcc tinyxml openmpi 
~~~

In the openLB directory (olb-1.7r0, depending on release version), open config.mk and change following flags

~~~
CXX := mpic++
CC  := gcc-(your version)
PARALLEL_MODE := MPI
USE_EMBEDDED_DEPENDENCIES := ON
~~~

for efficient use of MPI when running examples.

Then, in the openLB directory, do

~~~
mkdir myExamples
cd myExamples
mkdir reactingCylinder
~~~

In this folder, feel free to download the example provided (reactingCylinder.cpp) in this repository. To run this program, include the Makefile in the same folder and run 
~~~
make
mpirun -np X reactingCylinder
~~~
where X is your chosen number of processors. The results of the simulation is stored in the /tmp folder. For visualization, I recommend using ParaView to open the reactingCylinder.pvd file in the /tmp folder. 





