# SpinTwin 

Simulation module:

* This module contains the main simulation engine, 
and interfaces with the CUDA kernels for the Bloch simulation.

###
<details><summary><b>API</b></summary>

The main interface is:

**>> [timeSolution, stats] = spinTwin.runSimulation( ...
    simModel, pulseSequence, motionModel, simControl, dbgControl);**

</details>

###
<details><summary><b>Installation and Dependencies</b></summary>

Self-contained package in Matlab.
No dependencies.

Kernel run requires correct CUDA kernel .ptx file available.

</details>

###
<details><summary><b>Structure</b></summary>

The different folders are:

* **+fwdBloch**:
Module to perform the forward Bloch Simulation: from data generate signals.

* **+setup**:
Basic set up of standard variables and structures. 
 
* **+test**:
Has the testing functions. See below for details. 
Inside the **+test** folder, we have:

* **+test/+unit**:
Unit tests.

* **+test/+example**:
Scripts that allow run cases similar to the tests, 
allowing to vizualize the results and compare.
Check the individual scripts for more details.

* **+test/+performance**:
Scripts that allow run performance runs.
These run all the cases for a FID and for an RF.
It generates an output with errors and performance times, 
as well as plots with the performance results.

* **+test/+model**:
Auxiliary functions to generate the models.

* **+test/+seq**:
Auxiliary functions to generate the sequences.

</details>


###
<details><summary><b>Testing</b></summary>

The test are based on Matlab unit testing,
with parameterized runs.

They will run a set of different tests
for different configurations, which are under
**+test/**:

To run all the tests, call:

**>> [results] = spinTwin.runUnitTest();**

The console will show progress and will report fails.
To see more details, you can use the **results** structure
by doing:

**>> table(results);**

To run specific test suite, call:

**>> [results] = simulator.test.runTestSuite(<type>);**

where **type** is a string with the type of test.
Currently supported types are:
**'fid'** / **'flipAngle'** / **'sliceSelection'** / **'signalIntegration'** / **'diffusion'** / **'all'**.

For more details on the specifics of the tests, see the corresponding 
function documentation.


</details>

###
<details><summary><b>TODO</b></summary>

* Refine interfaces when more functionalities are added. 

</details>

###
<details><summary><b>References</b></summary>

* **[1]** 
* **[2]** 

</details>

