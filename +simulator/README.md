# simulator 

Simulation module:

* This module contains the main simulation engine, 
and interfaces with the CUDA kernels for the Bloch simulation.

###
<details><summary><b>API</b></summary>

The main functions interfaces to the kernel are:

* To run the standard analytical kernel:

**>> [solution] = simulator.bloch.kernel.runAnalytical(solution,pulseSequence.subSeq{jj},spinModel.slice{ii}.model,expControl);**

* To run the Phasor version, with accurate signal integration:

**>> [solution] = simulator.bloch.kernel.runPhasor(solution,pulseSequence.subSeq{jj},spinModel.slice{ii}.model,expControl);**

* To run the Self-Weighted Diffusion (SWD) version, with intra-voxel diffusion and signal integration:

**>> [solution] = simulator.bloch.kernel.runSWD(solution,pulseSequence.subSeq{jj},spinModel.slice{ii}.model,expControl);**


</details>


###
<details><summary><b>TODO</b></summary>

* Refine interfaces when more functionalities are added. 

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

* **+bloch**:
Module to simulate basic Bloch equations.
Includes Analytical and Phasor versions.

* **+diffusion**:
Module to simulate Bloch-Torrey, including 
self-diffusion effects. 

* **+motion**:
Module to simulate Bloch with motion in time.
 
* **+test**:
Has the testing functions. See below for details. 

* **+example**:
Has some scripts to visualize some features and 
results. Can be useful for debugging or verifying 
correctness or results. See below for details. 

</details>


###
<details><summary><b>Testing</b></summary>

The test are based on Matlab unit testing,
with parameterized runs.

They will run a set of different tests
for different configurations.

To run all the tests, call:

**>> [results] = simulator.test.runTestSuite();**

The console will show progress and will report fails.
To see more details, you can use the **results** structure
by doing:

**>> table(results);**

To run specific test suite, call:

**>> [results] = simulator.test.runTestSuite(<type>);**

where **type** is a string with the type of test.
Currently supported types are:
**'flipAngle'** / **'signalIntegration'** / **'diffusion'** / **'all'**.

For more details on the specifics of the tests, see the corresponding 
function documentation.


</details>

###
<details><summary><b>Examples</b></summary>

There are a set of scripts that work as examples, to illustrate 
some features, basically ploting some comparisons.
These can be useful to qualitatively asses the results or 
verify the correctness of the physics.

* Flip angle illustration, where the RF and slice selection can be 
modified to see the effects of the RF and final magnetization 
along the slice profile:

**>> simulator.example.flipAngleExample;**

* Signal integration, and how the Phasor performs the correct voxel
spatial integration, and the analytical needs higher resolution to 
converge:

**>> simulator.example.signalIntegrationExample;**

* Effect of a 180 refocusiong RF pulse, and how it affects the 
phase and its derivatives w.r.t. the spatial directions.
This is rllevant for the signal integration and the diffusion:

**>> simulator.example.rf180Example;**

* Effect of the diffusion, and comparison with the theoretical
formula, for a simplified Pulsed Gradient (PGR) reversed sequence,
with no RF, and for a Pulsed Gradient Spin Echo (PGSE), 
where a 180 RF is used for the reverse:

**>> simulator.example.diffusionPGRExample;**

**>> simulator.example.diffusionPGSEExample;**



</details>



###
<details><summary><b>References</b></summary>

* **[1]** 
* **[2]** 

</details>

