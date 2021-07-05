# data 

This folder contains the basic data structures.

Each of the folders is related to a structure,
and there is a single function that initializes 
the data structures with empty fields.

This is useful to have a common representation 
of all the structures, and potential fields, 
so that it can be used as a "living" doc.

###
<details><summary><b>Structure</b></summary>

The different folders are:

* **+anatomicalModel**: 
Anatomical model with discretization and properties.
To be fed to the spinModel generation.

* **+coilModel**:
Structure with all the available coils, the active one, 
and for each of the coils available, details, maps and 
related info.
To be fed to the spinModel generation.
To be used in SENSE reconstruction.

* **+spinModel**:
Structure with the model to simulate.
It contains the geometric information of the slices
given by the Front End, and the interpolated data
that is going to be fed to the simulator. 

* **+mrSystem**:
Structure with basic information of the 
MR system, including field stregnth, 
gradient strength, slew rate, etc... 

* **+pulseSequence**:
Structure with the pulse sequence data, 
as well as the sequence signals.
To be fed to the simulator, and used in 
the reconstruction. 

* **+reconData**:
Structure with the results of the simulation
(signals), information for the reconstruction,
and that will store the (re-ordered) K-space
and the generated Image space.

* **+expControl**:
Structure with controls and configuration
for the current experiment.
It will store all flags and values for the 
simulation, as well as functional mode, debugging, etc.

</details>

###
<details><summary><b>TODO</b></summary>

* Complete and incorporate more structures and more data. 

</details>

###
<details><summary><b>Installation and Dependencies</b></summary>

Self-contained package in Matlab.
No dependencies.

</details>

###
<details><summary><b>Testing</b></summary>

Not functional module, so no testing needed.

</details>


###
<details><summary><b>References</b></summary>

* **[1]** 
* **[2]** 

</details>

