CRISM MLM Production
Version 4.7.0
25 April 2017

*** For a list of recent changes, please see CHANGELOG_MLM.txt ***

-------------------------------------------------------------
The Program's expections about input
-------------------------------------------------------------

* SSA & DDR must have the same number of rows
* SSA & DDR must both have the full 640 columns (uncropped)
* SB & WA must be for the entire array on the instrument (uncropped)
* The selected SSA, DDR, SB, and WA files must have corresponding ENVI header
  files, whose names differ only in the extension being changed to .hdr

-------------------------------------------------------------
Running the code through the GUI 
-------------------------------------------------------------

The GUI can be opened by executing crism_gui in the Matlab main window, or in
many cases by double-clicking the file crism_gui.fig in this directory.

Options in the GUI parallel those in test5v5_compact.m as listed below. Some
of the options have tooltips that describe their function more clearly. This
help information is a work in progress.

Since for a given set of bands (L or S) the same SB and WA files are used for
every run, a button has been added to find those files quickly. If the
appropriate SB and WA files are found in the same directory as the user has
browsed to for the SSA file, the paths to those files are automatically
entered in the GUI.

An additional automatic feature is a button that sets the directory for output
to be the same as the input directory (i.e., the directory that the SSA the
user selected is in). This is the typical case.

-------------------------------------------------------------
Running the code by modifying the test5v5_compact.m script
-------------------------------------------------------------

All parameters to the algorithm are found in test5v5_compact.m in lines 23-94.
These are the only lines that need to be changed during routine use.

These parameters include, from top to bottom:
    - Penalty function settings
    - Switch for L vs S data
    - Setup of the scene (pixel size, etc) and its general geometry settings
    - Number of iterations to run
	- Whether to apply switchover (and if so, in what iteration and with what code)
    - Spectral offset in nm (to apply to bandpass centers given in wa file)
    - Data Saving controls (can choose separately whether to save the scene and/or sensor-space estimates during specified iterations)
    - Input file names
        - Datafiles: crism_iof (sensor-space data), ddr (derived data record), sb (spectral bandpasses), wa (bandpass centers)
        - Off-nadir-only datafiles: radius, topo, spice
    - Input file subsetting
    - Output file names (observation_name, which is used as a prefix in filenames when saving the scene or estimate of sensor-space measurements)
    
After setting these parameters, the code can be run by clicking the green
arrow run button (in the Matlab editor), or by executing test5v5_compact in
the Matlab main window.

-------------------------------------------------------------
Running the code using test5v5_compact as a function
-------------------------------------------------------------

The file test5v5_compact.m is a Matlab function that requires no arguments.
This allows the above use of it as if it were a script. It is also possible
to call it as a function, and give arguments that replace the some of the 
parameter values hard-coded into the test5v5_compact.m file. This is done
by giving test5v5_compact pairs of arguments, where the first is a string
specifying the name of an algorithm parameter and the second is the value
which that parameter should be set to in the run. For example, the call:

> test5v5_compact('beta_spec',0.3,'observation_base','ATO0002EC79_test');

runs the algorithm as if the script's text were edited to make beta_spec = 0.3
and observation_base = 'ATO0002EC79_test'.

This technique is useful for running the algorithm many times with sets of
parameters that are mostly the same (but with a few that are systematically
changed). It is also useful for submitting jobs to computing clusters
supporting MATLAB jobs, such as the Center for High Performance Computing
(CHPC) at the Washington University Medical School.

-------------------------------------------------------------
About the C++ MEX Components
-------------------------------------------------------------
Since version 4.6.0, some time-intensive components of the algorithm have had
C++ translations of the Matlab originals. This folder includes the C++ source,
in .cpp and .h files, as well as those files compiled for Windows in .mexw64
files. These files are kept up-to-date with each release. If these files were
missing, the Matlab equivalents would be used as a fallback.

The C++ versions make use of OpenMP, a system for specifying ways to parallelize
C++ code to spread computations across several processors. This ability in the
primary reason that translations were written in the first place (Matlab being
a more difficult environment to take advantage of this potential).

You might want to build MEX binaries from these source files if:
	A) you've made changes to the C++ code and want the modified version to be
	   run instead of the prepackaged version, or
	B) you have a non-Windows (or non-64-bit Windows) operating system and want
	   to be able to run the C++ version of these components rather than the
	   Matlab equivalent.
	
The steps to do this are roughly:
	1) Install a C++ compiler supported by your version of Matlab for compiling
	   MEX files (see http://www.mathworks.com/support/compilers/R2016a/index.html
	   or your Matlab version's equivalent).
	2) Use the command 'mex -setup' at the Matlab prompt to configure it to use
	   your compiler.
	3) Use the command 'make omp' to use the included 'make' utility to build
	   the relavent C++ files with OpenMP multiprocessing hooks. I can't
	   guarantee that the compiler flags I used in the MexMakeFile are
	   compatible with your compiler (these work under Microsoft Visual Studio
	   2015), and it may work better to call the 'mex' command directly with
	   your preferred compiler options.

-------------------------------------------------------------
A note on versioning
-------------------------------------------------------------
Version numbers since 4.0.0 done in accordance with Semantic Versioning 2.0.0
(semver.org). The format is MAJOR.MINOR.PATCH. A release that has non-backward
compatible changes will be the next MAJOR version. A release that adds features will
be the next MINOR version. A release that doesn't do either of these but only e.g. fixes
bugs increments PATCH.

-------------------------------------------------------------
Contact
-------------------------------------------------------------

Since so many things are still changing in this system, it is not unexpected
that there will be bugs, poorly-explained features, etc. For questions and
feedback at this early stage, it is probably most efficient for you to contact
me, Daniel Politte, at dvpolitte <at> wustl <dot> edu.
