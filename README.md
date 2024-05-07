# RDFanalysis
Data analysis with RDataFrame by ROOT CERN

UI
is used for data translation, pre-cleaning, calibration and analysis. The task is chosen by commenting/uncommenting one of the below functions. Geometry is chosen by changing runNo flag in constants.h file - 1,2,3 for each experimental data geometry and 11, 12, 13 for each simulated geometry.
	FUNCTIONS
	translator - converts original .root data to Root DataFrame (RDF) readable format. It runs TSelector defined in simplifier123.cxx file. No data or accuracy is lost, it's just a matter of the original tree structure that RDF couldn't use 
	cleaner - removes events which couldn't be used for the cross-section calculation. Those include wrong ToF range, no particle coincidence in the left and the right telescopes and corrupted MWPC events. 
	calibratorCaller - calibrates detectors, calcuates ToF and average energy loss in ToF F5
	analysis - reconstructs reaction vertex, scattering angles in LAB and CM, calculates missing masses of the particles. After applying 1D and 2D cuts, the data can be saved to smaller output files.

simplifier - small script for simplifying tree structure for ROOT DF

modifiedAnalysis
is used for calculating cross-sections from event angular distributions
	FUNCTIONS
	realDataAnalyzer - calculates event angular distribution in the CM based on files created with UI. Histogram is created and saved in .root file
	calculateEfficiency - calculated the detection efficiency based on the GEANT4 simulated data for both hydrogen isotopes and for each geometry. Saves output histogram to .root file
	calculateXsection - based on the data obtained in the above functions differential cross-section is calculated for each geometry and ion. Saves output histogram to .root file
	makeGraph - merges the cross-sections along the geometries and creates one for each isotope

	functions with "DT" suffix perform the same for the (d,t) reaction data


constants - file with physical constants and geometry info

TODO
-load .root data files lists to lists
-perform "translator", "cleaner", "calibratorCaller" and "analysis" on lists of files
-remove unused methods
-create function for calculating kinematics
