# README #
Sharp GUI ver0.4 Demo without sample electrophysiology data.
Mini GUI project for collaborating lab in 2016.

For finding signal peaks and bursts in electrophysiology recordings,
like EMG / sharp microelectrode / patch clamp

### Usage ###
	need 'plotpkmarkers.m', 'sharp_gui.m', and 'sharp_gui.fig' in same folder
	use MATLAB import data tool, to load 'SharpEmgRecording.xlsx' spreadsheet.
	After loading data into the base workspace, in the box next to Data Name, 
	type in the data set name, such as 'data_P18', then click 'Import Data', 
	which loads into current program workspace. 
	After setting parameters, click 'Find Peaks', which plots data, peaks and other statistics.

	Or use 'import_xls' to load two-column (time; voltage) spreadsheet.

### What is this repository for? ###
GUI to find signal peaks realiablely in electrophysiology recordings,
like EMG / sharp microelectrode / patch clamp
and display statistical/spectrum analysis.

### How do I get set up? ###
Clone repo.
Load recordings.
Use default parameters for myometrium cell sharp recordings.
Tune parameters for other recordings.

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###
* Repo owner: Harry
