This 02_data_analysis README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: 02_data_analysis

2. Author Information
	A. Principal Investigator Contact Information
		Name: Rainer Machne
		Institution: Heinrich Heine University - Institute for synthetic microbiology
		Address: 
		Email: machne@hhu.de

	B. Associate or Co-investigator Contact Information
		Name: Tobias Pfennig
		Institution: Heinrich Heine University - Institute for synthetic microbiology
		Address: 
		Email: tobias.pfennig@hhu.de

	C. Alternate Contact Information
		Name: 
		Institution: 
		Address: 
		Email: 

3. Date of data collection (single date, range, approximate date): 2020-12-31

4. Geographic location of data collection:

5. Information about funding sources that supported the collection of the data: 


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: 

2. Links to publications that cite or use the data: 

3. Links to other publicly accessible locations of the data: 

4. Links/relationships to ancillary data sets: 

5. Was data derived from another source? no
	A. If yes, list source(s): 

6. Recommended citation for this dataset: 


DATA & FILE OVERVIEW

1. File List: 
	RUNFIRST.R: R script meant to be executed when dataset is obtained first

	XXXXXXXX_analysis_SDYXXX_X.R: R scripts analyzing a specific study
	
	XXXXXXXX_analysis_SDYXXX_X_X.R: R scripts analyzing a specific sub-study
	
	XXXXXXXX_custom_functions_vX.R: R script containing custom function defined for the anaylsis
	
	data_analysis_outputs: folder containing outputs of data_analysis scripts, notably the linear models used for drymass estimation

2. Relationship between files, if important: 

3. Additional related data collected that was not included in the current data package: 

4. Are there multiple versions of the dataset? no
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? 
		ii. When was the file updated? 


METHODOLOGICAL INFORMATION

Run the RUNFIRST.R script before running any other script. For this first adjust the line
	
	INVESTIGATIONPATH <- "" #change to local investigation path
	
to point towards the location of the investiation folder.
This script also installs all needed packages for the analysis.