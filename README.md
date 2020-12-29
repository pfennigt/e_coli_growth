# e_coli_growth
An investigation of the growth of Escherichia coli strains W3110Z1 and W3110 using optical measurement of culture parameters and measurements of drymass and metabolites.


# Structure information

This project is structured in the style of ISA: https://isa-tools.org/index.html
	investigation: "High level concept to link related studies"
		The overall project of investigating the E. coli growth
		This uppermost folder level will also be called the investigation folder
		The data analysis scrips expect this folders path to be set as the INVESTIGATIONPATH variable

	study: "The central unit, containing information on the subject under study, its characteristic and any treatment applied. A study has assays"
		The different experiments performed during the investigation
		If multiple measurements are performed on the same samples or during the same protocol, they belong to the same study
		A study may be logically split into sub-studies, e.g. if different carbon sources or strain were used in parallel in the protocol. In this case the sub-studies are separately read and analyzed
		Can be found in folder 01_studies. The file studies_index.csv gives information to each study and possible sub-studies
		Studies are consecutively named with an identifier (e.g. "SDY001"), a short token for better understanding (e.g. "HighGA"), and a short string explining the experiment further (e.g. "gradient_GA_test")
		
		The folders containall relevant assays, supplementary data/ files, and an explanatory README.md file
		
	assay: "Test performed either on material taken from the subject or on the whole initial subject, which produce qualitative or quantitative measurements (data)"
		The single measurements performed within a study, leading directly to data
		Are consecutively named with an identifier (e.g. "ASY001") and a description of the performed measurement/ the collected data (e.g. "biolector_data")
		
		The folders contain all relevant data, metadata, supplementary data, a read_data.R file (or multiple for sub-studies), and an explanatory README.md file


Folder structure in the investigation folder:
	01_studies: Contains studies and their assays
	
	02_data_analysis: Data analysis performed on the study data using R
	
	03_figures: Relevant figures, especially figures resulting from datanalysis are saved here
	
	04_presentations: relevant presentations regarding the investigation
	
	05_posters: relevant posters regarding the investigation
	
	06_manuscripts: relevant documents regarding the investigation


# Dataset information

This e_coli_growth README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: e_coli_growth

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
		Name: Jonas Burmester
		Institution: Heinrich Heine University - Institute for synthetic microbiology
		Address: 
		Email: jonas.burmester@hhu.de

3. Date of data collection (single date, range, approximate date): February 2020 to October 2020

4. Geographic location of data collection: 51.19055739134027, 6.7904520065037675

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
	01_studies: Contains studies and their assays
	
	02_data_analysis: Data analysis performed on the study data using R
	
	03_figures: Relevant figures, especially figures resulting from datanalysis are saved here
	
	04_presentations: relevant presentations regarding the investigation
	
	05_posters: relevant posters regarding the investigation
	
	06_manuscripts: relevant documents regarding the investigation

2. Relationship between files, if important: 

3. Additional related data collected that was not included in the current data package: 

4. Are there multiple versions of the dataset? no
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? 
		ii. When was the file updated? 


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
	Information to used media and references can be found in 06_manuscripts/Studying_Multiple_Growth_Phases_of_E_coli_on_Minimal_Medium_with_Experiment_and_Theory_Tobias_Pfennig.pdf

	Experiments were performed on Escherichia coli strain W3110Z1 (genotype: laciq, PN25-tetR, SpR, IN(rrnD-rrnE)1, rph-1, (see: ATCC 39936); received from: Expressys (Dr. Rolf Lutz)) or strain W3110 (genotype: IN(rrnD-rrnE)1, rph-1, (see: ATCC 39936)).
	
	Typically the Bacteria were defrosted and held in preculture using M9 minimal medium with either glucose or acetate as carbon source. Then, the bacteria were pipetted into media with predefined conditions and grown in a BioLector Pro (m2p-labs GmbH, Baesweiler, Germany) bioreactor. There, the bioreactor recorded optically measured culture parameters: backscatter, riboflavin fluorescence, pH , O2 saturation, and NADH fluorescence.
	
	Depending on the experiment time series sampling was performed. Then, the dry mass or metabolite concentrations (acetate, glucose, or glycogen) were measured.
	
	
2. Methods for processing the data: depending on the study

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)
	
	BioLection (V.3.16.74.0, optional)
	
	Opentrons App for use with OT-2 robot (3.18.1)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: depending on the study

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig, Jonas Burmester (depending on the study)