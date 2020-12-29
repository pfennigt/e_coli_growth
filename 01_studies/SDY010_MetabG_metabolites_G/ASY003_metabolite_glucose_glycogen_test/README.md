This ASY003_metabolite_glucose_glycogen_test README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: SDY010_MetabG_metabolites_G/ASY003_metabolite_glucose_glycogen_test

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

3. Date of data collection (single date, range, approximate date): 2020-10-16

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
	20201016_metabolite_glucose_glycogen_test_depr.csv: DEPRICATED, CLARIOstar glucose and glycogen measurements of test samples with highly set gain
	
	20201016_metabolite_glucose_glycogen_test_METADATA.csv: metadata regarding the analyzed samples
	
	20201016_metabolite_glucose_glycogen_test_best_gain.csv: RECOMMENDED, CLARIOstar glucose and glycogen measurements of test samples with better gain settings
	
	read_data.R: reads metabolite measurements and metadata data into R

2. Relationship between files, if important: 
	20201016_metabolite_glucose_glycogen_test_best_gain.csv measures the same 96 well plate as 20201016_metabolite_glucose_glycogen_test_depr.csv but a better gain setting was used
	
	20201016_metabolite_glucose_glycogen_test_METADATA.csv contains the metadata of the samples in 20201016_metabolite_glucose_glycogen_test_best_gain.csv (and 20201016_metabolite_glucose_glycogen_test_depr.csv)
	
	read_data reads 20201016_metabolite_glucose_glycogen_test_best_gain.csv.csv and 20201016_metabolite_glucose_glycogen_test_METADATA.csv.csv into R

3. Additional related data collected that was not included in the current data package: 

4. Are there multiple versions of the dataset? no
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? 
		ii. When was the file updated? 


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
	Information to used media and references can be found in ../../../06_manuscripts/Studying_Multiple_Growth_Phases_of_E_coli_on_Minimal_Medium_with_Experiment_and_Theory_Tobias_Pfennig.pdf

	Escherichia coli strain W3110Z1 was defrosted from -80 °C freezer sample #DD13 (genotype: laciq, PN25-tetR, SpR, IN(rrnD-rrnE)1, rph-1, (see: ATCC 39936); received from: Expressys (Dr. Rolf Lutz)).
	
	The bacteria were plated on an agarose Petri dish containing 35 g l^-1 LB broth and 100 mg ml^-1 spectinomycin. The inoculated plate was left to incubate at 37 °C overnight. The bacteria were transferred onto new agarose plates and incubated at 37 °C multiple times to acchieve single, monoclonal colonies. The last plate was placed in refrigeration and used as the source for inoculations of liquid media.
	
	Liquid cultures were created by picking a colony from the source-plate and transferring it into 5 ml of minimal M9 medium with growth-supporting supplements, antibiotics, and 22 mM (0.4 % w/v) glucose ("M9G", Tab. S2). This liquid culture was then incubated at 250 rpm and 37 °C. Every two to seven days, 50 ul of this liquid culture were transferred into 5 ml of fresh M9G medium and continued to be incubated. One day before the measurement the culture was renewed as explained before.
	
	This experiment served as a test of protocols and accuracy, as well as the conversion of glycogen to glucose concentrations. The Sigma aldrich Glycogen Assay Kit MAK016 was used for glucose and glycogen quantification. The test samples analyzed were the supernatant of an overnight preculture (uen), the supernatant (ex) and cell pellet (exg) of 15 ml (???) of a preculture approximately in exponential phase, a sample of purified glycogen (glyc) and ???. Also, each four glucose (S1_glc to S4_glc) and glycogen standards (S1_glyc to S4_glyc) were measured. Cell pellets for glycogen measurements were first washed with MilliQ and then cooked in 400 µL KOH (30% w/v) for 2 h. Then glycogen was  precipitated with ice cold 74% ethanol and incubation at -20 °C for at least 2 h. The samples were washed twice with 1 ml of 70 % and 90 % ethanol and vacuum dried at 60 °C for 20 min. Afterwards the samples were resuspended in 50 ul kit buffer and treated and measured according to the kit instructions. Glucose measurements were performed by skipping the first hydrolysis step of the kit. The fluorescences (lambda_ex = 535-15 nm/lambda_em = 587-20 nm) were measured in a black, clear bottom 96 well plate.
	
	The same well plate was measured twice with different gain settings. First, a gain setting allowing for the measurement of all samples and standards was used (recommended, 20201016_metabolite_glucose_glycogen_test_best_gain.csv). Additionally, a gain setting adjusted to the highest sample concentration was tested (20201016_metabolite_glucose_glycogen_test_depr.csv), leading to some standard values to exceed the detection maxiumum (260000 AU in Raw Data).
	
2. Methods for processing the data: 

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: visual comparions of data curves with replicates and previous experiments, data is unchanged

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig


DATA-SPECIFIC INFORMATION FOR: 20201016_metabolite_glucose_glycogen_test_depr.csv

1. Number of variables: 4

2. Number of cases/rows: 22

3. Variable List: 
	Well: well of the 96 well plate
	
	Content: identifier for the well content  given by the CLARIOstar, corresponds to well_id in 20201016_metabolite_glucose_glycogen_test_METADATA.csv
	
	Raw Data (535-15/587-20): raw fluorescence data (em: 587 nm, bandpass: 15 nm; ex: 535 nm, bandpass: 20 nm) [AU] (recommended for analysis)
	
	Linear regression fit based on Raw Data (535-15/587-20): estimated glucose content in the digested well from a linear model of the standard wells [ug]

4. Missing data codes:

5. Specialized formats or other abbreviations used:


DATA-SPECIFIC INFORMATION FOR: 20201016_metabolite_glucose_glycogen_test_METADATA.csv

1. Number of variables: 8

2. Number of cases/rows: 13

3. Variable List: 
	well_id: identifier for the well content given by the CLARIOstar, corresponds to Content in 20201016_metabolite_glucose_glycogen_test_best_gain.csv and 20201016_metabolite_glucose_glycogen_test_depr.csv
	
	sample_id: identifier for the well content given during the experiment
	
	time_h: runtime of the biolector at sample collection [h]
	
	metabolite_type: measured metabolite
	
	sample_type: do the wells contain samples or standards
	
	standard_concentration_mug: set concentration in standard wells [ug]
	
	sample_dilution: dilution of the sample before measurement (diluted to x times concentration)
	
	well_volume_mul: liquid volume in the measured well [ul]

4. Missing data codes:

5. Specialized formats or other abbreviations used:


DATA-SPECIFIC INFORMATION FOR: 20201016_metabolite_glucose_glycogen_test_best_gain.csv

1. Number of variables: 4

2. Number of cases/rows: 22

3. Variable List: 
	Well: well of the 96 well plate
	
	Content: identifier for the well given by the CLARIOstar, corresponds to well_id in 20201016_metabolite_glucose_glycogen_test_METADATA.csv
	
	Raw Data (535-15/587-20): raw fluorescence data (em: 587 nm, bandpass: 15 nm; ex: 535 nm, bandpass: 20 nm) [AU] (recommended for analysis)

	Linear regression fit based on Raw Data (535-15/587-20): estimated glucose content in the digested well from a linear model of the standard wells [ug]

4. Missing data codes:

5. Specialized formats or other abbreviations used: