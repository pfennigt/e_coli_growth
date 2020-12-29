This ASY001_biolector_data README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: SDY001_HighGA_gradient_GA_test/ASY001_biolector_data

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

3. Date of data collection (single date, range, approximate date): 2020-02-18

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
	25-Ecoli_2020--2020-02-18-15-37-46.csv: raw data file exported from BioLector Pro

	25_Ecoli_2020_REDUCTION-1.csv: processed data from BioLection software

	25_Ecoli_2020_REDUCTION-1_layout.csv: wellplate layout description with applied factors
	
	plot_annotations_acetate.csv: manual annotations for points of interest in the BioLector Pro data for acetate samples
	
	plot_annotations_glucose.csv: manual annotations for points of interest in the BioLector Pro data for glucose samples
	
	read_data_acetate.R: reads acetate samples from BioLector Pro data and respective plot annotations into R
	
	read_data_glucose.R: reads glucose samples from BioLector Pro data and respective plot annotations into R

2. Relationship between files, if important: 
	25_Ecoli_2020_REDUCTION-1.csv is better readable, processed version of 25-Ecoli_2020--2020-02-18-15-37-46.csv

	25_Ecoli_2020_REDUCTION-1_layout.csv describes applied factors for each well in 25-Ecoli_2020--2020-02-18-15-37-46.csv and 25_Ecoli_2020_REDUCTION-1.csv
	
	plot_annotations_acetate.csv describes manually annotated points in acetate samples of 25-Ecoli_2020--2020-02-18-15-37-46.csv and 25_Ecoli_2020_REDUCTION-1.csv
	
	plot_annotations_glucose.csv describes manually annotated points in glucose samples of 25-Ecoli_2020--2020-02-18-15-37-46.csv and 25_Ecoli_2020_REDUCTION-1.csv
	
	read_data_acetate.R: reads acetate samples from 25_Ecoli_2020_REDUCTION-1.csv and plot_annotations_acetate.csv into R
	
	read_data_glucose.R: reads glucose samples from 25_Ecoli_2020_REDUCTION-1.csv and plot_annotations_glucose.csv into R

3. Additional related data collected that was not included in the current data package: 

4. Are there multiple versions of the dataset? no
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? 
		ii. When was the file updated? 


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
	Information to used media and references can be found in ../../../06_manuscripts/Studying_Multiple_Growth_Phases_of_E_coli_on_Minimal_Medium_with_Experiment_and_Theory_Tobias_Pfennig.pdf

	Escherichia coli strain W3110Z1 was defrosted from -80 °C freezer sample #DD13 (genotype: laciq, PN25-tetR, SpR, IN(rrnD-rrnE)1, rph-1, (see: ATCC 39936); received from: Expressys (Dr. Rolf Lutz)). 		
	
	The bacteria were first transferred into liquid LB medium (20 g l^-1 LB Broth (Lennox); Carl Roth GmbH + Co. KG, Karlsruhe, Germany) with 100 mg ml^1 spectinomycin and incubated at 250 rpm and 37 °C. After one change of medium, they were plated on an agarose Petri dish containing 35 g l^-1 LB broth and 100 mg ml^-1 spectinomycin. The inoculated plate was left to incubate at 37 °C overnight. The bacteria were transferred onto new agarose plates and incubated at 37 °C multiple times to acchieve single, monoclonal colonies. The last plate was placed in refrigeration and used as the source for inoculations of liquid media.
	
	Liquid cultures were created by picking a colony from the source-plate and transferring it into 5 ml of minimal M9 medium with growth-supporting supplements, antibiotics, and 22 mM (0.4 % w/v) glucose ("M9G", Tab. S2). Yeast extract was also added to the medium to a concentration of 1 g l^-1. This liquid culture was then incubated at 250 rpm and 37 °C. Every two to seven days, 50 ul of this liquid culture were transferred into 5 ml of fresh M9G medium and continued to be incubated. One day before the measurement the culture was renewed as explained before.
	
	5 ml of an 1 OD cell solution were prepared by first measuring the optical density at 600 nm wavelength of a 1:4 dilution of culture with MilliQ. An appropriate amount of cell culture was spun down and the cell pellet was resuspended in 5 ml of fresh M9 medium with growth-supporting additives and antibiotics ("M9S", Tab. S2) by pipetting. Each 50 ml of M9S, M9G, and 66 mM (0.4 % w/v) acetate ("M9A", Tab. S2) medium, as well as the 5 ml of 1 OD cell solution were placed in an OT-2 pipetting robot with P300 pipette module (Opentrons, New York, USA). A 48-well MTP-48-BOH well plate ("FlowerPlate"; m2p-labs GmbH, Baesweiler, Germany) was filled using the pipetting script ../ASY001_biolector_data_OTscript.py. The layout included seven-step gradients of both glucose and acetate with equal carbon-molarity ranging from 0 to 119.88 C-mmol l^-1. Each incoulated well was filled with 900 ul of appropriately mixed medium and 100 ul of cell solution. The 6 blank wells contained 1 ml of medium with different carbon sources and concentrations.
	
	The FlowerPlate was measured in a BioLector Pro (m2p-labs GmbH, Baesweiler, Germany) bioreactor at 37 °C with 900 rpm and the optical modules (product number and gain in parenthesis): scatter (E-OP-201, 3), riboflavin fluorescence (E-OP-227, 6), pH (E-OP-202, 7), O2 saturation (E-OP-203, 7), and NADH fluorescence (E-OP-405, 1).
	
2. Methods for processing the data: 
	25_Ecoli_2020_REDUCTION-1.csv was produced by reading 25-Ecoli_2020--2020-02-18-15-37-46.csv into the BioLection software and saving the data via Data Management > Save as... > 1 (No Reduction)
	
	Plot annotations in plot_annotations_acetate.csv and plot_annotations_glucose.csv were manually set during data analysis for better visibility.

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)
	
	BioLection (V.3.16.74.0, optional)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: visual comparions of data curves with replicates and previous experiments, data is unchanged

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig


DATA-SPECIFIC INFORMATION FOR: 25_Ecoli_2020_REDUCTION-1_layout.csv

1. Number of variables: 4

2. Number of cases/rows: 48

3. Variable List: 
	Ec_W3110Z1: used strain (Escherichia coli W3110Z1)
	
	Glc: concentration of glucose [C-mmol l^-1]
	
	Ace: concentration of acetate [C-mmol l^-1]
	
	aTc: concentration of anhydrotetracycline for indcution (NOT USED)

4. Missing data codes:

5. Specialized formats or other abbreviations used:


DATA-SPECIFIC INFORMATION FOR: 25_Ecoli_2020_REDUCTION-1.csv

1. Number of variables: 17

2. Number of cases/rows: 48

3. Variable List: 
	Biomass: backscattering/ 180 ° scattering value [AU]
	
	Riboflavine: riboflavin fluorescence [AU]
	
	Cali.pH(HP8): calibrated pH signal []
	
	pH(HP8): raw pH signal [AU]
	
	Cali.DO(Pst3): calibrated O2 concentration signal [% max]
	
	DO(Pst3): raw O2 concentration signal [AU]
	
	NADH - NADPH: pooled NADH and NADPH fluorescence signal [AU]
	
	Temperature: unknown, possibly measurement chamber temperature [°C]
	
	Temperature Down: unknown
	
	Temperature Water: unknown
	
	O2: O2 fraction in measurement chamber atmosphere [%]
	
	CO2: CO2 fraction in measurement chamber atmosphere [%] (NOT USED SINCE MODULE NOT INSTALLED)
	
	Shaker: shaking speed [rpm]
	
	Humidity: atmosphere humidity in measurement chamber [%rH]
	
	User Comment: manually added comments during the run (uses second time column)
	
	System Events: automatically noted events during the run (uses second time column)

4. Missing data codes:

5. Specialized formats or other abbreviations used: 


DATA-SPECIFIC INFORMATION FOR: plot_annotations_acetate.csv

1. Number of variables: 4

2. Number of cases/rows: 35

3. Variable List: 
	groups: acetate concentration for the annotated point [C-mmol l^-1]
	
	x: time position of the annotated point [h]
	
	pos: assigned position of the annotation in the plot (+: top, -: bottom)
	
	lab: numeric-label of the annotated point in the plot

4. Missing data codes: NA

5. Specialized formats or other abbreviations used:


DATA-SPECIFIC INFORMATION FOR: plot_annotations_glucose.csv

1. Number of variables: 4

2. Number of cases/rows: 28

3. Variable List: 
	groups: glucose concentration for the annotated point [C-mmol l^-1]
	
	x: time position of the annotated point [h]
	
	pos: assigned position of the annotation in the plot (+: top, -: bottom)
	
	lab: numeric-label of the annotated point in the plot

4. Missing data codes: NA

5. Specialized formats or other abbreviations used: