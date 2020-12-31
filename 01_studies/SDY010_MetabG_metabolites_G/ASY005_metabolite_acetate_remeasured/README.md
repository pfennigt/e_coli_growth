This ASY005_metabolite_acetate_remeasured README.md file was generated on 2020-12-30 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: SDY010_MetabG_metabolites_G/ASY005_metabolite_acetate_remeasured

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

3. Date of data collection (single date, range, approximate date): 2020-10-27

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
	20201027_metabolite_acetate_remeasured.csv: CLARIOstar acetate measurements of midrange BioLector samples
	
	20201027_metabolite_acetate_remeasured_METADATA.csv: metadata regarding the analyzed samples
	
	read_data.R: reads metabolite measurements and metadata data into R

2. Relationship between files, if important: 
	20201027_metabolite_acetate_remeasured_METADATA.csv contains the metadata of the samples in 20201027_metabolite_acetate_remeasured.csv
	
	read_data reads 20201027_metabolite_acetate_remeasured.csv and 20201027_metabolite_acetate_remeasured_METADATA.csv into R

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
	
	5 ml of an 1 OD cell solution were prepared by first measuring the optical density at 600 nm wavelength of a 1:4 dilution of culture with MilliQ. An appropriate amount of cell culture was spun down and the cell pellet was resuspended in 6 ml of fresh M9 medium with growth-supporting additives and antibiotics ("M9S", Tab. S2) by pipetting. Each 50 ml of M9S and M9G, as well as the 6 ml of 1 OD cell solution were used for culture praparation. A 48-well MTP-48-BOH well plate ("FlowerPlate"; m2p-labs GmbH, Baesweiler, Germany) was filled with 900 ul M9G medium per well resulting in a uniform glucose carbon-molarity of 119.88 C-mmol l^-1. Each inocoulated well was additionally filled with 100 ul of cell solution. The 3 blank wells contained additional 100 ul of M9S.
	
	The FlowerPlate was measured in a BioLector Pro (m2p-labs GmbH, Baesweiler, Germany) bioreactor at 37 °C with 900 rpm and the optical modules (product number and gain in parenthesis): scatter (E-OP-201, 3), riboflavin fluorescence (E-OP-227, 6), pH (E-OP-202, 7), O2 saturation (E-OP-203, 7), and NADH fluorescence (E-OP-405, 7).
	
	Every 15 min to 45 min, depending on how important a high resolution was deemed for the respective time range, a well of the FlowerPlate was sampled. The sample was spun down at 4 °C and both cell pellet and supernatant were frozen at -20 °C. Apart from the regular time series samples (numbered 1 to 19), five special samples were taken: A sample to assess the preculture glycogen contents (M) where 15 ml of preculture were spun down, two samples taken shortly before and after the first scatter drop (G1 and G2 respectively), and two samples taken shortly before and after the second scatter drop (G3 and G4 respectively).
	
	For acetate measurements the Sigma aldrich Acetate Colorimetric Assay Kit MAK086 was used. The acetate content in the midrange BioLector samples (9 to 15, G1, and G2) supernatant was estimated from ASY002_metabolite_acetate. The dilution of the samples for measurement was then chosen accordingly to keep the expected acetate content of the diluted samples within the kits detection range. Also, the dilutions were kept more similar, ranging between 1:25 and 1:50. Measurements were performed according to the kit instructions using in a clear 96 well plate and measuring the 450 nm absorbance in a CLARIOstar Plus (BMG LAB-TECH, Ortenberg, Germany) plate reader.
	
2. Methods for processing the data: 

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: visual comparions of data curves with replicates and previous experiments, data is unchanged

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig


DATA-SPECIFIC INFORMATION FOR: 20201027_metabolite_acetate_remeasured.csv

1. Number of variables: 4

2. Number of cases/rows: 14

3. Variable List: 
	Well: well of the 96 well plate
	
	Content: identifier for the well content  given by the CLARIOstar, corresponds to well_id in 20201027_metabolite_acetate_remeasured_METADATA.csv
	
	Raw Data (450): raw absorbance data at 450 nm [AU] (recommended for analysis)
	
	Linear regression fit based on Raw Data (450): estimated acetate content in the well from a linear model of the standard wells [nmol]

4. Missing data codes:

5. Specialized formats or other abbreviations used:


DATA-SPECIFIC INFORMATION FOR: 20201027_metabolite_acetate_remeasured_METADATA.csv

1. Number of variables: 8

2. Number of cases/rows: 14

3. Variable List: 
	well_id: identifier for the well content given by the CLARIOstar, corresponds to Content in 20201027_metabolite_acetate_remeasured.csv
	
	sample_id: identifier for the well content given during the experiment
	
	time_h: runtime of the biolector at sample collection [h]
	
	metabolite_type: measured metabolite
	
	sample_type: do the wells contain samples or standards
	
	standard_concentration_mmol_l: set concentration in standard wells [mmol l^-1]
	
	sample_dilution: dilution of the sample before measurement (diluted to x times concentration)
	
	well_volume_mul: liquid volume in the measured well [ul]
	
	wells: FlowerPlate wells sampled for the metabolite sample

4. Missing data codes:

5. Specialized formats or other abbreviations used: