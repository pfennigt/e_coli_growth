This SDY006_InocG_inoc_gradient_G README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: SDY006_InocG_inoc_gradient_G

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

3. Date of data collection (single date, range, approximate date): 2020-06-06

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
	ASY001_biolector_data: folder with BioLector Pro measurement data and layout file
	
	ASY001_biolector_data_OTscript.py: OT-2 script used for pipetting of the measured FlowerPlate

2. Relationship between files, if important: 
	ASY001_biolector_data BioLector Pro measurement data was obtained from measuring the FlowerPlate resulting from the OT-2 pipetting given by ASY001_biolector_data_OTscript.py

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
	
	Liquid cultures were created by picking a colony from the source-plate and transferring it into 5 ml of minimal M9 medium with growth-supporting supplements, antibiotics, and 22 mM (0.4 % w/v) glucose ("M9G", Tab. S2). Yeast extract was also added to the medium to a concentration of 1 g l^-1. This liquid culture was then incubated at 250 rpm and 37 °C. Every two to seven days, 50 ul of this liquid culture were transferred into 5 ml of fresh M9G medium and continued to be incubated. One day before the measurement the culture was renewed as explained before.
	
	5 ml of an 1 OD cell solution were prepared by first measuring the optical density at 600 nm wavelength of a 1:4 dilution of culture with MilliQ. An appropriate amount of cell culture was spun down and the cell pellet was resuspended in 10 ml of M9G by pipetting. 50 ml of M9 medium with growth-supporting additives and antibiotics ("M9S", Tab. S2) and 10 ml M9G, as well as the 10 ml of 1 OD cell solution were placed in an OT-2 pipetting robot with P300 pipette module (Opentrons, New York, USA). A 48-well MTP-48-BOH well plate ("FlowerPlate"; m2p-labs GmbH, Baesweiler, Germany) was filled using the pipetting script ASY001_biolector_data_OTscript.py. The layout included 15-step gradients of inoculation OD ranging from 0.02 to 0.5 and a uniform glucose carbon-molarity of 119.88 C-mmol l^-1. Each incoulated well was filled with 100 ul of M9S and 900 ul of a mixture of M9G and the cell solution with varying ratio. The 3 blank wells contained 1 ml of 119.88 C-mmol l^-1 medium.
	
	The FlowerPlate was measured in a BioLector Pro (m2p-labs GmbH, Baesweiler, Germany) bioreactor at 37 °C with 900 rpm and the optical modules (product number and gain in parenthesis): scatter (E-OP-201, 3), riboflavin fluorescence (E-OP-227, 6), pH (E-OP-202, 7), O2 saturation (E-OP-203, 7), and NADH fluorescence (E-OP-405, 7).
	
2. Methods for processing the data: 

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)
	
	BioLection (V.3.16.74.0, optional)
	
	Opentrons App for use with OT-2 robot (3.18.1)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: visual comparions of data curves with replicates and previous experiments, data is unchanged

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig