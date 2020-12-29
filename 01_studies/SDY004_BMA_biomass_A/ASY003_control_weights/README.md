This ASY003_control_weights README.md file was generated on 2020-12-29 by Tobias Pfennig


GENERAL INFORMATION

1. Title of Dataset: SDY004_BMA_biomass_A/ASY003_control_weights

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

3. Date of data collection (single date, range, approximate date): 2020-06-05

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
	20200605_control_weights_GA.csv: manually noted drymass weights in mg for each hourly pooled sample
	
	20200605_control_weights_GA_METADATA.csv: metadata regarding the analyzes samples
	
	read_data.R: reads drymass data and metadata into R

2. Relationship between files, if important: 
	20200605_control_weights_GA_METADATA.csv contains the metadata of the samples in 20200605_control_weights_GA.csv
	
	read_data reads 20200605_control_weights_GA.csv and 20200605_control_weights_GA_METADATA.csv into R

3. Additional related data collected that was not included in the current data package: 

4. Are there multiple versions of the dataset? no
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? 
		ii. When was the file updated? 


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
	Information to used media and references can be found in ../../../06_manuscripts/Studying_Multiple_Growth_Phases_of_E_coli_on_Minimal_Medium_with_Experiment_and_Theory_Tobias_Pfennig.pdf
	
	5 ml of pure 22 mM (0.4 % w/v) glucose ("M9G", Tab. S2) or 66 mM (0.4 % w/v) acetate ("M9A", Tab. S2) medium were pulled through pre-dried and pre-weighed ReliaDisc Non-sterile Membrane Filtration Media (Ahlstrom-Munksjo, Helsinki, Finland) cellulose-acetate filters with 0.2 mm pore size filters using a vaccuum. The filters were subsequently dried at 60 °C in an emptied PEQLAB PerfectBlot hybridization oven (VWR International GmbH, Darmstadt, Germany) and weighed after 1, 2, 3, 5, 8, and 11 days after cooling to room temperature. During each weighing, filter weights were measured in 3 to 6 cycles on an Analytical balance ABP 100-4M (KERN & SOHN GmbH, Balingen, Germany) with 1 mg verification value. This served to evaluate if the protocol and media themselves produced significant changes in filter weights.
	
	CAUTION: The non-centered placement of the filters on the scale - as used in this experiment - has been shown to influence the measured weights of the filters. It is revised to redo the experiment.
	
2. Methods for processing the data: 

3. Instrument- or software-specific information needed to interpret the data: 
	R(4.0.2), packages: data.table (1.12.8), deSolve (1.28), ggplot2 (3.3.1), growthrates (0.8.2), lattice (0.20-41), mugro (0.0.1), platexpress (0.1), plyr (1.8.6), pracma (2.2.9), pspline (1.0-18),  stringr (1.4.0), tidyr (1.0.3)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: controlled environment (37 °C, 900 rpm shaking, gas transfer with kLa = 230 h^-1), atmosphere unchanged

6. Describe any quality-assurance procedures performed on the data: visual comparions of data curves with replicates and previous experiments, data is unchanged

7. People involved with sample collection, processing, analysis and/or submission: Rainer Machne, Tobias Pfennig, Jonas Burmester


DATA-SPECIFIC INFORMATION FOR: 20200605_control_weights_GA.csv

1. Number of variables: 35

2. Number of cases/rows: 7

3. Variable List: 
	sample_id: type of applied medium (A: M9A, G: M9G) and replicate number separated by a "_"
	
	pre_1: measured weight of the filter without medium application in the first measurement cycle [mg]
	
	pre_2: measured weight of the filter without medium application in the second measurement cycle [mg]
	
	pre_3: measured weight of the filter without medium application in the third measurement cycles [mg]
	
	pre_4: measured weight of the filter without medium application in the fourth measurement cycle [mg]
	
	pre_5: measured weight of the filter without medium application in the fifth measurement cycle [mg]
	
	pre_6: measured weight of the filter without medium application in the sixth measurement cycle [mg]
	
	post1_1: measured weight of the filter with medium application after one day of drying in the first measurement cycle [mg]
	
	post1_2: measured weight of the filter with medium application after one day of drying in the second measurement cycle [mg]
	
	post1_3: measured weight of the filter with medium application after one day of drying in the third measurement cycle [mg]
	
	post1_4: measured weight of the filter with medium application after one day of drying in the fourth measurement cycle [mg]
	
	post1_5: measured weight of the filter with medium application after one day of drying in the fifth measurement cycle [mg]
	
	post1_6: measured weight of the filter with medium application after one day of drying in the sixth measurement cycle [mg]
	
	post2_1: measured weight of the filter with medium application after two days of drying in the first measurement cycle [mg]
	
	post2_2: measured weight of the filter with medium application after two days of drying in the second measurement cycle [mg]
	
	post2_3: measured weight of the filter with medium application after two days of drying in the third measurement cycle [mg]
	
	post2_4: measured weight of the filter with medium application after two days of drying in the fourth measurement cycle [mg]
	
	post2_5: measured weight of the filter with medium application after two days of drying in the fifth measurement cycle [mg]
	
	post2_6: measured weight of the filter with medium application after two days of drying in the sixth measurement cycle [mg]
	
	post3_1: measured weight of the filter with medium application after three days of drying in the first measurement cycle [mg]
	
	post3_2: measured weight of the filter with medium application after three days of drying in the second measurement cycle [mg]
	
	post3_3: measured weight of the filter with medium application after three days of drying in the third measurement cycle [mg]
	
	post3_4: measured weight of the filter with medium application after three days of drying in the fourth measurement cycle [mg]
	
	post5_1: measured weight of the filter with medium application after five days of drying in the first measurement cycle [mg]
	
	post5_2: measured weight of the filter with medium application after five days of drying in the second measurement cycle [mg]
	
	post5_3: measured weight of the filter with medium application after five days of drying in the third measurement cycle [mg]
	
	post8_1: measured weight of the filter with medium application after eight days of drying in the first measurement cycle [mg]
	
	post8_2: measured weight of the filter with medium application after eight days of drying in the second measurement cycle [mg]
	
	post8_3: measured weight of the filter with medium application after eight days of drying in the third measurement cycle [mg]
	
	post8_4: measured weight of the filter with medium application after eight days of drying in the fourth measurement cycle [mg]
	
	post8_5: measured weight of the filter with medium application after eight days of drying in the fifth measurement cycle [mg]
	
	post8_6: measured weight of the filter with medium application after eight days of drying in the sixth measurement cycle [mg]
	
	post11_1: measured weight of the filter with medium application after eleven days of drying in the first measurement cycle [mg]
	
	post11_1: measured weight of the filter with medium application after eleven days of drying in the second measurement cycle [mg]
	
	post11_1: measured weight of the filter with medium application after eleven days of drying in the third measurement cycle [mg]
	
4. Missing data codes:

5. Specialized formats or other abbreviations used:
	pre: weighed before medium application
	
	post: weighed after medium application


DATA-SPECIFIC INFORMATION FOR: 20200303_drymass_weights_G_METADATA.csv

1. Number of variables: 3

2. Number of cases/rows: 7

3. Variable List: 
	sample_id: type of applied medium (A: M9A, G: M9G) and replicate number separated by a "_"
	
	time_h: runtime of the biolector at sample collection [h]
	
	volume_ml: sample volume applied to the cellulose acetate filter [ml]
	
	wells: flowerplate wells pooled to create the respective sample, separated by ";"
	
4. Missing data codes:

5. Specialized formats or other abbreviations used: