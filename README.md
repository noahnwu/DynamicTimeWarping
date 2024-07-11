# DynamicTimeWarping

Thanks for taking a look at our repository. For any questions that arise, please reach out to the corresponding author of 
our manuscript, Theresa Walunas. 

## Files in this repository: 

* ManuscriptPipeline.R
	* This file, when executed sequentially, runs the main pipeline described on a simulated dataset. It is divided into 4 primary parts
			1. Generation of the Simulated Dataset
			1. Execution of our implementation of the DTW algorithm and a comparison with a previously published R package 
			1. Clustering of patients based on the DTW distance with the K-Medoids Algorithm 
			1. A simplified longitudinal visualizaiton of patients' cachexia sequences over pseudotime
			
* DTW.R 
	* This contains the code for our implementation of the DTW algorithm and is used in the ManuscriptPipeline.R script
	* We developed our own code for DTW due to our future goal of utilizing other distance metrics with different data types

* AddIndividualIDs.R
	* Contains a function referenced in the ManuscriptPipeline.R script and is used for simple data manipulation
