# DynamicTimeWarping

Thanks for taking a look at our repository. For any questions that arise, please reach out to the corresponding author of 
our manuscript, Theresa Walunas.

The goal of our implementation of DTW on EHR data was to test two hypotheses: 

1. Different longitudinal progression patterns of cancer cachexia exist 
2. Differences in longituindal progression portend differential outcomes

In this repository, we will use the original DTW code on a simulated dataset to illustrate its uses in clustering.

## Files in this repository: 

* ManuscriptPipeline.R
	* This file, when executed sequentially, runs the main pipeline described on a simulated dataset. It is divided into 4 primary parts: 
	
			1. Generation of a Simulated Toy Dataset
			2. Execution of our implementation of the DTW algorithm and a comparison with a previously published R package 
			3. Clustering of patients based on the DTW distance with the K-Medoids Algorithm 
			4. A simplified longitudinal visualizaiton of patients' cachexia sequences over pseudotime
			
* scripts/DTW.R 
	* This contains the code for our implementation of the DTW algorithm and is used in the ManuscriptPipeline.R script
	* We developed our own code for DTW due to our future goal of utilizing other distance metrics with different data types

* scripts/AddIndividualIDs.R
	* Contains a function referenced in the ManuscriptPipeline.R script and is used for simple data manipulation



## How to use this code. 

Once you have successfully cloned this repo, run the ManuscriptPipeline.R in the IDE of your choice to view the results of the analysis. You are 
also welcome to use our DTW implementation for your own projects found in DTW.R (please cite our Cachexia Study!)

When running the ManuscriptPipeline.R script, please be sure that your working directory 
is the base directory of the repository (i.e. the directory in which the ManuscriptPipeline.R script is located)