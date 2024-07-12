#This script runs through the DTW and Clustering Pipeline Utilized in 
#our manuscript 

#Here you will find 4 primary sections. Note that these must be run sequentially
#in order for them to function.

#1. Generation of a simulated dataset with three distinct groups of patients
#   based on their Fearon Consensus Criteria Levels

#2. Comparison of our DTW algorithm versus a previously published R package 
#   To ensure there were no errors in DTW development

#3. Clustering of patients using the K-Medoids Algorithm (Paritioning Around Medoids)
#   based on their DTW distance.

#4. Visualization of the clusters over pseudotime


#DTW package is that described in
#Toni Giorgino (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#The dtw Package. Journal of Statistical Software, 31(7), 1-24, doi:10.18637/jss.v031.i07.
required.packages = list("dtw", 
                         "tidyverse", 
                         "cluster")

sapply(required.packages, function(x) { 
  if(length(find.package(x)) == 0) { 
    install.packages(x)
    }
  })
#previously developed dtw package. function used here is called 'dtw'
#note that this was not the function used for our manuscript and is 
#included here only for the sake of comparison with our implementation of DTW
library(dtw)
#will use for some data manipulation
library(tidyverse)
#cluster package contains the K-Medoids Clustering Algorithm
#used in our manuscript
library(cluster)

#DTW algorithm used in our work. Function imported 
#from this script is named "DTW.Forrest" for clarity
#We developed our own implementation of DTW to allow for the use of 
#different distance metrics, a future direction in our work. However, 
#we first compared the output distance metrics from our algorithm 
#to a previously developed one to ensure that there were no errors in our 
#code
source("scripts/DTW.R")
#A data manipulation function
source("scripts/AddIndividualIDs.R")

######## 1. Generate Simulated Data ########

#Due to data privacy restrictions, we are unable to provide the exact dataset
#utilized in this work. However, a simulated dataset will be created 
#consisting of data similar to that used in the manuscript 
#for illustrative purposes. 

#All Cachexia Measures Calculated from 6-month Weight Loss and
#BMI resulted in one of three possible values according to the 2011 Fearon
#Consensus Criteria: 
#1 = No Cachexia 
#2 = Pre-Cachexia 
#3 = Cachexia 
#Sequences of 1s, 2s, and 3s were the final inputs into the DTW algorithm. 
FearonLevels = c(1, 2, 3)
set.seed(38212)
GroupSize = 100
MinSequenceLength = 2
MaxSequenceLength = 30
#We will artificially create three groups of patients in our simulated 
#data to show that these can be recapitulated following use of K-Medoids on 
#DTW distance measures.
#Generation of These Three Groups has been repeated manually to clearly 
#illustrate the differences in sequences
Mild = vector("list", length = GroupSize) %>%  
  lapply(., function(...) { 
    #Randomly Generated Length of Cachexia Measures
    SequenceLength = sample(MinSequenceLength:MaxSequenceLength, 
                            size = 1)
    #randomly generated sequence of 1s, 2s, and 3s
    PatientSequence = sample(FearonLevels, 
                             size = SequenceLength, 
                             #artificially creating patients with more or less
                             #severe measures
                             prob = c(0.8, 0.1, 0.1), 
                             replace = T)
    return(data.frame(TimePoint = 1:SequenceLength, 
                      CachexiaMeasures = PatientSequence))
    })

Moderate = vector("list", length = GroupSize) %>%  
  lapply(., function(...) { 
    #Randomly Generated Length of Cachexia Measures
    SequenceLength = sample(MinSequenceLength:MaxSequenceLength, 
                            size = 1)
    #randomly generated sequence of 1s, 2s, and 3s
    PatientSequence = sample(FearonLevels, 
                             size = SequenceLength, 
                             prob = c(0.1, 0.8, 0.1), 
                             replace = T)
    return(data.frame(TimePoint = 1:SequenceLength, 
                      CachexiaMeasures = PatientSequence))
  })

Severe = vector("list", length = GroupSize) %>%  
  lapply(., function(...) { 
    #Randomly Generated Length of Cachexia Measures
    SequenceLength = sample(MinSequenceLength:MaxSequenceLength, 
                            size = 1)
    #randomly generated sequence of 1s, 2s, and 3s
    PatientSequence = sample(FearonLevels, 
                             size = SequenceLength, 
                             prob = c(0.1, 0.1, 0.8), 
                             replace = T)
    return(data.frame(TimePoint = 1:SequenceLength, 
                      CachexiaMeasures = PatientSequence))
  })

#combined all groups into one dataframe
PatientData = bind_rows(
  AddIndividualIDs(Mild, "Mild"), 
  AddIndividualIDs(Moderate, "Moderate"), 
  AddIndividualIDs(Severe, "Severe")
)
#check that there are 3*GroupSize Patients in Our Dataset
length(unique(PatientData$PatientID)) == (3*GroupSize)

############## 2. Perform DTW #################

#Before proceeding with DTW, we must create a reference data frame that lists 
#all combinations of the 300 patients in our simulated dataset

Patients = unique(PatientData$PatientID)
Patient.Combos = combn(Patients, 2) %>% 
  t() %>% 
  as.data.frame() %>%  
  rename_with(.fn = ~c("Patient_1", "Patient_2"))


#Perform DTW with previously defined package (not used in this manuscript)
#And our code in the manuscript
DTW.Previous = vector("numeric", length = nrow(Patient.Combos))
DTW.Manuscript = vector("numeric", length = nrow(Patient.Combos))

#Brute Force Run DTW Across All Patient Combinations
#Took approximately 2 minutes on 80gb of RAM
for(i in 1:nrow(Patient.Combos)) { 
  #Get Sequence of Cachexia Measures of the two patients being compared
  Patient1 = PatientData$CachexiaMeasures[PatientData$PatientID == Patient.Combos$Patient_1[i]]
  Patient2 = PatientData$CachexiaMeasures[PatientData$PatientID == Patient.Combos$Patient_2[i]]
  
  #perform the previously published and our implementation of DTW
  DTW.Previous[i] = dtw(Patient1, Patient2)$normalizedDistance
  DTW.Manuscript[i] = DTW.Forrest(Patient1, Patient2)$normalized.distance
}


#compare the two distance measures to show that they are equivalent
plot(DTW.Previous, 
     DTW.Manuscript) 
cor(DTW.Manuscript, DTW.Previous) #1
sum(DTW.Manuscript == DTW.Previous) == length(DTW.Previous) #TRUE

############# 3. Cluster Patients Based on DTW Distance ##############

#We first must create a square distance matrix to implement k-medoids

PatientCombosDist = Patient.Combos %>% 
  mutate(Dist = DTW.Manuscript)


#Produces a symmetric Distance Matrix
PatientDistMatrix = PatientCombosDist %>% 
  bind_rows(
    #for the the symmetric triangle of the original distances
    PatientCombosDist %>% 
      rename_with(.fn = ~c("Patient_2", "Patient_1", "Dist")) %>%  
      select(Patient_1, Patient_2, Dist)
  ) %>%  
  bind_rows(
    #for the matrix diagonal. The distance of a patient sequence 
    #from itself will always be 0, which is shown below.
    data.frame(Patient_1 = Patients, 
               Patient_2 = Patients, 
               Dist = 0)
  ) %>%  
  arrange(Patient_1, Patient_2) %>%  
  pivot_wider(names_from = Patient_2, values_from = Dist) %>%  
  column_to_rownames("Patient_1")

#Simply to Show that the Normalized DTW Distance of two identical sequences 
#is 0
DTW.Forrest(PatientData$CachexiaMeasures[PatientData$PatientID == Patients[1]], 
            PatientData$CachexiaMeasures[PatientData$PatientID == Patients[1]])[["normalized.distance"]]

#Will try k = 2 to k= 10 and select the best solution based on 
#the maximum of the sillhouette value
KTry = 2:10
ClusterSolutions = vector("list", length(KTry))

for(i in 1:length(KTry)) { 
  cat("Trying K = ", KTry[i], "\n", sep = "")
 ClusterSolutions[[i]] = pam(as.dist(PatientDistMatrix), 
      diss = T, 
      k = KTry[i], 
      nstart = 20)
}

#Find the Best Cluster Number K by idenfiying the maximum sillhouette width
#With the specified parameters, Best k = 3
BestK = KTry[which.max(sapply(ClusterSolutions, function(x){
  x$silinfo$avg.width
}))]

BestCluster = pam(as.dist(PatientDistMatrix), 
                  diss = T, 
                  k = BestK, 
                  nstart = 20)

#Check cluster membership against artificially applied labels

BestCluster$clustering %>%  
  as.data.frame() %>%  
  rownames_to_column() %>%  
  rename_with(.fn = ~c("Patient", "Cluster")) %>%  
  mutate(PatientGroup = str_remove(Patient, "_\\d+$"), 
         Cluster = paste0("Cluster ", Cluster)) %>%  
  group_by(Cluster, PatientGroup) %>%  
  count() %>% 
  ungroup() %>% 
  group_by(PatientGroup) %>%  
  mutate(Percent = paste0((n/sum(n))*100, "%") ) %>%  
  select(-n) %>%
  pivot_wider(names_from = PatientGroup, 
              values_from = Percent) %>%  
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(is.na(.x), "0%", .x))) 

#note that cluster membership (for the most part) follows the labels applied to 
#patients at the beginning


###### 4. Visualize Cachexia Severity Over Time By Cluster #######
#We cannot entirely recreate our visualization. 
#However, we can compute the mean Fearon Level over the artificial time points
#Between our starting groups as well as our clusters to see their differences 
#over time

#Visualize by Cluster
PatientData %>%  
  left_join(BestCluster$clustering %>%  
              as.data.frame() %>%  
              rownames_to_column() %>%  
              rename_with(.fn = ~c("PatientID", "Cluster"))) %>%  
  group_by(TimePoint, Cluster) %>%  
  summarise(AvgFearon = mean(CachexiaMeasures), 
            SEFearon = sd(CachexiaMeasures)/sqrt(n()), 
            MaxFearon = AvgFearon + 1.96*SEFearon, 
            MinFearon = AvgFearon - 1.96*SEFearon) %>%  
  ggplot(., aes(TimePoint, AvgFearon))+
  geom_line(aes(color = factor(Cluster)))+
  geom_point(aes(color = factor(Cluster)))+
  geom_ribbon(aes(ymin = MinFearon, ymax = MaxFearon, 
                  fill = factor(Cluster)), 
              alpha = 0.25)+
  labs(color = "Cluster")+
  guides(fill = "none")+
  theme_bw()

#Visualize by Patient Group
PatientData %>%  
  mutate(PatientGroup = str_remove(PatientID, "_\\d+$")) %>%  
  group_by(TimePoint, PatientGroup) %>% 
  summarise(AvgFearon = mean(CachexiaMeasures), 
            SEFearon = sd(CachexiaMeasures)/sqrt(n()), 
            MaxFearon = AvgFearon + 1.96*SEFearon, 
            MinFearon = AvgFearon - 1.96*SEFearon) %>%  
  ggplot(., aes(TimePoint, AvgFearon))+
  geom_line(aes(color = PatientGroup))+
  geom_point(aes(color = PatientGroup))+
  geom_ribbon(aes(ymin = MinFearon, ymax = MaxFearon, 
                  fill = factor(PatientGroup)), 
              alpha = 0.25)+
  labs(color = "Cluster")+
  guides(fill = "none")+
  theme_bw()





