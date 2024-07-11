library(tidyverse)
AddIndividualIDs = function(SequenceList, 
                            GroupName, 
                            TimeColumn = "TimePoint") { 
  #Function that Concatenates Raw Cachexia 
  #Measure Sequences into a single dataframe
  #And adds patient ids
  return(
    SequenceList %>%  
      reduce(., bind_rows) %>%  
      mutate(PatientID = paste0(GroupName, "_", cumsum(!!sym(TimeColumn) == 1))) %>%  
      relocate(PatientID)
  )
  
}

