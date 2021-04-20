#Merge Picard output from "CollectRnaSeqMetrics"

#Identify all files
todolist <- list.files(pattern="*.txt")

mydata<- as.data.frame(matrix(nrow=1, ncol=31)) #creates blank matrix to bind all other data to

#Read tabular portion of output file to get column names
mydata_names <- readLines(file(todolist[1]), 7)[7]
mydata_names <- unlist(strsplit(mydata_names, split = "\t")[1])

names(mydata) <- c("Library", mydata_names)

#Setup a holding dataframe for histogram data (normalized coverage at percentiles of transcript length)
histogram_data <- as.data.frame(matrix(nrow=1, ncol=101))
names(histogram_data) <- c("Library", paste0("Position_", c(1:100)))

for (i in todolist) {
  #Make the name less messy
  library_name <- gsub("_+.*$", "", i) 
  con <- file(i)
  
  #Read data from appropriate lines - Picard-derived QC metrics
  workingData <- readLines(con, 8)[8]
  workingData <- unlist(strsplit(workingData, split = "\t")[1])
  
  #Merge data to holding dataframe
  workingData <- c(library_name, workingData)
  mydata <- rbind(mydata, workingData)
  
  #Read and merge histrogram data
  workinghistogramData <- readLines(con, 112)[12:112]
  workinghistogramData <- as.numeric(gsub(".*\\\t", "", workinghistogramData))
  workinghistogramData <- c(library_name, workinghistogramData)
  
  histogram_data <- rbind(histogram_data, workinghistogramData)
}

#Remove "placeholders" and save data

mydata <- mydata[-1,-c(29:31)]
write.csv(mydata,file="Picard_human_summary_output.csv")

histogram_data <- histogram_data[-1,]
write.csv(histogram_data,file="Picard_human_histogram_output.csv")
