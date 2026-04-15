library(seraphim)

# Set working directory
setwd("/Users/claramal/Desktop/Moncla_Lab/H5N5-mini-project/SERAPHIM/")
getwd()

# Step 1 - extract spatio-temporal information from trees

localTreesDirectory = "Extracted_trees"
nberOfTreesToSample = 851 
mostRecentSamplingDatum = 2025.8575342465754

  # Use the "postTreeExtractions" function (advised)
  
trees = readAnnotatedNexus("../BEAST/cont-phyl_dta/results/A6_concat_genome_DTA-contphyl_2025-12-28/A6_concat_genome_skyride_DTA-contphyl_2026-01-06_851combined_regionreplaced.trees")
# This takes a bit of time
dir.create(localTreesDirectory, showWarnings=F)
for (i in 1:length(trees))
{
  tab = postTreeExtractions(post_tre=trees[[i]], mostRecentSamplingDatum)
  write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
}
# This also takes a while, and doesn't tell you while its going, but you can check the folder

# Step 2 - estimate dispersal statistics

nberOfExtractionFiles = 851
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = FALSE
outputName = "./spread_statistics/A6_concat-genome_skyride_cont-phyl-dta_stats"
nberOfCores = 1
slidingWindow = 1

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, 
                 showingPlots, outputName, nberOfCores, slidingWindow)
# 3 extraction files takes 10 minutes...