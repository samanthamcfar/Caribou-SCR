library(secr)
library(rgdal)
 
setwd("~/") #set working directory to location of txt files

traps <- read.traps('detectorfile.txt', detector = 'polygon') #read txt file of detector locations
caps <- read.capthist('capthist.txt', 'detectorfile.txt', fmt = 'XY', detector = 'polygon') #create capture history from txt file of captures and txt file of detectors
data.shape <-readOGR(dsn="./",layer="LittleSmoky_Range") #study area 


###mu 1###

b <- 1 #value of binomN
m <- 1 #level of clustering (mu). Repeat run_simpop below for varying levels of mu.
from= 1 #first file
fem_total = 110 #to - will run from 1 to 110 simulations. Set number of simulations here. 


###ET UP SIMULATION LOOP###
run_simpop <- function(textfile, position){
  tpop.sex.ls.1.1.0003. = sim.popn(D = 0.0003, popn = caps, core = traps, buffer = 15000, model2D = "cluster", buffertype = "rect", poly = data.shape, details = list(mu = m, hsigma = 1)) #set up simulated population
  tch.ls.1.1.0003. = sim.capthist(traps, tpop.sex.ls.1.1.0003., detectfn = 'HEX', detectpar = list(g0 = 0.03, sigma = 1600), noccasions = 3, binomN = b) #set up simulated capture history
  tch.disc.ls.1.1.0003. = discretize(tch.ls.1.1.0003., spacing = 1500, outputdetector = 'proximity') #discretize polygon data to proximity detector data
  fit0 = list(capthist = tch.disc.ls.1.1.0003., buffer = 15000, detectfn = 'HEX', method = 'BFGS', trace = FALSE, CL = TRUE) #set parameters for running secr model
  file1_data_name_export <- paste("./tfit.ls.1.1.0003.", position, ".Rdata", sep = "") #create file names for export
  data_frame1_name <- paste("tfit.ls.1.1.0003.", position, sep = "")
  assign(data_frame1_name, fit0)
  save(list = data_frame1_name, file=file1_data_name_export)
  
}

###RUN SIMULATIONS###
for(i in from:fem_total){
  tryCatch({
    text_file_name <- paste("tp.ls.1.1.0003.", i, ".txt", sep = "")
    run_simpop(text_file_name, i)
  }, error=function(e){})
}
testfiles <-  list.files(path = "~/", pattern = "tfit.ls.1.1") #list all simulated capture histories created above
for (i in 1:length(testfiles)) load(testfiles[[i]])
testfileslist <- mget(ls(pattern = "tfit.ls.1.1"))

fits <- par.secr.fit(testfileslist, ncores = 16, save.intermediate = TRUE, prefix = "par", LB = TRUE) #run secr models on list of capture histories. Set 'save.intermediate' to save all output files.


#calculate density estimate of each simulation using derived()
sims <- list.files(pattern = "tfit.ls.1.1.") #list files with your simulated capture histories from folder
for (i in 1:length(sims)) load(sims[[i]])
simslist <- mget(ls(pattern = "tfit.ls.1.1."))
secrlist <- secrlist(simslist)

sims.ls.1.derived <- par.derived(secrlist, ncores=16)
save(sims.ls.1.derived, file = "~/sims.ls.1.derived.RData")





