### Simulates from the Pseudpod model by calling the Java file
# DirName = directory name relative to .jar file location
# v = "false" no GUI prompt
# Duration = duration of simulation (100000) - set to default inside the function (might want to change it later)
# Pars = parameters list in the following order: ACTIN_FK, GI_SPEED, bB, dB, diffB, KM, sA, BASAL_A, diffA, dA

PseudoSim <- function( DirName, v = "false", Pars = c(0.0015, 0.07, 0.0028, 0.013, 0.045, 0.16, 7e-5, 0.1, 0.025, 0.02) ) {

    # Build command to use with system() function
    com <- paste(   "java -jar MC.jar v false DURATION 100000",
                    "DIR_OUT", DirName,
                    "ACTIN_FK", Pars[1], 
                    "GI_SPEED", Pars[2], 
                    "bB", Pars[3], 
                    "dB", Pars[4], 
                    "diffB", Pars[5],
                    "KM", Pars[6], 
                    "sA", Pars[7], 
                    "BASAL_A", Pars[8], 
                    "diffA", Pars[9], 
                    "dA", Pars[10], 
                    sep=" ")
    
    # Run MC.jar file with parameters Pars and save results in folder DirName
    system( com )
    
}


# Reads in the relevant output .csv file GI LI LA S X Y 

PseudoRead <- function ( Path = "MS1", Output = "GI" )
{

    # CSV file name in the format: LHC/MSi/Cell1/cellOutput.csv  MSi indicated by path
    fpath <- paste( "~/ownCloud/Phd Diana/Cside 2018/Simulator/Init_LHS/", Path, "/Cell1/", sep="")
    fname <- paste( "cell", Output, ".csv", sep="")
    
    # For GI read in one column file
    if ( Output == "GI")
        csvfile <- read.csv( paste( fpath, fname, sep=""), header = F)
    
    # For other outputs import csv file with variable number of columns files
    else {
        # detects maximum number of elements in the files (in each row) and takes the max of that
        Rdim <- max( count.fields( paste( fpath, "cellX.csv", sep=""), sep = ","))
        
        # forms column names N1, ... N Rnum - this is needed to force import missing values at the end of 
        # the beginning rows.
        CN <- paste0( "N", seq_len(Rdim))
        
        # read in relevant csv file
        csvfile <- read.csv( paste( fpath, fname, sep=""), header = F, col.names = CN, fill = T )
        
        # transform into numeric matrix file from list/data frame
        csvfile  <- matrix( unlist( csvfile,), ncol = dim( csvfile)[2] , byrow = FALSE)
        }
        
    
    # return csv file
    return( csvfile)
}