"
Example:
find the test files at: tests/hypogeometric_test_files/
Rscript hypogeometric_enrichment_generic_test.R 10000 backround_input.tsv 200 foreground_input.tsv

args[1] is the total number of background objects (i.e. number of marbles)

args[2] is the background_input.tsv, space seperated in the format 'Annotation Number' , e.g.
cheese 2656
ham 2457
soup 2148

args[3] is the total number of foreground objects (i.e. number of white marbles)

args[4] is the foreground_input.tsv, space seperated in the format 'Annotation Number' , e.g.
cheese 76
ham 44
soup 31

E.g.
Input:
Rscript hypogeometric_enrichment_generic_test.R 10000 tests/hypogeometric_test_files/backround_input.tsv 200 tests/hypogeometric_test_files/foreground_input.tsv 

Output:
   Annot HgtTestMtCorrected FoldChange
1 cheese       0.0006937437  0.4307229
2    ham       2.4747251596 -0.1045991
3   soup       2.9620473074 -0.2783985

Author: Chris Rands 2018
"

args<-commandArgs(TRUE)

run_hypergeometric_test <- function(params)
{
    # Using the marbles analogy from wikipedia: https://en.wikipedia.org/wiki/Hypergeometric_distribution
    
    N <- params[1] # Number of marbles
    m <- params[2] # Number of white marbles
    n <- params[3] # Number of marbles drawn
    k <- params[4] # Number of white marbles drawn
    
    r <- N - m # Number of black marbles
    
    #Compute dhyper over range as want the probability that this many white marbles OR MORE would occur by chance
    #if marbles were drawn randomly (i.e. regardless of colour)
    
    hgt_result <- sum(dhyper(k:n, m, r, n, log = FALSE))
    
    # Note: this gives an identical ouput to these 2 approaches
    #hgt_result <- 1 - phyper(k-1, m, r, n, log = FALSE)    
    #hgt_result <- phyper(k-1, m, r, n, log = FALSE, lower.tail=FALSE)

    # Note: the raw result here is NOT corrected for multiple testing
    return(hgt_result)
}

wrapper_for_hypogeometric_test <- function(background_num, background, foreground_num, foreground)
{
    # Read in the data
    background_input_data <- read.table(background, sep = " ", blank.lines.skip = TRUE, col.names = c("Annot", "Number"))
    foreground_input_data <- read.table(foreground, sep=" ", blank.lines.skip = TRUE, col.names = c("Annot", "Number"))

    # Combine data frames
    merged_data <- merge(background_input_data,foreground_input_data, by = "Annot", suffixes = c("Background","Foreground"))
    
    # Calculate relevant parameters and run hypergeometric test
    N <- background_num # Number of marbles
    m <- foreground_num # Number of white marbles
    
    merged_data$HgtTest <- apply(merged_data[,c("NumberBackground","NumberForeground")], 1, function(x) run_hypergeometric_test(c(N,m,x[1],x[2])))
 
    # Apply Bonferroni mulitple testing correction; multiply by number of tests
    A <- nrow(merged_data) # Number of marbles drawn (number of tests)
    merged_data$HgtTestMtCorrected <- apply(merged_data["HgtTest"], 1, function(y) {y*A} )

    # Add fold change values and sort by lowested corrected hypergeometric test value    
    # Fold change is (B-A) / A where A is start value and B is final value
    # Here A is the proportion of all marbles that are white; and B is the proportion of drawn marbles that wre white
    merged_data$FoldChange <- apply(merged_data[,c("NumberBackground","NumberForeground")], 1, function(z) { ((z[2]/m)-(z[1]/N)) / (z[1]/N) } )
    merged_data_sorted <- merged_data[order(merged_data$HgtTestMtCorrected),]

    # Return final data frame
    return(merged_data_sorted[,c("Annot","HgtTestMtCorrected","FoldChange")])
}

wrapper_for_hypogeometric_test(as.integer(args[1]), args[2], as.integer(args[3]), args[4])
