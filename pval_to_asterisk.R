
# Function to convert p-value to significance asterisks

pval_to_asterisk <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return(NA)
}