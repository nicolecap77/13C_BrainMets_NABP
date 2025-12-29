
# Generate regression formulas for all possible combinations of independent variables (indep_vars)
# including  (up to) two way interaction terms. Random effect terms are not included here; they are 
# listed in the main script as an input to the lme() function

generate_regression_formulas = function(indep_vars, dep_vars){
  n = length(indep_vars)
  formulas = list()
  
  make_interaction = function(x,y) paste(x,y,sep=":")
  
  for (k in 1:n){
    main_effects = combn(indep_vars, k, simplify=FALSE)
    
    for (main_set in main_effects){
      if (length(main_set)>1){
        interactions = combn(main_set, 2, FUN = function(x) make_interaction(x[1],x[2]), simplify = TRUE)}
      else {
        interactions = character(0)}
      
      interaction_subsets = unlist(lapply(0:length(interactions), function(i) {combn(interactions, i, simplify=FALSE)}),
                                   recursive=FALSE)
      
      for (int_set in interaction_subsets){
        terms=c(main_set, int_set)
        formula_txt=paste(dep_vars, paste(terms,collapse="+")) 
        formulas[[length(formulas)+1]] = as.formula(formula_txt)}
    }
  }
  return(formulas)
}

# # Example:
# # Define independent variables
# indep_vars= c("ParticipantType","SegID","Sex","Age")
# 
# # Define dependent variable, including transformation and tilda sign
# dep_vars = "log10(Bicarbonate) ~"
# 
# # Pass independent and dependent variables to generate_regression_formulas() function
# regression_formulas = generate_regression_formulas(indep_vars, dep_vars)