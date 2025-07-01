##########################################  MASTER FILE TO GENERATE BIVARIATE DATA AND FIT MODELS ########################################### 
########## This file generates all the data underlying results in the 'Simulation' section of the paper.                           ##########
########## Results are saved in the 'Data' folder to then be used by the plotting function to generate all figures for simulation. ##########
#############################################################################################################################################

library(callr)

# Set up a list of input configurations you wish to run
input_list <- list(
  list(Bt_mode=TRUE , a=.5*1:5    , b=.5*1:5        ,c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9), mu1=1          , mu2=2                  , n=1000,dist="NO"),
  list(Bt_mode=FALSE, a=.5*1:5    , b=.5*1:5        ,c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9), mu1=1          , mu2=2                  , n=1000,dist="NO")
  #,
  #list(Bt_mode=TRUE , a=NA        , b=c(.2,.5,1,2,5),c=c(.2,.5,1,2,5)               , mu1=c(.5,1,2,5), mu2=c(.5,1,2,5)        , n=1000,dist="PO"),
  #list(Bt_mode=FALSE, a=NA        , b=c(.2,.5,1,2,5),c=c(.2,.5,1,2,5)               , mu1=c(.5,1,2,5), mu2=c(.5,1,2,5)        , n=1000,dist="PO"),
  #list(Bt_mode=TRUE , a=.1+.1*1:20, b=.1+.1*1:20    ,c=NA                           , mu1=10         , mu2=12                 , n=1000, dist="GA"),
  #list(Bt_mode=FALSE, a=.1+.1*1:20, b=.1+.1*1:20    ,c=NA                           , mu1=10         , mu2=12                 , n=1000, dist="GA"),
  #list(Bt_mode=TRUE , a=NA        , b=NA            ,c=c(.1,.25,.5,.75,.9), mu1=c(.1,.25,.5,.75,.9)  , mu2=c(.1,.25,.5,.75,.9), n=1000,dist="LO"),
  #list(Bt_mode=FALSE, a=NA        , b=NA            ,c=c(.1,.25,.5,.75,.9), mu1=c(.1,.25,.5,.75,.9)  , mu2=c(.1,.25,.5,.75,.9), n=1000,dist="LO")
  
  # Add more configurations as needed
)

##

# Get the script lines (excluding inputs)
script_lines <- readLines("Simulations - 1 - Generate and Fit Models.R")
script_to_run <- script_lines[12:length(script_lines)] # lines 12 onwards, this doesn't actually matter as previous lines are just comments
# 
# #Iterate through the list of input configurations and save results to /Data/ folder
# for(inputs in input_list) {
#   
#   # Reset environment except for 'inputs', 'input_list', etc.
#   rm(list = setdiff(ls(), c("inputs", "input_list", "script_lines", "script_to_run")))
#   gc()
#   
#   # Assign each input variable into the global environment
#   list2env(inputs, envir = .GlobalEnv)
#   
#   # Optionally: print which configuration is running
#   print(inputs)
#   
#   # Evaluate the rest of the original script
#   eval(parse(text = script_to_run))
#   
#   # Optionally: move/summarise/save output files, etc.
# }

###Iterate over the sets of inputs, calling a new R session each time to ensure isolation of results and avoid package interaction errors

for (inputs in input_list) {
  p <- callr::r_bg(
    function(inputs) {
      list2env(inputs, envir = .GlobalEnv)
      script_lines <- readLines("Simulations - 1 - Generate and Fit Models.R")
      eval(parse(text = script_lines[12:length(script_lines)]))
    },
    args = list(inputs = inputs),
    supervise = TRUE,
    stdout = "|",  # Pipe output so it can be polled
    stderr = "|"
  )
  
  # Poll and print output live
  while (!p$is_alive()) Sys.sleep(0.1)
  while (p$is_alive()) {
    out <- p$read_output_lines()
    if (length(out)) cat(out, sep = "\n")
    err <- p$read_error_lines()
    if (length(err)) cat(err, sep = "\n")
    Sys.sleep(0.2)
  }
  # Print remaining output after process ends
  out <- p$read_all_output()
  if (length(out)) cat(out, sep = "\n")
  err <- p$read_all_error()
  if (length(err)) cat(err, sep = "\n")
}