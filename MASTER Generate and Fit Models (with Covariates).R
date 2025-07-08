## This is a test script to run the models with covariates

library(callr)

# Set up a list of input configurations you wish to run
input_list <- list(
  #list(Bt_mode=FALSE, a=c(1,2)    , b=.1+.4*1:2    ,c=NA                           , mu1=1         , mu2=2                 , n=1000,dist="GA",x1=1,x2=1),
  #list(Bt_mode=FALSE, a=1         , b=2            ,c=c(0.25,0.5,0.75), mu1=1          , mu2=2                  , n=1000,dist="NO",x1=1,x2=1),
  #list(Bt_mode=FALSE, a=NA        , b=2            ,c=1               , mu1=1, mu2=1    , n=1000,dist="PO",x1=1,x2=1),
  #list(Bt_mode=FALSE, a=NA        , b=NA           ,c=.75, mu1=.25  , mu2=.75, n=1000,dist="LO",x1=1,x2=1)
  list(Bt_mode=FALSE, a=c(.5,2.5)    , b=c(.5,2.5)        ,c=c(.25,.75), mu1=1          , mu2=2                  , n=1000,dist="NO",x1=1,x2=1),
  list(Bt_mode=FALSE, a=NA        , b=c(.2,5),c=c(.2,5)               , mu1=1, mu2=2        , n=1000,dist="PO",x1=1,x2=1),
  list(Bt_mode=FALSE, a=c(.5,1,2), b=c(.5,1,2)    ,c=NA                           , mu1=1         , mu2=2                 , n=1000,dist="GA",x1=1,x2=1),
  list(Bt_mode=FALSE, a=NA        , b=NA            ,c=c(.25,.75), mu1=c(.25,.5)  , mu2=.75, n=1000,dist="LO",x1=1,x2=1)
)

for (inputs in input_list) {
  p <- callr::r_bg(
    function(inputs) {
      list2env(inputs, envir = .GlobalEnv)
      script_lines <- readLines("Simulations - 1 - Generate and Fit Models (with Covariates).R")
      eval(parse(text = script_lines[11:length(script_lines)]))
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

