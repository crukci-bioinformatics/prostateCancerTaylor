
dplyrConvert <- function(){

  data(taylor,package = "prostateCancerTaylor")
  library(tidyr)
  library(dplyr)
  geoData <- taylor
  
  zscore <- function(x) t(scale(t(x)))
  
  taylor <- tbl_df(data.frame(Probe=featureNames(geoData),zscore(exprs(geoData)))) %>%
    gather(geo_accession,Expression,- Probe)
  
  fd <- tbl_df(fData(geoData)) %>% rename(Probe = ID)
  
  pd <- tbl_df(pData(geoData))
  

  taylor <- taylor %>% full_join(pd) %>% full_join(fd)
  taylor
}
