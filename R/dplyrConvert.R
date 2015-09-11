
dplyrConvert <- function(){

  data(taylor,package = "prostateCancerTaylor")
  library(tidyr)
  library(dplyr)
  geoData <- taylor
  
  zscore <- function(x) t(scale(t(x)))
  
  iqrs <- apply(log2(exprs(geoData)), 1, IQR,na.rm=TRUE)
  
  taylor <- tbl_df(data.frame(Probe=featureNames(geoData),IQR=iqrs,zscore(exprs(geoData)))) %>%
    gather(geo_accession,Expression,- c(Probe,IQR))
  
  fd <- tbl_df(fData(geoData)) %>% rename(Probe = ID)
  
  pd <- tbl_df(pData(geoData))
  

  taylor <- taylor %>% full_join(pd) %>% full_join(fd)
  taylor
}
