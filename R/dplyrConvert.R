
dplyrConvert <- function(){

  data(taylor,package = "prostateCancerTaylor")
  library(tidyr)
  library(dplyr)
  library(Biobase)
  geoData <- taylor
  
  zscore <- function(x) t(scale(t(x)))
  
  iqrs <- apply(log2(exprs(geoData)), 1, IQR,na.rm=TRUE)
  
  taylor <- tbl_df(data.frame(Probe=featureNames(geoData),IQR=iqrs,log2(exprs(geoData)))) %>%
    gather(geo_accession,Expression,- c(Probe,IQR))
  
  taylor.z <- tbl_df(data.frame(Probe=featureNames(geoData),zscore(log2(exprs(geoData))))) %>%
    gather(geo_accession,Z,- Probe)
  
  fd <- tbl_df(fData(geoData)) %>% rename(Probe = ID)
  
  pd <- tbl_df(pData(geoData))
  

  taylor <- taylor %>% full_join(pd) %>% full_join(fd) %>% full_join(taylor.z)

  ##selecting most variable probe for each gene

    varProbes <- taylor %>% group_by(Gene) %>% 
    summarise(Probe = Probe[which.max(IQR)])
  
  taylor <- inner_join(taylor, varProbes,by="Probe") %>% rename(Gene = Gene.x) %>% select(-c(Gene.y))
  
}
