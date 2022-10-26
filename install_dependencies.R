if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#Install CRAN dependencies
cran_pkgs <- c("tidyverse","parallel","gtools","combinat")
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0)
{
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst)
  {
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repos="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}
