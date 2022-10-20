
#Universal Bioconductor package installation function
  install.bioc <- function(pkg){
    vers <- getRversion()
    if (vers < 4.1){
      "You need R >= 4.1"
    }else{
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }

#Install Bioconductor dependencies
bioc_pkgs <- c()
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.bioc(pkg)
    print("Installed!!!")
  }
}

#Install CRAN dependencies
cran_pkgs <- c("tidyverse","parallel")
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}
