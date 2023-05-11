#' Parentage analysis (using SNPPT.exe)
#' 
#' Function that simplifies using SNPPT.exe to calculate parentage based on SNP data
#' Requires that snppt.exe is installed in the working directory of this file!
#' More info at: https://swfsc.noaa.gov/textblock.aspx?Division=FED&ParentMenuId=54&id=16021
#'
#' Download snppit here: https://github.com/eriqande/snppit
#' Extract the zip-file, rename the "snppit-master" folder to "snppit", and place the folder in the same working directory as this script.
#' Then you're ready! Usually, to run snppit, you need to prepare a specially formatted genpop file with your offspring and parent SNP data,
#' and then run snppit via the command line. This r function saves you some time by doing the genpop formatting and command line stuff for you in R.
#' All you need to do is to supply the function with two datasets: parents and offspring, and the function will return the resulting parentage analysis ready to use in R.
#' The function takes two datasets as parameters, one for offspring and one for parents. \cr
#' Both datasets must be formated as folows:
#' \itemize{
#' \item Individuals as rows
#' \item SNP data on columns, one column pr SNP. Important: any column that is not named "ID", (or for parents, "sex" or "population"), will be interpreted as a SNP.
#' \item SNP genotype written numerically: 1=homozygous for allele a, 2=heterozygous, 3=homozygous for allele b
#' \item One column named "ID" that contains some unique ID for each individual,
#' \item Only for parents: One column named "sex" that's either "M" or "F"
#' \item Only for parents: One column named "population" that gives that individual's population
#' \item Only for parents: One column named "group" that gives the parents spawning group (foursome) must come after "sex"
#' }
#'
#'
#' @param df_offspring Data frame contaning offspring ID's and SNP genotypes (see below)
#' @param df_parents Data frame containing parent ID's, sex, population, and SNP genotypes (see above)
#' @param projectName Optional. A name that will be used for files created during the proces
#' @param overwrite Optional. If the script should overwrite any older snppit analysis (or read the old one again instead of re-doing it)
#' @param useGroups Optional. If parent grups should be used, only parents within the same group will be considered as parents together.
#' @export
#' @examples parentage_data <- run_snppt(offspring, parents, "Project_oct2019")
run_snppit <- function(df_offspring, df_parents, projectName="project1",overwrite=F, useGroups=F, full=F){
  require(glue)
  require(tidyverse)
  
  oldwd = getwd()
  setwd(paste(oldwd,"/snppit",sep=""))

  filename = glue("snppit_output_ParentageAssignments - {projectName}.csv")
  filename_full = glue("snppit_output_TrioPosteriors - {projectName}.csv")

  Sys.sleep(1)

  # First: check if this analysis has maybe already been done, if so, ask the user if she wants to skip it and load data from the previous run
  if (file.exists(filename) & overwrite==F){
    message("{filename} already exsist, loading that instead of doing new SNPPT run.")
    message("Set argument overwrite=T if you want to do a new analysis and overwrite old file.")
    if (!full){
      data_snppit = read_delim(filename)
      setwd(oldwd)
      return(data_snppit)
    }
    else{
      data_snppit = read_delim(filename_full)
      setwd(oldwd)
      return(data_snppit)
    }

  }

  check_columns <- function(tb,cols,name){
    missing <- cols[!cols %in% names(tb)]
    if (length(missing) > 0) stop(glue::glue("{name} missing columns {missing}"))
    }

  # Check that the parent set has columns "ID", "Sex", "population" and "group"
  if(useGroups==T)  check_columns(df_parents,c("ID","sex","population","group"),"df_parents")
  else check_columns(df_parents,c("ID","sex","population"),"df_parents")

  check_columns(df_offspring,c("ID"),"df_offspring")

  # Convert the genotype data in the dataset from normal numeric to snppt numeric type
  if (useGroups==T)
  {
  df_parents = df_parents %>%
    renameGenotypes(LUT=c("1"="1 1","2"="1 2","3"="2 2"), not_genotypes=c("ID","sex","population","group"))
  df_parents["group"][is.na(df_parents["group"])] <- "?"
  }
  else
  {
  df_parents = df_parents %>%
    renameGenotypes(LUT=c("1"="1 1","2"="1 2","3"="2 2"), not_genotypes=c("ID","sex","population"))
  }

  df_offspring = df_offspring %>%
    renameGenotypes(LUT=c("1"="1 1","2"="1 2","3"="2 2"), not_genotypes=c("ID"))

  # remove any column that is not equal between parents and offspring dataset (except population and sex and group)

  if (useGroups==T)  df_parents <- df_parents %>% select( "ID","population","sex","group", df_offspring %>% names() %>% one_of() )
  else               df_parents <- df_parents %>% select( "ID","population","sex", df_offspring %>% names() %>% one_of() )

  df_offspring <- df_offspring %>% select( "ID", df_parents %>% names() %>% one_of() )

  #Filename to use for snpptfile
  snpptfile_name = paste("snpptfile -", projectName)

  #Write SNPPIT file based on parent and offspring data
  message("Attempting to write SNPPT settings file...")
  SNPPITfile(snpptfile_name,df_parents, df_offspring, parentGroup=useGroups)
  message("SNPPT settings file written!")

  # Run snppit
  message("Starting SNPPT analysis...")
  location = paste("'",getwd(),"'",sep="") %>% str_replace_all("/","\\\\")
  system2("powershell", args=c("Set-location",location))
  system2("powershell", args=c(glue("./snppit -f '{snpptfile_name}.txt' --max-par-miss 70")))

  # Obtain snppit results and clean them
  message("Trying to obtain SNPPT results...")
  data_snppit = read.table("snppit_output_ParentageAssignments.txt", head=T, comment.char = "") %>%
    rename(ID_offspring=Kid, ID_pa=Pa, ID_ma=Ma, population=PopName)
  readr::write_delim(data_snppit, filename, delim="\t")

  data_snppit_full = read.table("snppit_output_TrioPosteriors.txt", head=T, comment.char = "")  %>%
    rename(ID_offspring=Kid, ID_pa=Pa, ID_ma=Ma)
  readr::write_delim(data_snppit_full, filename_full, delim="\t")

  # Return results
  message("results ready!")
  message(glue("Raw data can be found in {location}/{filename} and {filename_full}"))
  setwd(oldwd)
  if (full) return(data_snppit_full)
  else return(data_snppit)
}

#' SNPPITfile
#' @keywords internal
SNPPITfile = function(name,parents,offspring,parentSex=T,parentPopulation=T,parentGroup=T){
  # ASSUMES the following names for columns:
  # IDs:        "ID"         (string, no whitespace)
  # Sex:        "sex"        (F/M/?)
  # Population: "population" (string, no whitespace)
  # Family:     "group"     number or string, as long as they're unique
  # And no other columns other than loci!

  printLines = c()

  loci = names(parents)
  loci = loci [! loci %in% c("ID")]

  if (parentPopulation)
  {
    loci = loci[! loci %in% c("population")]
  }

  if (parentSex)
  {
    loci = loci[! loci %in% c("sex")]
    printLines = c(printLines, "POPCOLUMN_SEX")
  }

  if (parentGroup)
  {
    loci = loci[! loci %in% c("group")]
    printLines = c(printLines, "POPCOLUMN_SPAWN_GROUP")
  }

  printLines = c(paste("NUMLOCI", length(loci)),"MISSING_ALLELE *",printLines)



  for( i in loci)
  {
    printLines = c(printLines, paste(i, "0.005"))
  }

  snppitFile = file(paste(name,".txt",sep=""), "w")
  writeLines(printLines, snppitFile)

  if (parentPopulation)
  {
    for (i in unique(parents$population))
    {
      writeLines(paste("POP",i),snppitFile)

      write.table(parents[which(parents$population==i),] %>% select(-c("population")), snppitFile, sep="\t", quote=F, row.names = F, col.names=F)
    }
  }
  else
  {
    writeLines("POP Parentpool_0",snppitFile)
    write.table(parents, snppitFile, sep="\t", quote=F, row.names = F, col.names=F)
  }


  writeLines("OFFSPRING Offspring_pool ?",snppitFile)
  write.table(offspring, snppitFile, sep="\t", quote=F, row.names = F, col.names=F)

  close(snppitFile)
}



#' SNPPITfile
#' @keywords internal
format.snpptFile = function(name,df_parents,df_offspring,parentSex=T,parentPopulation=T,parentGroup=T){
  # ASSUMES the following names for columns:
  # IDs:        "ID"         (string, no whitespace)
  # Sex:        "sex"        (F/M/?)
  # Population: "population" (string, no whitespace)
  # Family:     "group"     number or string, as long as they're unique
  # And no other columns other than loci!

  # Initiate the text file
  snppitFile = file(paste(name,".txt",sep=""), "w")
  # Make a vector to contain all the lines in the text file
  printLines = c()

  # Create a list of loci - based on the column names in the parent dataset - then removing columns like ID, sex, population etc
  loci = names(df_parents)
  # clean non-snps from the loci list
  loci = loci[! loci %in% c("ID", "population", "sex", "group")]

  # Write POP_COLUMN_SEX
  if ("sex" %in% names(df_parents))
  {
    printLines = c(printLines, "POPCOLUMN_SEX")
  }

  # Write POPCOLUMN_SPAWN_GROUP
  if ("group" %in% names(df_parents))
  {
    printLines = c(printLines, "POPCOLUMN_SPAWN_GROUP")
  }

  # Write the number of loci.
  # Write the MISSING_ALLELE specifyer (we're using *)
  printLines = c(paste("NUMLOCI", length(loci)),"MISSING_ALLELE *", printLines)

  # Write in the ist of all the loci
  for( i in loci)
  {
    printLines = c(printLines, paste(i, "0.005"))
  }

  # Write the lines we have so far
  writeLines(printLines, snppitFile)

  if (! "population" %in% names(df_parents)) df_parents$population = "POP Parentpool_0"

  # For each population:
  for (i in unique(df_parents$population))
  {
    # Write the populatin name first
    writeLines(paste("POP",i),snppitFile)
    # Then write inn al the parents and their snps
    write.table(df_parents[which(parents$population==i),] %>% select(-c("population")), snppitFile, sep="\t", quote=F, row.names = F, col.names=F)
  }

  # Write the offspring pool (all the same pool)
  writeLines("OFFSPRING Offspring_pool ?",snppitFile)
  # write inn al the offspring
  write.table(df_offspring, snppitFile, sep="\t", quote=F, row.names = F, col.names=F)

  # close/save the file
  close(snppitFile)
}


#' Rename genotypes based on a lookup table
#'
#' In a dataframe, rename genotype columns
#' @keywords internal
renameGenotypes = function(dataframe, LUT, not_genotypes=c()) {
  for (i in names(dataframe %>% select(-c(not_genotypes)))) {
    dataframe <- dataframe %>% renameGenotype(i, LUT)
  }
  dataframe
}


#' renameGenotype
#' @keywords internal
renameGenotype = function(dataframe, column, LUT=c("1"="1 1","2"="1 2","3"="2 2")){
  genotype <- dataframe[[column]]

  col = LUT[genotype]
  col[is.na(col)] = "* *"


  dataframe[column] = col

  return(dataframe)
}
