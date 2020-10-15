#' Make Database
#'
#' This function takes three inputs of CSV files and transforms them into
#' a RSQLite database. The three inputs are: \cr
#' 1. gtseq: Output of the GTseq pipeline \cr
#' 2. metadata: Metadata about the samples \cr
#' 3. allele_info: Information of the loci being used \cr
#'
#' @param gtseq Comma separated file that is the output of the GT-seq pipeline.
#' Sample IDs must match the Sample IDs in the `metadata` CSV file. The Allele IDs
#' from the column names must match the `site_id` found in the `allele_info` CSV file.
#' @param metadata Comma separated file that includes two columns: `Sample` and `Sex`.
#' The Sample IDs must match the Sample IDs from the `gtseq` CSV file. Other metadata
#' columns may be present.
#' @param allele_info Comma separated file with information about the loci of interest.
#' The `site_id` must match the Allele IDs from the GT-seq output. the `advantage` column
#' is required if the `type` `all_alleles` is specified.
#' @return An RSQLite database saved to memory
#' @export
#' @examples
#' \dontrun{
#' DBs <- make_database(gtseq = geno, metadata = meta_data,
#'                      allele_info = allele_info)
#' }
#' @import dplyr
#' @import RSQLite
#' @import stringr
#' @import tidyverse

make_database <- function(gtseq, metadata, allele_info){
  # Start Database
  DB_Test <- RSQLite::dbConnect(RSQLite::SQLite(), ":memory:")
  ## Get original DF
  dbWriteTable(DB_Test, "original_df", gtseq)
  ## Get the experimental data
  ### Need gtseq dataframe
  Experimental_data <-
    gtseq %>%
    select(Sample:IFI) %>%
    rename_at(vars(starts_with("X")),
              list(~(str_replace(., "X", "Percent")))) %>%
    rename_with(~gsub(".", "_", .x, fixed = TRUE))
  ### Check data types
  if (!all(sapply(Experimental_data, class) ==
           c("character", "integer", "integer",
             "numeric", "numeric", "numeric"))) {
    message("Not all columns are expected type. Do you have proper input data?")
  }
  ### Make table - Experimental_data
  s<- sprintf("create table %s(%s, primary key(%s))", "Experimental_data",
              paste(colnames(Experimental_data), collapse = ", "),
              colnames(Experimental_data)[1])
  dbSendStatement(DB_Test, s)
  dbWriteTable(DB_Test, "Experimental_data",
               Experimental_data, append = TRUE,
               row.names = FALSE,
  )
  ### Make table - MetaData
  s.m <- sprintf("create table %s(%s, primary key(%s))", "meta_data",
                 paste(colnames(metadata), collapse = ", "),
                 colnames(metadata)[1])
  dbSendStatement(DB_Test, s.m)
  dbWriteTable(DB_Test, "meta_data",
               metadata, append = TRUE,
               row.names = FALSE,
  )
  ### Make table - allele_info
  s.a <- sprintf("create table %s(%s, primary key(%s))", "allele_info",
                 paste(colnames(allele_info), collapse = ", "),
                 colnames(allele_info)[4])

  dbSendStatement(DB_Test, s.a)
  dbWriteTable(DB_Test, "allele_info",
               allele_info, append = TRUE,
               row.names = FALSE,
  )
  return(DB_Test)
}
