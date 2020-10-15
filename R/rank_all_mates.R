#' Rank all mates
#'
#' This function takes the output of `get_all_rankings` and ranks all possible matings so
#' that the sum of the score of all mate pairings is optimized. \cr
#' This function uses a permutation without replacement and makes a matrix of the size
#' of all possibilities. It is not recommended to use more than 10x10 samples.
#'
#' @param females Vector of sample IDs of the females to compare to the males. The Sample IDs
#' must match the ones found in the `DB`.
#' @param males Vector of sample IDs of the males to compare to the females. The Sample IDs
#' must match the ones found in the `DB`.
#' @param ranked_list List output from `get_all_rankings`
#' @param Output_head Integer value for the top number of mating pairings. Default:
#' top 5 mating pairings
#' @return Dataframe
#' @export
#' @examples
#' \dontrun{
#' all_mates <- rank_all_mates(DB = DBs, female = females,
#'                               males = males, type = "all_alleles",
#'                               bonus=NULL, weighted_alleles=NULL)
#' }
#' @export
#' @import dplyr
#' @import gtools
#' @import tidyverse

rank_all_mates <- function(females, males,
                           ranked_list, output_head = 5){
  # get permutation combinations
  ekk.3 <- gtools::permutations(
    n=length(males),r=length(females), v=males) %>%
    as.data.frame() %>%
    `colnames<-`(females)
  # Get socres
  tops <- sapply(1:ncol(ekk.3), function(i){
    df.t <- ranked_list[[colnames(ekk.3)[i]]]
    oo <- df.t[match(ekk.3[,i], df.t$male),
               "rank_score"]
    return(oo)
  }) %>% `colnames<-`(females) %>%
    data.frame()
  tops$sums <- rowSums(tops)
  # add scores
  outs <- ekk.3 %>%
    add_column(sums = tops$sums) %>%
    arrange(desc(sums)) %>%
    head(n=output_head)
  tops <- tops %>%
    arrange(desc(sums)) %>%
    head(n=output_head)
  # prepare for final output
  tips.1 <- outs %>%
    select(-sums) %>%
    t() %>%
    `colnames<-`(paste0("male_set_", colnames(.)))
  tips.2 <- tops %>%
    select(-sums) %>%
    t() %>%
    `colnames<-`(paste0("male_set_", colnames(.), "_score"))
  # combine for final output
  final_df <- data.frame(females = colnames(tops)[-ncol(tops)],
                         tips.1,tips.2)
  final_df <- final_df[,order(colnames(final_df))]
  return(final_df)
}
