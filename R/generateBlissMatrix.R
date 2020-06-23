#' A function to generate bliss matrices for anchor / library drug combinations
#'
#' This function, given an input dataframe containing info on library dose, anchor dose, and inhibition, generates a list of bliss matrices for each combination. The matrix is intended to be later used as input to a heatmap function.
#' @param drug.data The input anchor / drug data.
#' @export
#' @examples
#' generateBlissMatrix(drug.data = my_data)

generateBlissMatrix <- function(drug.data) {
  library_drugs <- unique(drug.data$Compound)
  anchor_drugs <- unique(drug.data$Anchor)
  bliss_list <- list()
  for (drug in seq(library_drugs)) {
    # create bliss matrix for each lib / anchor combo & save in list
    # then call plotting function & save those results too
  }
}
