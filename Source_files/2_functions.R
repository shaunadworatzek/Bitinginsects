

# function to find medoid sequence in a cluster
find_medoid <- function(distmat, members) {
  mat <- as.matrix(distmat)
  members <- intersect(members, rownames(mat)) # only keep matching names
  submat <- mat[members, members, drop = FALSE]
  sums <- rowSums(submat)
  medoid <- names(which.min(sums))
  return(medoid)
}

#binary trait matrix

get_binary_trait_matrix <- function(seq_data) {
  
  mosmetadata <- kbimp_mosmetadata %>%
    filter(SampleID %in% names(seq_data)) %>%
    select(SampleID, ExactSite) %>%
    column_to_rownames(var = "SampleID")
  
  binary_matrix_trait <- model.matrix(~ ExactSite - 1, data = mosmetadata)
  
  binary_matrix_trait <- binary_matrix_trait %>% as.data.frame() %>%
    rownames_to_column(var = "SampleID") 
  
  return(binary_matrix_trait)
  
}


# Create a function to calculate hulls for each group
calculate_hull <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]  # Select convex hull points
}
