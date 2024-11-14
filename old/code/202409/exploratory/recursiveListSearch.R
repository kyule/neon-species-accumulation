results <- list(
  A = list(x = list(out = 1, other = 2)),
  B = list(y = list(out = 3, other = 4)),
  C = list(z = list(out = 5, other = 6)),
  D = list(x = list(not_out = 10, other = 20))
)


find_out <- function(lst) {
    # Check if the list has an element named "out"
    if ("out" %in% names(lst)) {
      return(lst$out)
    } else {
      # Recursively search all elements of the list
      return(unlist(lapply(lst, find_out)))
    }
  return(NULL)
}

# Extracting all "out" items from the nested list
out_items <- find_out(results)

print(out_items)