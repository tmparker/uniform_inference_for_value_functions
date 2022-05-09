# Check the parts directory for missing runs, since sometimes they get dropped
# on the big computer.

ss <- 500
part <- 1:500

# Find which jobs have been dropped.
parts <- paste0("ss", ss, "part", sprintf("%03d", part), ".rda")
done <- match(list.files("./parts", pattern = paste0("ss", ss, "*")), parts)
done <- na.omit(done)
redo <- part[-done]
redo <- redo + 500 * (ss == 500) + 1000 * (ss == 1000) - 1

# Print the job numbers for those that remain to be done in a way that can be
# easily pasted into the run_simulations.sh script.
jmp <- which(diff(redo) != 1)

if (length(jmp) == 0) {
  if (length(redo) == 0) {
    for_slurm <- "all done!\n"
  } else {
    for_slurm <- paste0(redo[1], "-", redo[length(redo)])
  }
}

if (length(jmp) > 0) {
  if (jmp[length(jmp)] != redo[length(redo)]) {
    jmp <- c(jmp, length(redo))
  }
  contig <- vector(mode = "list", length = length(jmp))
  contig[[1]] <- seq(redo[1], redo[1] + jmp[1] - 1)
  for (i in 2:length(jmp)) {
    contig[[i]] <- seq(redo[jmp[i - 1] + 1], redo[jmp[i]])
  }
  fun <- function(s) {
    if (length(s) == 1) { paste0(s, ",", collapse = "")
    } else { paste0(s[1], "-", s[length(s)], ",", collapse = "") }
  }
  all_one <- paste(sapply(contig, fun), sep = "", collapse = "")
  for_slurm <- substr(all_one, 1, nchar(all_one) - 1) # no comma
}

cat(for_slurm)

