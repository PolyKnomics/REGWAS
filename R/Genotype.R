Genotype <- R6::R6Class("Genotype", list(genotypefile = "", indexfile = "", chunks = "", setGenotype = function(genotypefile) {
  genotypefile <- normalizePath(genotypefile)
  indexfile <- paste0(genotypefile, ".bgi")
  if (!file.exists(indexfile)) {
    stop(paste0("could not find the index file: ", indexfile, " this can be generated with bgenix -index -g ", genotypefile))
  }
  self$genotypefile <- genotypefile
  self$indexfile <- indexfile
}, getNumberChunks = function() {
  return(nrow(self$chunks))
}, determ_chunks = function(nphenotypes = 1500, maxvar = 1e+05, mb = 512, region = NA) {
  bytespergenotype <- (3 * 8) + 4 # (3 values of 8 bytes+4 bytes for ploidy)
  bytes_per_variant <- bytespergenotype * nphenotypes
  memory_availble_byte <- mb * 1024^2
  chunksize <- min(maxvar, floor(memory_availble_byte / bytes_per_variant))
  regions <- data.frame(chromosome = character(0), start = integer(0), end = integer(0), stringsAsFactors = FALSE)
  a <- DBI::dbConnect(RSQLite::SQLite(), self$indexfile)

  position <- " "
  # count number of variants per chromosome
  if (is.na(region)) {
    chromosome_overview <- DBI::dbGetQuery(a, "SELECT chromosome,COUNT(*)FROM Variant GROUP BY  chromosome;")
  } else {
    chr_region <- unlist(strsplit(region, ":", 1))
    chr <- chr_region[1]
    if (length(chr_region) == 1) {
      # contain only a chromosome, not a start and end
      chromosome_overview <- DBI::dbGetQuery(a, paste0("SELECT chromosome,COUNT(*)FROM Variant WHERE chromosome==\"", chr, "\";"))
    } else {
      positionsplit <- unlist(strsplit(chr_region[2], "-", 1))
      region_start <- positionsplit[1]
      region_stop <- positionsplit[2]
      position <- paste0("AND position BETWEEN ", region_start, " AND ", region_stop, " ")
      print(paste0("SELECT chromosome,COUNT(*)FROM Variant WHERE chromosome==\"", chr, "\"", position, ";"))
      chromosome_overview <- DBI::dbGetQuery(a, paste0("SELECT chromosome,COUNT(*)FROM Variant WHERE chromosome==\"", chr, "\"", position, ";"))
    }
  }
  n_variants <- sum(chromosome_overview[, 2])
  n_chr <- nrow(chromosome_overview)
  print(paste0("found a total of ", n_variants, " variants in ", n_chr, " chromosomes"))

  if (n_variants == 0) {
    chromosome_overview <- DBI::dbGetQuery(a, "SELECT chromosome,COUNT(*)FROM Variant GROUP BY  chromosome;")

    stop(paste0("Could not find variants. The following chromosomes are availble:", paste(chromosome_overview$chromosome)))
  }


  for (i in seq_len(n_chr)) {
    chr_name <- chromosome_overview[i, 1]
    nchunks <- floor(chromosome_overview[i, 2] / chunksize)
    n_variants_chr <- chromosome_overview[i, 2] - 1
    if (nchunks < 2) {
      endpos <- c(n_variants_chr)
      startpos <- 0
    } else {
      endpos <- rev(seq(n_variants_chr, 0, by = -chunksize))[-1]
      startpos <- c(0, (endpos + 1)[1:(length(endpos) - 1)])
    }
    print(paste0("found in chr ", chr_name, "  ", n_variants_chr + 1, " variants"))
    for (j in seq_len(length(startpos))) {
      startsql <- paste0("SELECT position FROM Variant where chromosome==\"", chr_name, "\"", position, " LIMIT 1 OFFSET ", startpos[j], ";")
      endsql <- paste0("SELECT position FROM Variant where chromosome==\"", chr_name, "\"", position, " LIMIT 1 OFFSET ", endpos[j], ";")
      startq <- DBI::dbGetQuery(a, startsql)
      endq <- DBI::dbGetQuery(a, endsql)
      # check if the next start is not same position as current end position
      if (length(startpos) != j) {
        nextsql <- paste0("SELECT position FROM Variant where chromosome==\"", chr_name, "\"", position, " LIMIT 1 OFFSET ", endpos[j] + 1, ";")
        nextq <- DBI::dbGetQuery(a, nextsql)
        while (nextq$position == endq$position) {
          endpos[j] <- endpos[j] + 1
          startpos[j + 1] <- startpos[j + 1] + 1
          nextsql <- paste0("SELECT position FROM Variant where chromosome==\"", chr_name, "\"", position, " LIMIT 1 OFFSET ", endpos[j] + 1, ";")
          nextq <- DBI::dbGetQuery(a, nextsql)

          # add one
        }
      }

      regions <- rbind(regions, list(chromosome = chr_name, start = startq$position, end = endq$position), stringsAsFactors = FALSE)
    }
  }


  self$chunks <- regions
  print(regions)
  DBI::dbDisconnect(a)
}))
