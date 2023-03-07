genotypeid <- c("sample_001", "sample_002", "sample_003", "sample_004", "sample_005", "sample_005", "sample_006", "sample_007", "sample_008", "sample_008")

tstart <- c(0, 0, 0, 0, 0, 6, 0, 0, 0, 5)
tstop <- c(1, 4, 7, 10, 6, 10, 14, 18, 5, 18)
status <- c(0, 0, 0, 0, 1, 0, 0, 0, 1, 0)
sex <- c(0, 1, 0, 1, 1, 1, 0, 1, 0, 0)
age <- c(1, 2, 3, 4, 5, 5, 6, 7, 8, 8)
sourcepheno <- data.frame(genotypeid, tstart, tstop, status, sex, age)
write.table(sourcepheno, file = "small_phenotype.tsv", row.names = F)
