set.seed(5)
# generate basic phenotype information per sample
genotypeid <- sprintf("sample_%03d", c(1:500))
sex <- round(runif(500, 0, 1))
age <- round(runif(500, 18, 90))
basicpheno <- data.frame(genotypeid, sex, age)
rownames(basicpheno) <- genotypeid

# create for 99% of samples id suvival data
selected <- sample(genotypeid, as.integer(length(genotypeid) * 0.99))
tstart <- rep(0, length(selected))
tstop <- floor(runif(length(selected), min = 1, max = 11))
tstatus <- rep(0, length(selected))

# create for 99% of samples with survival data aditionalsuvival data
secondselected <- sample(selected, floor(length(selected) * 0.9))
tstatus[match(secondselected, selected)] <- 1
tstart2 <- tstop[match(secondselected, selected)]
tstop2 <- tstart2 + floor(runif(length(secondselected), min = 1, max = 5))
tstatus2 <- rep(0, length(secondselected))

# merge first and second round of survival data
sel_all <- c(selected, secondselected)
status <- c(tstatus, tstatus2)
tstart <- c(tstart, tstart2)
tstop <- c(tstop, tstop2)
length(tstop)
length(tstart)
length(status)
# selected basic phenodata in order of the survialdata merge this with survival data
sourcepheno <- data.frame(basicpheno[sel_all, ], tstart, tstop, status)

write.table(sourcepheno, file = "big_phenotype.tsv", row.names = F)
