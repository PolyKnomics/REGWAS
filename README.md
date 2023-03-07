# RE-GWAS

## Introduction
RE-GWAS is a tool for running recurrent-event genome-wide association
studies. The main program is written in R and can be installed as an R
package. The `Container` directory contains the necessary files to
wrap the R package into a Singularity/Apptainer container.

## Running a recurrent event GWAS
To run an recurrent event GWAS analysis using the container, the
following needs to be done:
- Copy the `regwas.sif` pacakge to a location where you will run the
  analysis (e.g. your PC or a compute cluster).
- Start an interactive session on that system. For example, on the
  Dutch national Cartesius super computer the following starts an
  interactive session of 60 minutes: `srun -t 60 -N 1 --pty bash -il`.
- Start R from within the container:  `singularity run regwas.sif` or
  `./regwas.sif`
- Load the R package:
  ```R
  library(REGWAS)
  ```
- Create an object with the genotypes in `bgen` format as first
  argument and the phenotype file as second argument. Note, the `bgen`
  file has to be indexed!
- The following function returns a `Regwas` object and performs
  several checks on load:
  ```R
  regwas <- Regwas$new("example.8bits.bgen", "phenofile.tsv")
  ```
  It is possible to change the default algorithm by calling a function
  of the object.
  - robust Cox PH: `regwas$set_default_model()` (this is the default)
    and doesn't need to be called explicitly)
   - frailty with gamma distribution: `regwas$set_frailty.gamma()`
   - frailty with gaussian distribution: `regwas$set_frailty.gaussian()`
   - coxme `regwas$set_coxme()`

  To change the genetic model call one of the following functions:
  - additive model: `regwas$setModel.additive()`
  - 2df model: `regwas$setModel.2df()`

  The `FlipMAF()` function can be used to flip the reference and effect
  allele according to the Minor Allele Frequency, such that the minor
  allele is the effect/predictor allele:
  ```R
  regwas$setFlipMaf(TRUE)
  ```
  This functionality can be turned like this:
  ```
  regwas$setFlipMaf(FALSE)
  ```
  This function only affects the 2df model. It has no effect on other
  genetic models.

  The current settings can be printed using:
  ```R
  print(regwas)
  ```
- Next, run the regression analysis on all chunks using the
  `regwasworker()` function. This function has two arguments, the
  `Regwas` object and the number of the chunk that needs to be
  processed. It generates a `data.frame` with the results as output.
  The number of chunks in the `Regwas` object can be found using the
  `regwas$g$getNumberChunks()` function. The
  `regwas$gwaswriter$write()` function is used to store the chunks. It
  has two arguments: a `data.frame` with the results (i.e. the output
  of `regwasworker()`) and the chunk number, so it can later put those
  in the correct order. This puts all of the above in a single piece
  of code:
   ```R
   for (chunknr in c(1:regwas$g$getNumberChunks())){
       regwas$gwaswriter$write(regwasworker(regwas, chunknr))
   }
   ```
- For parallel execution:
  ```R
  # Load package for parallel computation
  library(doParallel)
  # Detect the number of cores or set manually
  ncores <- detectCores()
  # Register the cores to use
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  # Excute in parallel
  foreach (chunknr=c(1:regwas$g$getNumberChunks()), .packages='REGWAS') %dopar% {
      regwas$gwaswriter$write(regwasworker(regwas, chunknr))
  }
  ```

Once the regression has been done on all chunks, the chunks can be
combined into a single file. The file name can be queried using
`regwas$gwaswriter$final_filename`. The function that merges the
temporary files is called `regwas$gwaswriter$mergeresults(nchunks)`,
where `nchunks` is the number of chunks (obtained from
`regwas$g$getNumberChunks()`).

It is also possible to call `qctools`, `bgenix`, etc. in the
container. To start a shell in the container, run `singularity shell
regwas.sif`, after which `R`, `qctools` and `bgenix` can be called.

# Install as R package
This packages has some dependencies. Most of them are availble on
CRAN, only the `rbgen` packages needs to be installed from another
source. Installing the dependencies can be done from within R using:
```R
install.packages( c("RSQLite", "Rcpp", "R6", "DBI", "coxme") )
install.packages("http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz",
                 repos = NULL, type = "source" )
```
Then copy the `REGWAS.tar.gz` packages to the current working direcory
and execute from within R:
```R
install.packages("REGWAS.tar.gz", repos = NULL, type="source")
```

## Development tips
Packages required for testing and development:
```R
install.packages(c("lintr", "testthat", "styler"))
```

To load the code during development run
```R
devtools::load_all()
```
To run the test:
```R
devtools::test()
```
Clean up stying issues with (best done in a separate commit to differentiate
'beauty' changes from functional changes
```
styler:::style_active_pkg()
```

Use lintr to detect errors:
```
lintr:::addin_lint_package()
```

## License
The code written for this project is licensed under the MIT license
(see the `LICENSE` file). Note, however, that this package depends on
other packages/libraries, which may have a different license.

## Acknowledgements
The development of this package was commissioned by Dr. Tessel
Galesloot, Department of Health Evidence, Radboud University Medical
Centre, Nijmegen, The Netherlands. The development was done by
PolyKnomics BV, The Netherlands (https://www.polyknomics.com).
