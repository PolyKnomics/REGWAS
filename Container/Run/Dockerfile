FROM clearlinux/r-base@sha256:8a59136e92d9b960927441fa4dea4d108a827e90aba54ab121e1dd174d1165db
ENV MKL_NUM_THREADS="N"
ENV MKL_DYNAMIC="FALSE"
ENV OMP_NUM_THREADS=1

RUN clrtrust generate
RUN swupd bundle-add curl which python3-basic R-extras

RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages("Rcpp", repos="https://mran.microsoft.com/snapshot/2020-11-30")'
RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" )'
RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages("survival", repos="https://mran.microsoft.com/snapshot/2020-11-30")'
RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages("RSQLite", repos="https://mran.microsoft.com/snapshot/2020-11-30")'
RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages("doParallel", repos="https://mran.microsoft.com/snapshot/2020-11-30")'
    
RUN curl http://code.enkre.net/bgen/tarball/release/bgen.tgz -o bgen.tar.gz && tar -xvf bgen.tar.gz && mv gavinband-* bgen && cd bgen && ./waf configure && ./waf && ./waf install && cd .. && rm -rf bgen.tar.gz bgen
#RUN swupd bundle-remove python3-basic 
RUN curl https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.6-Ubuntu16.04-x86_64.tgz -o qctool.tgz && tar -xvf qctool.tgz && cp /qctool*/qctool /usr/bin/qctool && rm -rf qctool*

RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages("coxme", repos="https://mran.microsoft.com/snapshot/2020-30-11")'

RUN swupd clean

CMD R
