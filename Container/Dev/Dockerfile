FROM r-base:latest

RUN apt-get update -qq 	&& apt-get dist-upgrade -y

RUN apt-get update -qq && apt-get -y --no-install-recommends install curl libxml2-dev libsqlite-dev libmariadbd-dev libmariadbclient-dev  libpq-dev libssh2-1-dev unixodbc-dev  libsasl2-dev libssl-dev libcurl4-openssl-dev build-essential libxml2-dev libssl-dev

#RUN install2.r --error --deps TRUE devtools
RUN apt-get -y install r-cran-devtools
RUN install2.r --error --deps TRUE Rcpp RSQLite
RUN R --vanilla -e 'options(Ncpus = parallel::detectCores());install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" )'
RUN install2.r --error --deps TRUE coxme
RUN install2.r --error --deps TRUE covr
RUN install2.r --error --deps TRUE doParallel
RUN apt-get -y install git
