FROM rocker/tidyverse:latest


RUN mkdir /usr/src/app/
RUN mkdir /usr/src/app/data
WORKDIR /usr/src/app
RUN apt-get update && apt-get upgrade -y
# Install R packages
RUN install2.r --error devtools
RUN install2.r --error argparser
RUN install2.r --error BiocManager
RUN apt-get update && apt-get install libz-dev
RUN apt-get -y install libbz2-dev


# RUN apt-get install libbz2-dev
ADD install_gitlab.R /Users/christianbouwens/Documents/resources/containers/dndscv-master/install_gitlab.R
RUN Rscript /Users/christianbouwens/Documents/resources/containers/dndscv-master/install_gitlab.R

# RUN R -e "devtools::install_github('im3sanger/dndscv')"
# RUN installGithub.r im3sanger/dndscv 
RUN apt-get -y install vim

COPY ./data/ /usr/src/app/data
COPY ./prp_dndscv.R /usr/src/app/prp_dndscv.R
COPY ./run_dndscv.R /usr/src/app/run_dndscv.R


CMD ["RSCRIPT", "/usr/src/app/run_dndscv.R"]
