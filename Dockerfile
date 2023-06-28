FROM rocker/tidyverse:latest


RUN mkdir /usr/src/app/
RUN mkdir /usr/src/app/data
WORKDIR /usr/src/app

# Install R packages
RUN install2.r --error devtools
RUN install2.r --error argparser

RUN installGithub.r im3sanger/dndscv 

COPY ./data/ /usr/src/app/data
COPY ./prp_dndscv.R /usr/src/app/prp_dndscv.R
COPY ./run_dndscv.R /usr/src/app/run_dndscv.R


CMD ["RSCRIPT", "/usr/src/app/run_dndscv.R"]

