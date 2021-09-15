FROM ubuntu:20.04
MAINTAINER Gut-COVID project "balazs.bohar@earlham.ac.uk"

USER root
RUN apt-get update

ADD install_base_layer.sh /
ADD install_python.sh /
ADD install_r.sh /
ADD install_r_packages.R /

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV PYTHONPATH /home/gutcovid/

RUN /bin/bash install_base_layer.sh
RUN /bin/bash install_python.sh
RUN /bin/bash install_r.sh
RUN Rscript install_r_packages.R

ADD /pipeline /home/gutcovid/
WORKDIR /home/gutcovid
