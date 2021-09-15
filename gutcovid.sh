#!/bin/bash

chmod 777 deploy/*.sh
docker build deploy/ -t gutcovid:latest
docker rm -f gutcovid
docker run -it --rm --name gutcovid gutcovid:latest /bin/bash
