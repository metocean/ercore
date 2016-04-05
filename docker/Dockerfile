FROM metocean/ops-base
MAINTAINER Rosa Trancoso <rosa@metocean.co.nz>

RUN echo "-----------------Install Ercore-----------------" &&\
    git clone git@github.com:metocean/ercore.git $REPOS/ercore &&\
    cd $REPOS/ercore &&\
    pip install -e .


