FROM metocean/ops-base
MAINTAINER Rosa Trancoso <rosa@metocean.co.nz>

RUN pip install shapely

RUN echo "-----------------Install Ercore-----------------" &&\
    git clone git@github.com:metocean/ercore.git $REPOS/ercore &&\
    cd $REPOS/ercore &&\
    pip install -e . &&\
    ln -s $REPOS/ercore/erbio /usr/local/bin/ &&\
    ln -s $REPOS/ercore/eroil /usr/local/bin/ &&\
    ln -s $REPOS/ercore/ersar /usr/local/bin/ &&\
    ln -s $REPOS/ercore/erall /usr/local/bin/ &&\



