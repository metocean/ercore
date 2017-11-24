FROM metocean/ops-base
MAINTAINER Simon Weppe <simon@metocean.co.nz>
RUN pip install shapely

RUN echo "-----------------Install Ercore-----------------" &&\
    git clone -b ercore_opensrc https://github.com/metocean/ercore.git $REPOS/ercore &&\
    cd $REPOS/ercore &&\
    pip install -e . &&\
    ln -s $REPOS/ercore/erbio /usr/local/bin/ &&\
    ln -s $REPOS/ercore/eroil /usr/local/bin/ &&\
    ln -s $REPOS/ercore/ersar /usr/local/bin/ &&\
    ln -s $REPOS/ercore/erall /usr/local/bin/