
    cd $REPOS/ercore/docker
    docker build -t metocean/ercore:latest .
    docker pull metocean/ercore:latest


Note:

This image is very big because it inherits everything from ops-base and doesn't need to (only needs git, ssh credentials and cdms2) but i was in a hurry (mpi mpi mpi) 
Feel free to change it please
