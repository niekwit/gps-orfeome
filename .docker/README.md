# Dockerfiles for each version

This directory contains Dockerfiles associated with each release. 

The Docker image derived from this file contains all Conda environments for each rule, i.e. the whole workflow is run in one image.

These images are shared via [Docker Hub](https://hub.docker.com/repository/docker/niekwit/crispr-screens/general) and are generated as follows (from directory with workflow code):

```shell
$ snakemake --containerize > Dockerfile
$ sudo docker build -t niekwit/gps-orfeome:{VERSION} .
$ sudo docker login
$ sudo docker push niekwit/gps-orfeome:{VERSION}
```