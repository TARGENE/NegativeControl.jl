FROM julia:1.9.0-bullseye

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/Amsterdam

ENV JULIA_DEPOT_PATH=/opt

RUN apt-get update && apt-get install -y procps

# Import project, build and precompile

COPY . /NegativeControl 

WORKDIR /NegativeControl 

# Precompile and test project
RUN julia -q --project -e'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile(); Pkg.test()'
