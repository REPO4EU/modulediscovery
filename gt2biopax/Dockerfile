FROM mambaorg/micromamba:1.5.1-alpine
# use alpine such that ps is available
LABEL image.author.name "Quirin Manz"
LABEL image.author.email "quirin.manz@tum.de"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

#COPY --chown=$MAMBA_USER:$MAMBA_USER paxtools-5.3.0.jar /tmp/paxtools-5.3.0.jar
#RUN wget -O paxtools-5.3.0.jar https://sourceforge.net/projects/biopax/files/paxtools/paxtools-5.3.0.jar/download

#update PATH such that singularity exec activate the environment, see: https://github.com/mamba-org/micromamba-docker/issues/155
ENV PATH="${PATH}:${MAMBA_ROOT_PREFIX}/bin"
