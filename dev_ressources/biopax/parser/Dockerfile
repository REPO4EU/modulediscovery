FROM mambaorg/micromamba:1.5.1-alpine
# use alpine such that ps is available
LABEL image.author.name "Quirin Manz"
LABEL image.author.email "quirin.manz@tum.de"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="${PATH}:${MAMBA_ROOT_PREFIX}/bin"
