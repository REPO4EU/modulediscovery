FROM amazoncorretto:20.0.2-alpine
# use alpine such that ps is available
LABEL image.author.name "Quirin Manz"
LABEL image.author.email "quirin.manz@tum.de"

RUN apk update && \
    apk upgrade && \
    apk add bash

ENV validator_version "6.0.0-SNAPSHOT"

RUN wget -O biopax-validator-${validator_version}-all.zip https://sourceforge.net/projects/biopax/files/validator/biopax-validator-${validator_version}-all.zip/download
RUN unzip biopax-validator-${validator_version}-all.zip
RUN mv distr-${validator_version} /biopax-validator
