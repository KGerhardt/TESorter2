# syntax=docker/dockerfile:1
#
# TEsorter2 with every cascade engine: pyhmmer/HMMER, nail, and BATH.
#
# Build from the repository root:
#   docker build -f container/Dockerfile -t tesorter2:2.0.0 .
#
# BATH is compiled from source here, pinned to v2.0.0-rc7 with the easel commit
# of the bioconda bath recipe (the release used in the TEsorter2 benchmarks). Override with --build-arg to
# change it.

ARG BATH_COMMIT=b65b6bb9cece997a782c2bdd7290be9df14e0d99
ARG EASEL_COMMIT=74e45a88acaa6355eb330799a24e9d6edc798dc4

# ---- stage 1: build BATH ----
FROM debian:bookworm-slim AS bath
ARG BATH_COMMIT
ARG EASEL_COMMIT
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential autoconf git ca-certificates python3 \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /build
RUN git clone https://github.com/TravisWheelerLab/BATH.git \
 && cd BATH && git checkout "${BATH_COMMIT}" \
 && git clone https://github.com/TravisWheelerLab/easel.git \
 && (cd easel && git checkout "${EASEL_COMMIT}") \
 && autoconf && ./configure && make -j"$(nproc)" \
 && mkdir -p /opt/bath/bin \
 && cp src/bathsearch src/bathconvert /opt/bath/bin/ \
 && cp LICENSE /opt/bath/LICENSE \
 && printf 'BATH %s\neasel %s\n' "${BATH_COMMIT}" "${EASEL_COMMIT}" > /opt/bath/VERSION

# ---- stage 2: runtime ----
FROM mambaorg/micromamba:1.5.10-bookworm-slim
LABEL org.opencontainers.image.title="TEsorter2" \
      org.opencontainers.image.version="2.0.0" \
      org.opencontainers.image.source="https://github.com/KGerhardt/TEsorter2" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

COPY --chown=$MAMBA_USER:$MAMBA_USER container/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean -a -y

COPY --from=bath /opt/bath /opt/bath

ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/tesorter2
RUN pip install --no-deps --no-cache-dir /tmp/tesorter2 && rm -rf /tmp/tesorter2

# Pre-convert the bundled protein databases to BATH format beside the
# originals, while site-packages is still writable: the image is read-only
# under Apptainer, and converting at run time costs minutes per database.
COPY --chown=$MAMBA_USER:$MAMBA_USER container/preconvert_bath.py /tmp/preconvert_bath.py
RUN BATH_BIN_DIR=/opt/bath/bin python /tmp/preconvert_bath.py \
 && rm /tmp/preconvert_bath.py

# Plain PATH (not only the micromamba entrypoint) so `apptainer exec` and
# `docker run --entrypoint` see the tools too.
ENV PATH=/opt/conda/bin:/opt/bath/bin:$PATH \
    BATH_BIN_DIR=/opt/bath/bin

WORKDIR /data
CMD ["tesorter2", "--help"]
