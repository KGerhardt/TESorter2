# TEsorter2 container

Docker and Apptainer images with TEsorter2 and every search engine of its
cascade: pyhmmer/HMMER, nail (+MMseqs2) and BATH. BATH is not on conda, so it
is compiled from source, pinned to the commit used in the TEsorter2
benchmarks (BATH `4d2d388`, easel `365a87a`; see `/opt/bath/VERSION` inside
the image).

## Build

From the repository root:

```bash
docker build -f container/Dockerfile -t tesorter2:2.0.0 .
```

Apptainer, on a machine that has the Docker image:

```bash
apptainer build tesorter2_2.0.0.sif container/tesorter2.def
```

On an HPC system without Docker, build the image elsewhere, push it to a
registry, and pull it with `apptainer pull docker://<registry>/tesorter2:2.0.0`.

## Run

```bash
# Docker: mount the working directory, run as your own user so the
# outputs are not owned by the container's user
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data tesorter2:2.0.0 \
    tesorter2 sequences TEs.fa -p 8

# TEsorter-compatible command line (as called by EDTA)
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data tesorter2:2.0.0 \
    TEsorter TEs.fa -db rexdb-plant -p 8

# Apptainer
apptainer exec tesorter2_2.0.0.sif tesorter2 sequences TEs.fa -p 8
```
