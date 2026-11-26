# SUSS Protein Atlas v5.3 — reproducible core container (Linux x86_64)
# Includes: conda env (Snakemake, Foldseek, FoldMason, MAFFT, FastTree, Rate4Site,
# BLAST+, fpocket, FreeSASA, igraph/Leiden, BioPython) + US-align + P2Rank (Java 17)
# + FoldTree + ESM-Scan.
#
# License-gated or large external resources (InterProScan, DeepTMHMM, Foldseek DBs)
# are mounted at runtime via `-v` and configured in config.yaml.
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        openjdk-17-jre-headless \
        unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/suss-atlas
COPY environment.yml /opt/suss-atlas/environment.yml

RUN micromamba create -y -n suss-atlas -f /opt/suss-atlas/environment.yml \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/suss-atlas/bin:/opt/tools/bin:/opt/p2rank:${PATH}

# 1) US-align (single-file C++ build from pylelab/USalign)
RUN mkdir -p /opt/tools/bin /opt/tools/src \
    && git clone --depth 1 https://github.com/pylelab/USalign.git /opt/tools/src/USalign \
    && g++ -static-libgcc -static-libstdc++ -O3 -ffast-math -lm \
        -o /opt/tools/bin/USalign /opt/tools/src/USalign/USalign.cpp \
    && rm -rf /opt/tools/src/USalign

# 2) P2Rank (Java 17 already installed above)
ARG P2RANK_VERSION=2.4.2
RUN curl -fsSL "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/p2rank_${P2RANK_VERSION}.tar.gz" \
        -o /tmp/p2rank.tar.gz \
    && mkdir -p /opt/p2rank \
    && tar -xzf /tmp/p2rank.tar.gz -C /opt/p2rank --strip-components=1 \
    && chmod +x /opt/p2rank/prank \
    && rm -f /tmp/p2rank.tar.gz

# 3) FoldTree (DessimozLab/fold_tree)
RUN git clone --depth 1 https://github.com/DessimozLab/fold_tree.git /opt/tools/fold_tree

# 4) ESM-Scan (optional GPU/CPU PyTorch mutational tolerance; install torch separately if GPU desired)
RUN git clone --depth 1 https://github.com/xuebingwu/ESM-Scan.git /opt/tools/ESM-Scan

COPY . /opt/suss-atlas
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["snakemake", "--help"]
