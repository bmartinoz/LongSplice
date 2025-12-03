# =============================================================================
# STAGE 1 · BUILDER
# =============================================================================
FROM debian:bullseye-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

ARG CHOPPER_VERSION=v0.7.0
ARG RESTRANDER_VERSION=v1.1.1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential cmake git curl wget make zlib1g-dev \
      libbz2-dev liblzma-dev ca-certificates pkg-config && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ---- Rust
ENV RUSTUP_HOME=/opt/rustup CARGO_HOME=/opt/cargo PATH="/opt/cargo/bin:${PATH}"
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --no-modify-path && \
    chmod -R a+w "$RUSTUP_HOME" "$CARGO_HOME"

# ---- Chopper
RUN git clone --depth 1 --branch ${CHOPPER_VERSION} https://github.com/wdecoster/chopper.git /opt/chopper && \
    cd /opt/chopper && cargo build --release && install -m 0755 target/release/chopper /usr/local/bin/chopper

# ---- Restrander
RUN git clone --depth 1 --branch ${RESTRANDER_VERSION} https://github.com/mritchielab/restrander.git /opt/restrander && \
    make -C /opt/restrander -j"$(nproc)" && \
    install -m 0755 /opt/restrander/restrander /usr/local/bin/restrander && \
    mkdir -p /usr/local/share/restrander/config && \
    cp -a /opt/restrander/config/*.json /usr/local/share/restrander/config/

# =============================================================================
# STAGE 2 · RUNTIME
# =============================================================================
FROM mambaorg/micromamba:1.5.8

# Binaries compiled in builder
COPY --from=builder /usr/local/bin/ /usr/local/bin/
# Configs de restrander
COPY --from=builder /usr/local/share/restrander/ /usr/local/share/restrander/

# ---- Create environment from environment.yml
COPY environment.yml /tmp/environment.yml
ARG ENV_REFRESH=0
RUN echo "ENV_REFRESH=${ENV_REFRESH}" && \
    micromamba create -y -n longsplice -f /tmp/environment.yml && \
    micromamba clean -a -y

# Environment PATH and Restrander config path
ENV MAMBA_ROOT_PREFIX=/opt/conda \
    CONDA_PREFIX=/opt/conda/envs/longsplice \
    PATH=/opt/conda/envs/longsplice/bin:/opt/conda/bin:/usr/local/bin:/usr/bin:/bin \
    MAMBA_PROC_LOCK_FALLBACK=1 \
    RESTRANDER_CONFIG_DIR=/usr/local/share/restrander/config

# ---- OS packages + Nextflow
USER root
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# Chromium + system libraries needed for Kaleido/Plotly
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      curl ca-certificates procps git chromium \
      libnss3 libatk-bridge2.0-0 libcups2 libxcomposite1 libxdamage1 \
      libxfixes3 libxrandr2 libgbm1 libxkbcommon0 libpango-1.0-0 libcairo2 \
      libasound2 fonts-liberation && \
    rm -rf /var/lib/apt/lists/*

# Caches/plots and permissions + headless matplotlib backend
ENV BROWSER_PATH=/usr/bin/chromium \
    CHROME_PATH=/usr/bin/chromium \
    CHROMIUM_PATH=/usr/bin/chromium \
    CHOREOGRAPHER_BROWSER_PATH=/usr/bin/chromium \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/.cache/matplotlib \
    TMPDIR=/tmp \
    MPLBACKEND=Agg
RUN mkdir -p /tmp/.cache/mamba/proc /tmp/.cache/matplotlib && chmod -R 1777 /tmp /tmp/.cache

# Nextflow
ENV NXF_VER=25.04.6
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/nextflow && \
    chmod +x /usr/local/bin/nextflow

# FIX NanoPlot/Kaleido/Choreographer (points to chromium in the choreographer package)
RUN set -eux; \
    for SP in /opt/conda/envs/longsplice/lib/python*/site-packages; do \
      mkdir -p "$SP/choreographer/cli/browser_exe"; \
      ln -sf /usr/bin/chromium "$SP/choreographer/cli/browser_exe/chrome"; \
      chmod -R 0777 "$SP/choreographer"; \
    done

# ---- Specific conda packages/fixes
RUN micromamba run -n longsplice micromamba install -y -c bioconda -c conda-forge pybedtools ncls pipettor && \
    micromamba clean -a -y

# ---- (Re)install FLAIR from Git within the same env
ARG FLAIR_GIT_REF=master
ENV FLAIR_GIT_REF=${FLAIR_GIT_REF}
RUN /opt/conda/envs/longsplice/bin/python -m pip install --no-cache-dir --no-deps --upgrade --force-reinstall \
      "flair-brookslab @ git+https://github.com/BrooksLabUCSC/flair.git@${FLAIR_GIT_REF}"

# Wrapper for 'transcriptome' if not exposed in the CLI
RUN bash -lc '\
  set -euo pipefail; \
  FLAIR="/opt/conda/envs/longsplice/bin/flair"; \
  echo "[FLAIR] versión CLI:"; "$FLAIR" --version || true; \
  if ! "$FLAIR" transcriptome -h >/dev/null 2>&1; then \
    echo "[FLAIR] Añadiendo wrapper de transcriptome"; \
    mv "$FLAIR" "${FLAIR}.real"; \
    printf "%s\n" "#!/usr/bin/env bash" \
      "if [ \"\${1-}\" = \"transcriptome\" ]; then" \
      "  shift" \
      "  exec /opt/conda/envs/longsplice/bin/python -m flair.modules.transcriptome \"\$@\"" \
      "else" \
      "  exec \"\$(dirname \"\$0\")/flair.real\" \"\$@\"" \
      "fi" > "$FLAIR"; \
    chmod 0755 "$FLAIR"; \
    printf "%s\n" "#!/usr/bin/env bash" \
      "exec /opt/conda/envs/longsplice/bin/python -m flair.modules.transcriptome \"\$@\"" \
      > "/opt/conda/envs/longsplice/bin/flair-transcriptome"; \
    chmod 0755 "/opt/conda/envs/longsplice/bin/flair-transcriptome"; \
  fi; \
  echo "[FLAIR] Probar ayuda transcriptome:"; \
  "$FLAIR" transcriptome -h | head -n 1 || true \
'

# ---- SUPPA from repo + wrapper
RUN git clone --depth 1 https://github.com/comprna/SUPPA.git /usr/local/share/SUPPA && \
    (head -n1 /usr/local/share/SUPPA/suppa.py | grep -q '^#!' || \
      sed -i '1i #!/usr/bin/env python3' /usr/local/share/SUPPA/suppa.py) && \
    sed -i 's/\r$//' /usr/local/share/SUPPA/suppa.py && \
    chmod 0755 /usr/local/share/SUPPA/suppa.py && \
    ln -sf /usr/local/share/SUPPA/suppa.py /usr/local/bin/suppa.py && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec python3 /usr/local/bin/suppa.py "$@"' > /usr/local/bin/suppa && \
    chmod 0755 /usr/local/bin/suppa

# ---- Ensure R libs in the same env (includes argparse for diffSplice_drimSeq.R)
ENV R_LIBS_USER=/opt/conda/envs/longsplice/lib/R/library
RUN /opt/conda/envs/longsplice/bin/Rscript -e '\
  pkgs <- c("argparse","DESeq2","DRIMSeq","BiocParallel","data.table","optparse","qqman"); \
  suppressPackageStartupMessages({ \
    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); \
  }); \
  if (!requireNamespace("argparse", quietly=TRUE)) install.packages("argparse", repos="https://cloud.r-project.org"); \
  for (p in c("DESeq2","DRIMSeq","BiocParallel")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE); \
  for (p in c("data.table","optparse","qqman")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org"); \
  q(status=0)' || true

# ---- Pipeline scripts
COPY scripts/make_validated_samplesheet.R /usr/local/bin/
COPY scripts/make_contrastsheet.R       /usr/local/bin/
COPY bin/run_bambu.r                    /usr/local/bin/
COPY bin/run_dexseq.r                   /usr/local/bin/
COPY bin/run_deseq2.r                   /usr/local/bin/
COPY bin/suppa_split_file.R             /usr/local/bin/

RUN chmod +x /usr/local/bin/make_validated_samplesheet.R \
             /usr/local/bin/make_contrastsheet.R \
             /usr/local/bin/run_bambu.r \
             /usr/local/bin/run_dexseq.r \
             /usr/local/bin/run_deseq2.r \
             /usr/local/bin/suppa_split_file.R

# ---- Normal user + final verification
USER mambauser
WORKDIR /data
RUN micromamba run -n longsplice bash -lc '\
  echo "--- Verificación final ---"; \
  echo -n "Python:        "; python --version || true; \
  echo -n "Java:          "; java -version 2>&1 | head -n1 || true; \
  echo -n "FLAIR:         "; flair --version || true; \
  echo -n "Minimap2:      "; minimap2 --version || true; \
  echo -n "Chopper:       "; chopper --version || true; \
  echo -n "Restrander:    "; restrander 2>&1 | head -n1 || true; \
  echo -n "Samtools:      "; samtools --version | head -n1 || true; \
  echo -n "Bedtools:      "; bedtools --version | head -n1 || true; \
  echo -n "NanoPlot:      "; NanoPlot --version 2>/dev/null || echo OK; \
  echo -n "MultiQC:       "; multiqc --version | sed "s/multiqc, version //g" || true; \
  echo -n "Nextflow:      "; /usr/local/bin/nextflow -version | head -n1 || true; \
  echo -n "deepTools:     "; bamCoverage --version 2>&1 | awk "{print \$NF}" || true; \
  echo -n "bedToBigBed:   "; bedToBigBed 2>&1 | head -n1 || true; \
  echo -n "R pkgs (quick): "; Rscript -e "cat(all(vapply(c(\"argparse\",\"DESeq2\",\"DRIMSeq\",\"BiocParallel\",\"data.table\",\"optparse\",\"qqman\"), requireNamespace, logical(1), quietly=TRUE)), \"\\n\")"; \
  echo "--- Done ---" \
'
# =============================================================================
