FROM condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive
# a C++ compiler (with OpenMP) for the SSAF-T extension, git/curl for the CLIs
RUN apt-get update && apt-get install -y --no-install-recommends build-essential git curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Workflow engine, bash and Node.js (for the Codex CLI) in the base env
RUN INSTALL="conda install -y"; \
    command -v mamba >/dev/null 2>&1 && INSTALL="mamba install -y"; \
    command -v micromamba >/dev/null 2>&1 && INSTALL="micromamba install -y -n base"; \
    $INSTALL -c conda-forge -c bioconda "snakemake-minimal>=9,<10" bash "nodejs=22" \
    && conda clean -afy

RUN if [ -x /opt/conda/bin/bash ]; then \
        ln -sf /opt/conda/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/bin/bash /bin/sh; \
    fi

RUN pip install uv

# SSAF-T environment: Python 3.11, the pinned lock file (CUDA 13 stack), then SSAF-T itself with --no-deps; without AVX-512 on the build
# host the C extension is built with -Dsimd=disabled
RUN CREATE="conda create -y"; command -v mamba >/dev/null 2>&1 && CREATE="mamba create -y"; \
    $CREATE -n ssaft -c conda-forge python=3.11 && conda clean -afy
COPY vendor/SSAF-T /opt/SSAF-T
RUN uv pip install --no-cache --python /opt/conda/envs/ssaft/bin/python -r /opt/SSAF-T/requirements-lock.txt \
        --index-strategy unsafe-best-match \
        --extra-index-url https://pypi.nvidia.com \
        --extra-index-url https://download.pytorch.org/whl/cu130
# CeLLMan needs pynndescent < 0.6 (0.6 breaks its descend with scanpy 1.11); SSAF-T's lock has 0.6.0
RUN uv pip install --no-cache --python /opt/conda/envs/ssaft/bin/python "pynndescent==0.5.13"
RUN SIMD=""; grep -q avx512f /proc/cpuinfo || SIMD="--config-settings=setup-args=-Dsimd=disabled"; \
    uv pip install --no-cache --python /opt/conda/envs/ssaft/bin/python --no-deps $SIMD /opt/SSAF-T \
    && /opt/conda/envs/ssaft/bin/python -c "import ssaft, ssaft.utils, spatialdata_io; print('ssaft', ssaft.__version__)"

# CeLLMan (annotation skills, lead and expert agents, helper, cellman-kg MCP server)
COPY vendor/CeLLMan /opt/CeLLMan

# AI CLIs, pinned, self-update off; the login is the server user's, mounted at run time
ENV DISABLE_AUTOUPDATER=1
RUN curl -fsSL https://claude.ai/install.sh | bash -s 2.1.296 \
    && cp -L /root/.local/bin/claude /usr/local/bin/claude && chmod 755 /usr/local/bin/claude \
    && claude --version
RUN npm install -g @openai/codex@0.165.0 && npm cache clean --force && codex --version

# no bytecode caches written next to the scripts (the pipeline folder is mounted from the host)
ENV SSAFT_PYTHON=/opt/conda/envs/ssaft/bin/python \
    CELLMAN_DIR=/opt/CeLLMan \
    NUMPY_MADVISE_HUGEPAGE=0 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /pipeline
# the pipeline folder last: runs mount it over /pipeline, so editing it needs no rebuild
COPY . /pipeline/

CMD ["snakemake", "--help"]
