# STARBASE frontend (Dash + thin Flask). Build from repo root:
#   docker build -t starbase-frontend:latest .
#
# All DB/BLAST/HMMER compute is handled by the backend service.
# This image serves the Dash UI and proxies data calls via backend_client.

FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE frontend (Dash UI, no direct DB or BLAST)"

ARG IPSTACK_API_KEY
ARG SECRET_KEY

ENV USER=starbase
ENV HOME=/home/$USER
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV SECRET_KEY=$SECRET_KEY

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/

# curl/wget: miniforge + healthcheck; build-essential: some conda C-ext builds
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl wget build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER $USER

RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase-frontend" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

COPY --chown=$USER:$USER frontend/environment.yaml frontend/environment.yaml
RUN mamba env create -y -f frontend/environment.yaml && \
    mamba clean -afy

ENV PATH=$HOME/miniconda/envs/starbase-frontend/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase-frontend

# Node.js + blasterjs (frontend BLAST visualisation component)
USER root
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Tailscale 
RUN curl -fsSL https://pkgs.tailscale.com/stable/tailscale_1.80.2_amd64.tgz -o /tmp/tailscale.tgz \
    && tar xzf /tmp/tailscale.tgz -C /tmp \
    && mv /tmp/tailscale_1.80.2_amd64/tailscale /tmp/tailscale_1.80.2_amd64/tailscaled /usr/local/bin/ \
    && rm -rf /tmp/tailscale.tgz /tmp/tailscale_1.80.2_amd64
USER $USER

RUN mkdir -p $HOME/src/database/logs \
             $HOME/src/database/cache \
             $HOME/.tailscale && \
    chmod -R 755 $HOME/src/database/logs \
                 $HOME/src/database/cache \
                 $HOME/.tailscale

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

COPY --chown=$USER:$USER ./ ./
RUN chmod +x start-script.sh

USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]
