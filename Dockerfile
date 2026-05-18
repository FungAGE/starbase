# STARBASE frontend (Dash + thin Flask). Build from repo root:
#   docker build -t starbase-frontend:latest .
#
# All DB/BLAST/HMMER compute is handled by the backend service.
# This image serves the Dash UI and proxies data calls via backend_client.

FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE frontend (Dash UI, no direct DB or BLAST)"

ARG IPSTACK_API_KEY
ARG MAINTENANCE_TOKEN

ENV USER=starbase
ENV HOME=/home/$USER
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV MAINTENANCE_TOKEN=$MAINTENANCE_TOKEN

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/

# curl/wget: miniforge + healthcheck; build-essential: some conda C-ext builds
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl wget build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

COPY environment.yaml .
RUN mamba env create -y -f environment.yaml && \
    mamba clean -afy

ENV PATH=$HOME/miniconda/envs/starbase/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Node.js + blasterjs (frontend BLAST visualisation component)
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

RUN chown -R $USER:$USER $HOME

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

COPY ./ ./
RUN chmod +x start-script.sh && \
    chown -R $USER:$USER $HOME/src

USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]
