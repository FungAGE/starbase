# STARBASE compute backend (FastAPI). Build from repository root:
#   docker build -f backend/Dockerfile -t starbase-backend:latest .
#
# Runtime expects SQLite + BLAST data on disk (see DATA_DIR / src/config/settings.py).
# Optional: BACKEND_API_KEY, CELERY_*, REDIS_* for workers (see backend/docker-compose when added).

FROM python:3.11
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE compute API (FastAPI: SQL, BLAST, HMMER)"

ARG IPSTACK_API_KEY
ARG MAINTENANCE_TOKEN

ENV USER=starbase
ENV HOME=/home/$USER
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV MAINTENANCE_TOKEN=$MAINTENANCE_TOKEN
ENV BACKEND_PORT=8001

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/

# curl/wget: miniforge + healthcheck; build-essential: some pip/conda builds
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl wget build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

COPY backend/environment.yaml backend/environment.yaml
RUN mamba env create -y -f backend/environment.yaml && \
    mamba clean -afy

ENV PATH=$HOME/miniconda/envs/starbase-backend/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase-backend

RUN mkdir -p $HOME/src/database/db && \
    chown -R $USER:$USER $HOME && \
    chmod -R 777 $HOME/src/database/db

COPY ./ ./

RUN chmod +x backend/start-script.sh && \
    chmod +x start_celery.py && \
    chmod +x manage_celery.py && \
    chown -R $USER:$USER $HOME/src $HOME/backend

USER $USER

# Port fixed here; if you change BACKEND_PORT at runtime, override healthcheck in compose/k8s.
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://127.0.0.1:8001/health || exit 1

EXPOSE 8001

ENTRYPOINT ["./backend/start-script.sh"]
