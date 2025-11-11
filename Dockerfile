# Inherit from base image with common setup
FROM ghcr.io/genecore/starbase:base
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE is a database and toolkit for exploring large transposable elements in Fungi"

ARG IPSTACK_API_KEY
ARG MAINTENANCE_TOKEN

# Set secrets from build args
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV MAINTENANCE_TOKEN=$MAINTENANCE_TOKEN

# Add healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

# Copy application code (changes most frequently, so do this last)
COPY ./ ./
RUN chmod +x start-script.sh && \
    chmod +x start_celery.py && \
    chmod +x manage_celery.py && \
    # Ensure all directories and files are owned by starbase user
    chown -R $USER:$USER $HOME/src

# Switch to user
USER $USER

EXPOSE 8000

ENTRYPOINT ["./start-script.sh"]