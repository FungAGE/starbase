from flask import Flask
from src.api.telemetry import telemetry_routes
from src.api.cache import cache_routes
from src.api.routes import blast_routes, accession_routes, error_handlers

def register_routes(server: Flask, limiter):
    """Register all API routes with the server."""
    server.register_blueprint(telemetry_routes, url_prefix="/api/telemetry")
    server.register_blueprint(cache_routes, url_prefix="/api/cache")
    server.register_blueprint(blast_routes, url_prefix="/api/blast")
    server.register_blueprint(accession_routes, url_prefix="/api/accession")
    server.register_blueprint(error_handlers)