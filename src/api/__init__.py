from flask import Flask
from src.api.cache import cache_routes
from src.api.routes import (
    accession_routes,
    backend_routes,
    blast_routes,
    error_handlers,
)
from src.telemetry.routes import telemetry_routes


def register_routes(server: Flask, limiter):
    """Register all blueprints with the Flask server."""
    server.register_blueprint(accession_routes)
    server.register_blueprint(backend_routes)
    server.register_blueprint(cache_routes)
    server.register_blueprint(blast_routes)
    server.register_blueprint(telemetry_routes)
    server.register_blueprint(error_handlers)
