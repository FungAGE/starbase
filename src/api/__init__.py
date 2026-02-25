from flask import Flask
from src.api.cache import cache_routes
from src.api.routes import blast_routes, accession_routes, error_handlers
from src.api.synteny_routes import synteny_routes
from src.api.submission_routes import submission_routes
from src.telemetry.routes import telemetry_routes


def register_routes(server: Flask, limiter):
    """Register all blueprints with the Flask server."""
    server.register_blueprint(accession_routes)
    server.register_blueprint(cache_routes)
    server.register_blueprint(blast_routes)
    server.register_blueprint(synteny_routes)
    server.register_blueprint(submission_routes)
    server.register_blueprint(telemetry_routes)
    server.register_blueprint(error_handlers)
