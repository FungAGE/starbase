import sqlite3
from datetime import datetime
import warnings

warnings.filterwarnings("ignore")

from src.components.sql_engine import telemetry_session_factory, telemetry_connected

# use the telemetry_engine to log requests
def log_request(ip_address, endpoint):
    """Log request details to telemetry database."""
    if not telemetry_connected:
        return

    session = telemetry_session_factory()
    
    try:
        query = """
        INSERT INTO requests (ip_address, endpoint, timestamp)
        VALUES (:ip, :endpoint, :timestamp)
        """
        
        session.execute(
            query,
            {
                "ip": ip_address,
                "endpoint": endpoint, 
                "timestamp": datetime.now()
            }
        )
        session.commit()
        
    except Exception as e:
        session.rollback()
    finally:
        session.close()
