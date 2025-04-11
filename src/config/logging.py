import logging
import warnings

logging.basicConfig(
    level=logging.ERROR,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
    ]
)

warnings.filterwarnings("ignore")

def get_logger(name: str):
    """Get a configured logger for a specific module."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.ERROR)

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger