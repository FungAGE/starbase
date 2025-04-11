import signal
from contextlib import contextmanager
import logging

logger = logging.getLogger(__name__)


class TimeoutError(Exception):
    pass


@contextmanager
def timeout(seconds):
    """Context manager for adding timeout functionality to operations.

    Args:
        seconds (int): Number of seconds to wait before timing out

    Raises:
        TimeoutError: If the operation takes longer than specified seconds

    Example:
        with timeout(5):
            # Operation that should timeout after 5 seconds
            long_running_operation()
    """

    def timeout_handler(signum, frame):
        raise TimeoutError("Operation timed out")

    # Set the timeout handler
    original_handler = signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(seconds)

    try:
        yield
    finally:
        # Restore the original handler and disable alarm
        signal.alarm(0)
        signal.signal(signal.SIGALRM, original_handler)
