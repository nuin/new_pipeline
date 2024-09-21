# lib/db_logger.py

from tinydb import TinyDB, Query
from datetime import datetime
import time
from functools import wraps
from pathlib import Path


def get_sample_db(datadir: Path, sample_id: str) -> TinyDB:
    db_path = datadir / "BAM" / sample_id / f"{sample_id}_log.json"
    return TinyDB(db_path)


def log_to_db(db: TinyDB, action: str, status: str, duration: float, error: str = None):
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "action": action,
        "status": status,
        "duration": duration,
        "error": error
    }
    db.insert(log_entry)


def timer_with_db_log(db: TinyDB):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            try:
                result = func(*args, **kwargs)
                end_time = time.time()
                duration = end_time - start_time
                log_to_db(db, func.__name__, "SUCCESS", duration)
                return result
            except Exception as e:
                end_time = time.time()
                duration = end_time - start_time
                log_to_db(db, func.__name__, "ERROR", duration, str(e))
                raise
        return wrapper
    return decorator