# lib/db_logger.py

from pathlib import Path
from tinydb import TinyDB
from datetime import datetime
import functools
import time


def log_to_db(db: TinyDB, message: str, level: str, program: str, sample_id: str, run_id: str):
    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "message": message,
        "level": level,
        "program": program,
        "sample_id": sample_id,
        "run_id": run_id
    }
    db.insert(log_entry)


def timer_with_db_log(db):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time

            # Extract sample_id and run_id from db.name
            db_path = Path(db.name)
            sample_id = db_path.stem.split('_')[0]
            run_id = db_path.parent.parent.parent.name

            log_to_db(db,
                      f"Function {func.__name__} executed in {execution_time:.2f} seconds",
                      "INFO",
                      func.__module__,
                      sample_id,
                      run_id
                      )
            return result

        return wrapper

    return decorator


def get_sample_db(datadir: Path, sample_id: str) -> TinyDB:
    db_path = datadir / "BAM" / sample_id / f"{sample_id}_pipeline_logs.json"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return TinyDB(str(db_path))


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass