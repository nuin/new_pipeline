# lib/db_logger.py

import functools
import time
from datetime import datetime
from pathlib import Path
from typing import Union
from tinydb import TinyDB


def log_to_db(db: Union[TinyDB, str], message: str, level: str, program: str, sample_id: str, run_id: str):
    """
    :param db: The database instance (TinyDB) or file path (str) where the log entry will be stored.
    :param message: The message to be logged.
    :param level: The severity level of the log message.
    :param program: The name of the program generating the log message.
    :param sample_id: The sample identifier related to the log entry.
    :param run_id: The run identifier for tracking the log entry.
    :return: None
    """
    if isinstance(db, str):
        db = TinyDB(db)

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
    """
    A decorator function that logs the execution time of a given function to a database.

    Args:
        db (TinyDB): The database instance to log to.

    Returns:
        function: A decorator function that can be used to wrap another function.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time

            # Extract sample_id and run_id from db.name
            db_path = Path(db.name)
            sample_id = db_path.stem.split("_")[0]
            run_id = db_path.parent.parent.parent.name

            log_to_db(
                db,
                f"Function {func.__name__} executed in {execution_time:.2f} seconds",
                "INFO",
                func.__module__,
                sample_id,
                run_id,
            )
            return result

        return wrapper

    return decorator


def get_sample_db(datadir: Path, sample_id: str) -> TinyDB:
    """
    Retrieves a TinyDB instance for a specific sample.

    :param datadir: The directory where the data is located.
    :type datadir: Path
    :param sample_id: The identifier of the sample.
    :type sample_id: str
    :return: A TinyDB instance for the specified sample.
    :rtype: TinyDB
    """

    db_path = datadir / "BAM" / sample_id / f"{sample_id}_pipeline_logs.json"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return TinyDB(str(db_path))


if __name__ == "__main__":
    # Add any testing or standalone execution code here
    pass
