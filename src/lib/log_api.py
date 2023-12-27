import requests
from datetime import datetime

LOG_API = "http://10.106.108.24:8001/logs"



def log_to_api(message, level, program, sample_id, run_id):


  timestamp = datetime.now().isoformat()
  requests.post(LOG_API, json={"timestamp": timestamp,
                               "message": message,
                               "sample_id": sample_id,
                               "run_id": run_id,
                                "level": level,
                                "program": program,
  })
  print("logged to api")

  # "timestamp": "2023-12-27T22:52:13.255Z",
  # "message": "string",
  # "sample_id": "string",
  # "run_id": "string",
  # "level": "string",
  # "program": "string",
  # "dagster_id": "string"
