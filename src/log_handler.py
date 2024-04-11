import json
import logging
from datetime import datetime

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


class CustomHttpHandler(logging.Handler):
    def __init__(self, url: str, silent: bool = True):
        """
        Initializes the custom http handler
        Parameters:
            url (str): The URL that the logs will be sent to
            token (str): The Authorization token being used
            silent (bool): If False the http response and logs will be sent
                           to STDOUT for debug
        """
        self.url = url
        self.silent = silent

        # sets up a session with the server
        self.MAX_POOLSIZE = 100
        self.session = session = requests.Session()
        session.headers.update(
            {
                "Content-Type": "application/json",
            }
        )
        self.session.mount(
            "http://",
            HTTPAdapter(
                max_retries=Retry(
                    total=1, backoff_factor=0.5, status_forcelist=[403, 500]
                ),
                pool_connections=self.MAX_POOLSIZE,
                pool_maxsize=self.MAX_POOLSIZE,
            ),
        )

        super().__init__()

    def emit(self, record):
        """
        This function gets called when a log event gets emitted. It recieves a
        record, formats it and sends it to the url
        Parameters:
            record: a log record
        """
        logEntry = self.format(record)
        response = self.session.post(self.url, data=logEntry)

        if not self.silent:
            print(logEntry)
            print(response.content)


# create formatter - this formats the log messages accordingly
formatter = logging.Formatter(
    json.dumps(
        {
            "timestamp": datetime.utcnow().isoformat(),
            "message": "%(message)s",
            "level": "%(levelname)s",
            "program": "%(program)s",
            "sample_id": "%(sample_id)s",
            "run_id": "%(run_id)s",
        }
    )
)


# should be added to script
log = logging.getLogger("")
log.setLevel(logging.INFO)
httpHandler = CustomHttpHandler(url="http://10.106.108.24:8001/logs", silent=False)
httpHandler.setLevel(logging.INFO)
httpHandler.setFormatter(formatter)
log.addHandler(httpHandler)


# log.info("This is a test", extra={'sample_id': sample_id, 'run_id': run_id, 'level': 'INFO', 'program': 'test.py'})
