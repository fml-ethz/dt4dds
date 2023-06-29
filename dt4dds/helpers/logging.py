import logging
import tqdm
import sys

class TqdmLoggingHandler(logging.Handler):

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)  






def default_logging(level='INFO'):
    
    logger = logging.getLogger()
    logger.setLevel(str(level).upper())
    
    handlers = [TqdmLoggingHandler()]
    for handler in handlers:
        handler.setLevel(str(level).upper())

    # create formatter
    formatter = logging.Formatter('[%(asctime)s][%(name)s][%(levelname)s] %(message)s')

    # add formatter
    for handler in handlers:
        handler.setFormatter(formatter)

    # add ch to logger
    for handler in handlers:
        logger.addHandler(handler)