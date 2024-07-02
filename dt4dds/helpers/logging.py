import logging
from rich.logging import RichHandler


def default_logging(level='INFO'):
    
    logger = logging.getLogger()
    logger.setLevel(str(level).upper())
    
    handlers = [RichHandler()]
    for handler in handlers:
        handler.setLevel(str(level).upper())

    # create formatter
    formatter = logging.Formatter('%(message)s')

    # add formatter
    for handler in handlers:
        handler.setFormatter(formatter)

    # add ch to logger
    for handler in handlers:
        logger.addHandler(handler)