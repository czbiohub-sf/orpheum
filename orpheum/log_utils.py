import logging
import sys


def get_logger(logger_name, level=logging.ERROR):

    # create logger for prd_ci
    logger = logging.getLogger(logger_name)

    # create formatter and add it to the handlers
    logging.basicConfig(
        stream=sys.stdout,
        format="%(name)s - %(asctime)s %(levelname)s: %(message)s",
        level=level,
    )
    logger.addHandler(logging.StreamHandler(sys.stdout))

    return logger
