import logging


def get_logger(file_name, level=logging.INFO):
    logging.basicConfig(
        format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
    logger = logging.getLogger(file_name)
    logger.setLevel(level)
    return logger
