import logging

def create(logger_name="mdgenesis", logfile="mdgenesis.log"):
    logger = logging.getLogger(logger_name)

    logger.setLevel(logging.DEBUG)

    # handler that writes to logfile
    logfile_handler = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile_handler.setFormatter(logfile_formatter)
    logger.addHandler(logfile_handler)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger

def clear_handlers(logger):
    """clean out handlers in the library top level logger

    (only important for reload/debug cycles...)
    """
    for h in logger.handlers:
        logger.removeHandler(h)

class NullHandler(logging.Handler):
    def emit(self, record):
        pass
