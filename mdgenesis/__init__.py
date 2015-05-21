__all__ = ['Analysis']

import logging
# see the advice on logging and libraries in
# http://docs.python.org/library/logging.html?#configuring-logging-for-a-library
class NullHandler(logging.Handler):
    def emit(self, record):
        pass
h = NullHandler()
logging.getLogger("mdgenesis").addHandler(h)
del h

def start_logging(logfile="mdgenesis.log"):
    """Start logging of messages to file and console.

    The default logfile is named `mdgenesis.log` and messages are
    logged with the tag *mdgenesis*.
    """
    import log
    core.log.create("mdgenesis", logfile=logfile)
    logging.getLogger("mdgenesis").info("mdgenesis STARTED logging to %r", logfile)

def stop_logging():
    """Stop logging to logfile and console."""
    import log
    logging.getLogger("mdgenesis").info("mdgenesis STOPPED logging")
    core.log.clear_handlers(logger)  # this _should_ do the job...

from batchanalysis import BatchAnalysis
