import datetime
import logging
import logging.handlers
import os


def setup(name="Benchmark", path=""):
    # format config
    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] (%(module)s) : %(message)s"
    )

    # handler config
    if path:
        date_time = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        file_name = os.path.join(path, "benchmark_" + date_time + ".log")
        handler = logging.FileHandler(file_name)
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    # logger config
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger
