import datetime
import logging
import logging.handlers
import os


def setup(
    logger_name="Benchmark",
    file_name_prefix="benchmark",
    path="",
    append_timestamp=False,
):
    # format config
    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] (%(module)s) : %(message)s"
    )

    # handler config
    if path:
        file_name = file_name_prefix
        if append_timestamp:
            date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            file_name += "_" + date_time
        file_name += ".log"
        handler = logging.FileHandler(os.path.join(path, file_name), mode="w")
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    # logger config
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger
