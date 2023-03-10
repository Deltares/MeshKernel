import datetime
import logging
import os


class Singleton(type):
    _instance = None

    def __call__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instance


class Logger(metaclass=Singleton):
    __name = "Benchmark"

    def __init__(
        self,
        path="",
        file_name_prefix="benchmark",
        append_timestamp=False,
    ):
        """
        Logger constructor.
        Sets up a logging instance with name "Benchmark".
        To log in a terminal simply use: logger = log.setup().
        """

        # format configuration
        formatter = logging.Formatter(
            fmt="[%(asctime)s] [%(levelname)s] (%(module)s) : %(message)s"
        )

        # handler configuration
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

        # logger configuration
        logger = logging.getLogger(self.__name)
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)

    @classmethod
    def get(cls):
        """
        Gets the logging instance
        """
        return logging.getLogger(cls.__name)
