import json
import logging

logger = logging.getLogger("Benchmark")


class JSONReader:
    def __init__(self, file_names):
        self.__data = {}
        for i, file_name in enumerate(file_names):
            logger.info("Loading file #%d : %s", i + 1, file_name)
            with open(file_name, "r") as f:
                self.__data[i] = json.loads(f.read())
            f.close()

    """
    Displays the contents of a file given the index
    """

    def display(self, i):
        str = json.dumps(self.__data[i], indent=2)
        print(str)
