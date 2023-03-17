import json
import sys
from collections import namedtuple

from logger import Logger

log = Logger.get()

FileMeta = namedtuple("FileMeta", ["path", "pretty_name"])


class JSONReader:
    """
    Reads Google benchmark JSON results
    """

    def __init__(self, file_names):
        """
        Constructs the JSON reader
        """
        # list  containing data read from input files
        self.__data = list()
        # unique id assigned to each file (in the order it is parsed): consists of a path and a pretty name
        self.__ids = list()
        # signals if more than one file has been parsed
        self.__has_multiple_contenders = False
        # list of the matching experiments in all of the input files, equiv to all the experiments if a single file is loaded
        self.__matches = list()
        # dict that maps the experiments to the set of parameters it was performed with
        self.__families = dict()
        # a dictionary of all measurements, maps the name of an experiment to the attributes of interest
        self.__measurements = dict()
        # attributes of interest in each measurement
        self.__attributes = {
            "real_time": {"pretty_name": "Real Time", "unit": "ms"},
            "cpu_time": {"pretty_name": "CPU Time", "unit": "ms"},
            "max_bytes_used": {"pretty_name": "Maximum Bytes Used", "unit": "byte"},
            "total_allocated_bytes": {
                "pretty_name": "Total Allocated Bytes",
                "unit": "byte",
            },
        }

        self.__load_data(file_names)
        self.__lookup_matches()
        self.__lookup_families()
        self.__store_measurements()

    def __load_data(self, file_names):
        """
        Loads json objects into a list of py dictionaries
        """
        for i, file_name in enumerate(file_names):
            log.info("Loading file {} : {}".format(i + 1, file_name))
            with open(file_name, "r") as f:
                self.__data.append(json.load(f))
            self.__ids.append(FileMeta(file_name, "benchmark_" + format(i + 1, "02d")))
            log.info(
                "Loaded {}, pretty name = {}".format(
                    self.__ids[-1].path, self.__ids[-1].pretty_name
                )
            )
        log.info("All files loaded")
        if len(self.__data) > 1:
            self.__has_multiple_contenders = True
        else:
            log.info("Only a baseline measurement is available")

    def __lookup_matches(self):
        """
        Finds benchmarks with matching names in all input files
        """
        log.info("Looking up matching measurements")
        # create a list of lists of benchmark names
        names = list()
        for data in self.__data:
            node = data["benchmarks"]
            name = list()
            for child in node:
                name.append(child["name"])
            names.append(name)

        # filter the matches
        if len(names) == 1:
            self.__matches = names[0]
        else:
            for i in range(1, len(names)):
                self.__matches = [item for item in names[i - 1] if item in names[i]]
            if self.__matches:
                log.info("Found matches: {}".format(self.__matches))
            else:
                log.error("Failed to find matches in all parsed json files")
                sys.exit()

    sep = "/"

    def __lookup_families(self):
        """
        Looks up families of measurements in matching measurements
        """
        log.info("Looking up families of measurements")
        # get unique prefixes (must preserve ordering)
        match_prefixes = set()
        for match in self.__matches:
            match_prefixes.add(match.split(self.sep, 1)[0])
        # store applicable unique args per prefix
        for match_prefix in match_prefixes:
            match_args = list()
            for match in self.__matches:
                if match.startswith(match_prefix):
                    match_args.append(match.split(self.sep, 1)[1])
            self.__families[match_prefix] = match_args
        log.info("Found families: {}".format(self.__families))

    def __store_measurements(self):
        """
        Stores matching measurements
        """
        log.info("Storing matching measurements")
        for match in self.__matches:
            self.__measurements[match] = {}
            for attribute in self.__attributes:
                values = list()
                for i, data in enumerate(self.__data):
                    node = data["benchmarks"]
                    for child in node:
                        if match == child["name"]:
                            values.append(child[attribute])
                            break
                self.__measurements[match][attribute] = values
        log.info("All matching measurements stored")

    @staticmethod
    def join_family(prefix, args):
        """
        Joins a benchmark name with benchmark parameters
        """
        return prefix + JSONReader.sep + args

    def ids(self):
        """
        Returns file ids (pair consisting of path and pretty name)
        """
        return self.__ids

    def num_experiments(self):
        """
        Returns the number of experiments (number of parsed JSON results)
        """
        return len(self.__data)

    def has_multiple_contenders(self):
        """
        Signals if more than one contender are available
        """
        return self.__has_multiple_contenders

    def keys(self):
        """
        Returns keys of the experiments
        """
        return list(self.__measurements)

    def attributes(self):
        """
        Returns the attributes
        """
        return self.__attributes

    def families(self):
        """
        Returns the families
        """
        return self.__families

    def measurements(self):
        """
        Returns the measurements
        """
        return self.__measurements

    def measurement(self, i, key, attribute):
        """
        Returns a measurement
        """
        return self.__measurements[key][attribute][i]

    def log_file_content(self, file_index):
        """
        Displays the contents of a file given the index
        """
        log.info(
            'Contents of file "{}":\n{}'.format(
                self.__ids[file_index].pretty_name,
                json.dumps(self.__data[file_index], indent=2),
            )
        )

    def log_node_content(self, file_index, node):
        """
        Displays the contents of a node given the file index and node name
        """
        log.info(
            'Contents of node "{}" in file "{}":\n{}'.format(
                node,
                self.__ids[file_index].pretty_name,
                json.dumps(self.__data[file_index][node], indent=2),
            )
        )
