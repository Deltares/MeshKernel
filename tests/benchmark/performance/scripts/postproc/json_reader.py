import json
import logging
import sys
from collections import namedtuple

log = logging.getLogger("Benchmark")

FileMeta = namedtuple("FileMeta", ["path", "pretty_name"])

AttributeMeta = namedtuple("AttributeMeta", ["name", "unit"])


class JSONReader:
    def __init__(self, file_names):
        # list  containing data read from input files
        self.__data = list()
        # unique id (NamedPair) assigned to each file (in the order it is parsed)
        self.__ids = list()
        # signals if more than one fole has been parsed
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
            log.info("Loading file %d : %s", i + 1, file_name)
            with open(file_name, "r") as f:
                self.__data.append(json.load(f))
            f.close()
            self.__ids.append(FileMeta(file_name, "benchmark_" + format(i + 1, "02d")))
            log.info(
                "Loaded %s, pretty name = %s",
                self.__ids[-1].path,
                self.__ids[-1].pretty_name,
            )
        log.info("All files loaded")
        if len(self.__data) > 1:
            self.__has_multiple_contenders = True
        else:
            log.info("Only a baseline measurment is available")

    def __lookup_matches(self):
        """
        Finds benchmarks with matching names in all inpout files
        """
        log.info("Looking up matching measurements")
        # create a list of lists of benchamrk names
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
                log.info("Found matches: %s", self.__matches)
            else:
                log.error("Failed to find matches in all parsed json files")
                sys.exit()

    sep = "/"

    def __lookup_families(self):
        log.info("Looking up families of measurements")
        # get unique prefixes (must preserve ordering)
        match_prefixes = list()
        for match in self.__matches:
            prefix = match.split(self.sep, 1)[0]
            if prefix not in match_prefixes:
                match_prefixes.append(prefix)
        # store applicable unique args per prefix
        for match_prefix in match_prefixes:
            match_args = list()
            for match in self.__matches:
                if match.startswith(match_prefix):
                    match_args.append(match.split(self.sep, 1)[1])
            self.__families[match_prefix] = match_args
        log.info("Found families: %s", self.__families)

    def __store_measurements(self):
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
        return prefix + JSONReader.sep + args

    def ids(self):
        return self.__ids

    def has_multiple_contenders(self):
        return self.__has_multiple_contenders

    def keys(self):
        return list(self.__measurements)

    def attributes(self):
        return self.__attributes

    def families(self):
        return self.__families

    def measurements(self):
        return self.__measurements

    def measuremenet(self, key, attribute):
        return self.__measurements[key][attribute]

    def display_contents(self, i):
        """
        Displays the contents of a file given the index
        """
        log.info(
            'Contents of "%s":\n%s',
            self.__file_names[i],
            json.dumps(self.__data[i], indent=2),
        )

    def display_node_contents(self, i, node):
        """
        Displays the contents of a node given the file index and node name
        """
        log.info(
            'Contents of node "%s" in "%s":\n%s',
            node,
            self.__file_names[i],
            json.dumps(self.__data[i][node], indent=2),
        )
