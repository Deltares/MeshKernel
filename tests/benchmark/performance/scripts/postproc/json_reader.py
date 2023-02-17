import json
import logging
import sys

log = logging.getLogger("Benchmark")


class JSONReader:
    def __init__(self, file_names):
        self.__file_names = file_names
        self.__data = []
        self.__matches = []
        self.__load_data()
        self.__lookup_matches()

    def __load_data(self):
        """
        Loads json objects into a list of py dictionaries
        """
        for i, file_name in enumerate(self.__file_names):
            log.info("Loading file #%d : %s", i + 1, file_name)
            with open(file_name, "r") as f:
                self.__data.append(json.load(f))
            f.close()

    def __lookup_matches(self):
        """
        Finds benchmarks with matching names in all inpout files
        """
        # create a list of lists of benchamrk names
        names = []
        for data in self.__data:
            node = data["benchmarks"]
            name = []
            for child in node:
                name.append(child["name"])
            names.append(name)

        # filter the matches
        for i in range(1, len(names)):
            self.__matches = set(names[i - 1]).intersection(names[i])

        if self.__matches:
            log.info("Found matches: %s", self.__matches)
        else:
            log.error("Failed to find matches in all parsed json files")
            sys.exit()

    def get_data(self, attributes):
        dico = {}
        for attribute in attributes:
            assert isinstance(attribute, str), "Attribute must be a string"
            for match in self.__matches:
                value = []
                for i, data in enumerate(self.__data):
                    node = data["benchmarks"]
                    for child in node:
                        if match == child["name"]:
                            value.append(child[attribute])
                dico[match, attribute] = value
        return dico

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
