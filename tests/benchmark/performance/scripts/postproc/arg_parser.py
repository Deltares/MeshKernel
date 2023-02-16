import argparse
import os


class ArgParser:
    def __init__(self):
        self.__args = self.__parse_args()
        self.__file_names = self.__process_file_names()
        self.__work_directory = self.__process_work_directory()

    def __parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--file_names",
            "-f",
            help="<Required> Path to JSON benchmark result file(s)",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--work_directory",
            "-w",
            help="<Required> Directory where the post-processed JSON results are to be saved",
            type=str,
            required=True,
        )
        return parser.parse_args()

    def __process_file_names(self):
        file_names = [
            file_name.strip() for file_name in self.__args.file_names.split(",")
        ]
        for file_name in file_names:
            if not os.path.isfile(file_name):
                raise FileNotFoundError(file_name + "does not exist")
        return file_names

    def __process_work_directory(self):
        if not os.path.isdir(self.__args.work_directory):
            print("Creating the work directory: ", self.__args.work_directory)
            os.mkdir(self.__args.work_directory)
        return self.__args.work_directory

    def file_names(self):
        return self.__file_names

    def work_dir(self):
        return self.__work_directory
