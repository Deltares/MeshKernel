import argparse
import os

"""
Command line argument parser
"""


class ArgParser:
    """
    Command line argument parser
    """

    def __init__(self):
        self.__args = self.__parse_args()
        self.__file_names = self.__process_file_names()
        self.__work_directory = self.__process_work_directory()
        self.__log_file_name_prefix = self.__process_log_file_name_prefix()

    def __parse_args(self):
        """
        Creates the parser and adds the required and optional arguments
        """
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
            help="<Required> Work directory where the post-processed JSON results are to be saved",
            type=str,
            required=True,
        )
        parser.add_argument(
            "--log_file_name_prefix",
            "-l",
            help="<Optional> prefix of the log file name which is saved in the work directory",
            type=str,
            required=False,
        )
        return parser.parse_args()

    def __process_file_names(self):
        """
        Processes the required file names argument
        """
        file_names = [
            file_name.strip() for file_name in self.__args.file_names.split(",")
        ]
        for file_name in file_names:
            if not os.path.isfile(file_name):
                raise FileNotFoundError(file_name + "does not exist")
        return file_names

    def __process_work_directory(self):
        """
        Processes the required work directory argument
        """
        work_directory = self.__args.work_directory
        if not os.path.isdir(work_directory):
            print("Creating the work directory: ", work_directory)
            os.mkdir(self.__args.work_directory)
        return work_directory

    def __process_log_file_name_prefix(self):
        """
        Processes the optional log file name prefix argument
        """
        if self.__args.log_file_name_prefix:
            return self.__args.log_file_name_prefix
        else:
            return "benchmark"

    def file_names(self):
        """
        Gets the list of file names
        """
        return self.__file_names

    def work_dir(self):
        """
        Gets the work directory
        """
        return self.__work_directory

    def log_file_name_prefix(self):
        """
        Gets the optional log file name prefix
        """
        return self.__log_file_name_prefix
