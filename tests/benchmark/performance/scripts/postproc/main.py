import log
from arg_parser import ArgParser
from json_reader import JSONReader
from matplotlib.pyplot import show as show_plots
from plotter import Plotter

if __name__ == "__main__":
    arg_parser = ArgParser()

    # to log to a file
    # logger = log.setup("Benchmark", arg_parser.work_dir())
    # to log in terminal
    logger = log.setup()

    logger.info("main message")

    result = JSONReader(arg_parser.file_names())

    n_files = len(arg_parser.file_names())
    for i in range(0, n_files):
        result.display(i)
