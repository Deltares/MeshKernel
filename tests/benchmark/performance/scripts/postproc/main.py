import log
from arg_parser import ArgParser
from json_reader import JSONReader
from matplotlib.pyplot import show as show_plots
from plotter import Plotter

if __name__ == "__main__":
    arg_parser = ArgParser()

    # to log to a file
    logger = log.setup(
        logger_name="Benchmark",
        file_name_prefix=arg_parser.log_file_name_prefix(),
        path=arg_parser.work_dir(),
        append_timestamp=False,
    )
    # to log in terminal
    ##logger = log.setup()

    logger.info("main message")

    result = JSONReader(arg_parser.file_names())

    # n_files = len(arg_parser.file_names())
    # for i in range(0, n_files):
    #    result.display_node_contents(i, "benchmarks")

    dico = result.get_data(("real_time", "cpu_time"))
    for key in dico:
        print(key)
        for y in dico[key]:
            print(y)
