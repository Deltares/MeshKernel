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
    # logger = log.setup()

    json_reader = JSONReader(arg_parser.file_names())

    attributes = json_reader.attributes()

    families = json_reader.families()

    fig_id = 0
    plotter = Plotter(json_reader, arg_parser.work_dir())

    if json_reader.has_multiple_contenders():
        for family in families:
            plotter.plot(fig_id, family, attributes)
            plotter.save(fig_id, family)
            fig_id += 1

    show_plots()
