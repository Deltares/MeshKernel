import log
from arg_parser import ArgParser
from json_reader import JSONReader
from plotter import Plotter

import matplotlib.pyplot as plt

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

    plotter = Plotter(json_reader, arg_parser.work_dir())

    if json_reader.has_multiple_contenders():
        for family in families:
            fig_id = plotter.plot(family, attributes)
            plotter.save(fig_id, family)

    # plt.close("all")

    plotter.display(0)
    plotter.close(0)
    plotter.display(1)

    # plotter.delete(0)

    plt.show()
