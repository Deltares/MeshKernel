import logger
from arg_parser import ArgParser
from json_reader import JSONReader
from plotter import Plotter

# import matplotlib.pyplot as plt

if __name__ == "__main__":
    # parse command line arguments
    arg_parser = ArgParser()

    # set up the logger now that the work dir is known
    logger = logger.setup(
        logger_name="Benchmark",
        file_name_prefix=arg_parser.log_file_name_prefix(),
        path=arg_parser.work_dir(),
        append_timestamp=False,
    )

    logger.info("Arg - Files: {}".format(arg_parser.file_names()))
    logger.info("Arg - Work directory: {}".format(arg_parser.work_dir()))
    logger.info("Arg - Log file prefix : {}".format(arg_parser.log_file_name_prefix()))

    # read the JSON results
    json_reader = JSONReader(arg_parser.file_names())
    attributes = json_reader.attributes()
    families = json_reader.families()

    # plot and save the results
    plotter = Plotter(arg_parser.work_dir(), json_reader)
    for family in families:
        # plot experiments
        fig_id = plotter.generate(
            family, attributes, Plotter.XMode.Experiments, ordinate_scale="linear"
        )
        plotter.save(fig_id, family + "_experiments")

        # plot measurements
        fig_id = plotter.generate(
            family, attributes, Plotter.XMode.Measurements, ordinate_scale="log"
        )
        plotter.save(fig_id, family + "_measurements")

    # plotter.display(0)
    # plotter.close(0)
    # plotter.display(1)
    # plotter.delete(0)
    # plt.show()

    # json_reader.log_file_content(0)
    # json_reader.log_node_content(0, "benchmarks")
