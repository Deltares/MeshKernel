from arg_parser import ArgParser
from html_reporter import HTMLReporter
from json_reader import JSONReader
from logger import Logger
from plotter import Plotter

if __name__ == "__main__":
    # parse command line arguments
    arg_parser = ArgParser()

    # create logger instance
    Logger(
        path=arg_parser.work_dir(), file_name_prefix=arg_parser.log_file_name_prefix()
    )

    log = Logger.get()

    log.info("Arg - Files: {}".format(arg_parser.file_names()))
    log.info("Arg - Work directory: {}".format(arg_parser.work_dir()))
    log.info("Arg - Log file prefix : {}".format(arg_parser.log_file_name_prefix()))

    # read the JSON results
    json_reader = JSONReader(arg_parser.file_names())

    # plot and save the results
    plotter = Plotter(arg_parser.work_dir(), json_reader)
    families = json_reader.families()
    attributes = json_reader.attributes()
    for family in families:
        # plot experiments only when multiple contenders are available
        if json_reader.has_multiple_contenders():
            fig_id = plotter.generate(
                family, attributes, Plotter.XMode.Experiments, ordinate_scale="linear"
            )
            plotter.save(fig_id, family + "_experiments")

        # plot measurements
        fig_id = plotter.generate(
            family, attributes, Plotter.XMode.Measurements, ordinate_scale="log"
        )
        plotter.save(fig_id, family + "_measurements")

    # write html report
    html_reporter = HTMLReporter(json_reader.report(), plotter.report())
    html_reporter.write(arg_parser.work_dir())
