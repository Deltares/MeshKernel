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

    logger.info("main message")

    json_reader = JSONReader(arg_parser.file_names())

    print("Has multiple contenders? ", json_reader.has_multiple_contenders())

    # keys
    print("keys: ", json_reader.keys())
    attributes = json_reader.attributes()
    print("Attributes: ", attributes)
    print(
        "BM_RTree/5000/5000: cpu_time = ",
        list(json_reader.measuremenet("BM_RTree/5000/5000", "cpu_time")),
    )

    # benchmarks = json_reader.measurements()
    # for benchmark in benchmarks:
    #     print(benchmark)
    #     for measurement in benchmarks[benchmark]:
    #         print(measurement, ":", benchmarks[benchmark].get(measurement))
    #         # print(measurement, ":", benchmarks[benchmark][measurement])
    #         # print(measurement, ":", benchmarks[benchmark].items())

    measurements = json_reader.measurements()
    # print(measurements)
    families = json_reader.families()
    print("Families: ", families)
    for family in families:
        print(family)
        for arg in families[family]:
            member = json_reader.join_family(family, arg)
            print(member)
            for attribute in attributes:
                print(attribute, measurements[member][attribute])
