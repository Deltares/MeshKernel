import os

from airium import Airium
from logger import Logger

log = Logger.get()


class HTMLReporter:
    def __init__(self, json_report, plotter_report):
        self.__page = Airium(base_indent="  ")
        self.__page("<!DOCTYPE html>")
        with self.__page.html(lang="en", style="font-family: Arial"):
            with self.__page.head():
                self.__page.meta(charset="utf-8")
                self.__page.title(_t="MeshKernel Benchmark")

            with self.__page.body():
                # summary
                contexts, families = json_report

                self.__page.h1(_t="Summary")

                self.__page.h2(_t="Benchmarks")
                with self.__page.table(border=True, align="center"):
                    with self.__page.tr(klass="header_row", align="left"):
                        self.__page.th(_t="Pretty Name")
                        self.__page.th(_t="Path")
                        self.__page.th(_t="Build Type")
                        self.__page.th(_t="Date")
                        self.__page.th(_t="Host Name")
                        self.__page.th(_t="CPUs")
                        self.__page.th(_t="Freq/CPU (MHz)")

                    for context in contexts:
                        with self.__page.tr():
                            self.__page.td(_t=context.pretty_name)
                            self.__page.td(_t=context.path)
                            self.__page.td(_t=context.library_build_type.capitalize())
                            self.__page.td(_t=context.date)
                            self.__page.td(_t=context.host_name)
                            self.__page.td(_t=context.num_cpus)
                            self.__page.td(_t=context.mhz_per_cpu)

                # itemize the families and measurements
                self.__page.h2(_t="Measurements")
                with self.__page.ol():
                    for family, measurements in families.items():
                        with self.__page.li():
                            self.__page(family)
                            with self.__page.ul():
                                for measurement in measurements:
                                    with self.__page.li():
                                        self.__page(measurement)

                # results
                self.__page.h1(_t="Results")

                for graphic_path in plotter_report:
                    with self.__page.p(align="center"):
                        self.__page.img(
                            src=graphic_path,
                            alt="Image not found!",
                        )

    def write(self, work_dir):
        # write the file
        html_file = os.path.join(work_dir, "report.html")
        with open(html_file, "wb") as f:
            f.write(bytes(self.__page))
        log.info("HTML report was written to {}".format(html_file))

    def display(self):
        print(str(self.__page))
