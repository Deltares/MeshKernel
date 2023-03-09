import logging
import os
import pickle as pkl
from enum import IntEnum, unique

import matplotlib.pyplot as plt
import numpy as np
from json_reader import JSONReader

log = logging.getLogger("Benchmark")


class Plotter:
    """
    Plots, saves, displays JSON benchmark results
    """

    def __init__(self, work_dir, json_reader):
        self.__json_ids = json_reader.ids()
        self.__json_num_experiments = json_reader.num_experiments()
        self.__json_attributes = json_reader.attributes()
        self.__json_families = json_reader.families()
        self.__json_measurements = json_reader.measurements()
        self.__figures = dict()
        self.__figures_dir = self.__create_directory(work_dir, "figures")
        self.__pickles_dir = self.__create_directory(work_dir, "bin")
        self.__figure_id = 0

    @staticmethod
    def __create_directory(path, dir_name):
        dir = os.path.join(path, dir_name)
        if not os.path.isdir(dir):
            log.info("Creating the work directory: {}".format(dir))
            os.mkdir(dir)
        return dir

    def __exists(self, figure_id):
        """
        Checks if fig exists in dict
        """
        return True if figure_id in self.__figures.keys() else False

    @staticmethod
    def __is_open(figure_id):
        """
        Checks if the figure with a given id is open
        """
        return plt.fignum_exists(figure_id)

    def __dump(self, figure_id, figure_handle, file_name):
        """
        Dumps a figure to a binary file
        """
        self.__figures[figure_id] = dict()
        bin_path = os.path.join(self.__pickles_dir, (file_name + ".bin"))
        self.__figures[figure_id]["binary"] = bin_path
        # open binary stream in write mode and dump the figure handle
        with open(bin_path, "wb") as binary_stream:
            pkl.dump(figure_handle, binary_stream)
        log.info("Dumped {}".format(bin_path))

    def __load(self, figure_id):
        """
        Loads a figure from a binary file
        """
        bin_path = self.__figures[figure_id]["binary"]
        # open binary stream in read mode and load the figure handle
        with open(bin_path, "rb") as binary_stream:
            binary = pkl.load(binary_stream)
        log.info("Loaded {}".format(bin_path))
        return binary

    @staticmethod
    def __grid_dim(n_figures):
        """
        Given a number of sub-figures, returns the nearest integer
        (dimension) required to fit the sub-figures on a square grid
        """
        return int(np.ceil(np.sqrt(n_figures)))

    """
    private enumeration for accessing label size
    """

    @unique
    class __LabelSizeIndex(IntEnum):
        AXIS = 0
        TICK = 1
        LEGEND = 2
        TITLE = 3

    @unique
    class XMode(IntEnum):
        Experiments = 0
        Measurements = 1

    __file_name_suffix = {
        XMode.Experiments: "_experiments",
        XMode.Measurements: "_measurements",
    }

    def generate(
        self,
        family,
        attributes,
        mode: XMode,
        ordinate_scale="log",
        font_size=(10, 8, 10, 16),
    ) -> int:
        """
        Plots a figure given a family of benchmarks and set of benchmark attributes.
        The mode parameter specifies whether experiments or measurements are shown on the x-axis.
        Optional:
        - The scale of the y-axis can be specified (see matplotlib.axes.Axes.set_yscale).
          For ex. ordinate_scale="log".
        - The font size of the axes labels, the legend and tick labels can be set using the parameter
          font_size = (axes_label_font_size tick_labels_font_size, legend_font_size, title_font_size)
        Returns a unique figure id.
        """
        self.__figure_id += 1
        if self.__exists(self.__figure_id):
            log.debug("Figure ID must be unique")
        n_attributes = len(attributes)
        if not n_attributes > 0:
            log.fatal("At least one attribute must be provided")
        figure_handle = plt.figure(self.__figure_id)
        grid_dim = self.__grid_dim(n_attributes)
        grid = figure_handle.add_gridspec(grid_dim, grid_dim)
        axes = np.ravel(grid.subplots())

        # delete unused axes
        for axis in axes[n_attributes:]:
            figure_handle.delaxes(axis)
        axes = np.delete(axes, np.s_[n_attributes:])

        # x data
        if mode == self.XMode.Experiments:
            x_data = [id.pretty_name for id in self.__json_ids]
        elif mode == self.XMode.Measurements:
            x_data = self.__json_families[family]
        else:
            log.error("Invalid mode")
            return -1

        for i_attribute, attribute in enumerate(attributes):
            # get axis of subplot
            axis = axes[i_attribute]

            # plot
            if mode == 0:
                for arg in self.__json_families[family]:
                    measurement = JSONReader.join_family(family, arg)
                    y_data = self.__json_measurements[measurement][attribute]
                    axis.plot(x_data, y_data, label=arg, marker="o")
            else:
                for experiment in range(self.__json_num_experiments):
                    y_data = list()
                    for arg in self.__json_families[family]:
                        measurement = JSONReader.join_family(family, arg)
                        y_data.append(
                            self.__json_measurements[measurement][attribute][experiment]
                        )
                    axis.plot(
                        x_data,
                        y_data,
                        label=self.__json_ids[experiment].pretty_name,
                        marker="o",
                    )

            # set y-axis label abd its font size
            y_label = (
                self.__json_attributes[attribute]["pretty_name"]
                + " ["
                + self.__json_attributes[attribute]["unit"]
                + "]"
            )
            axis.set_ylabel(y_label, fontsize=font_size[self.__LabelSizeIndex.AXIS])

            # set y-axis scale
            axis.set_yscale(ordinate_scale)

            # set the font size of x-axis and y-axis labels
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(font_size[self.__LabelSizeIndex.TICK])

            # activate grid
            axis.grid()

        # set title and its font size
        figure_handle.suptitle(family, fontsize=font_size[self.__LabelSizeIndex.TITLE])

        # set legend: all figures have the same legends, get those of the last axis
        handles, labels = axis.get_legend_handles_labels()
        figure_handle.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1, 0.75),
            fontsize=font_size[self.__LabelSizeIndex.LEGEND],
        )

        # set the layout
        figure_handle.tight_layout()

        # dump to binary file and close
        file_name = family + self.__file_name_suffix[mode]
        self.__dump(self.__figure_id, figure_handle, file_name)

        # all done
        plt.close(figure_handle)

        return self.__figure_id

    __default_format = "png"

    def save(self, figure_id, file_name, fmt=__default_format, res="figure"):
        """
        Saves a figure given an id, file name, and optionally a format (default is png)
        """
        if self.__exists(figure_id):
            # set path of graphic file
            path = os.path.join(self.__figures_dir, (file_name + "." + fmt))
            # store the path
            self.__figures[figure_id]["graphic"] = path
            # load the associated binary file and save
            figure_handle = self.__load(figure_id)
            figure_handle.savefig(path, format=fmt, dpi=res, bbox_inches="tight")
            plt.close(figure_handle)
            log.info("Saved {}".format(path))
        else:
            log.warning("Request to save figure {} is ignored".format(figure_id))

    def delete(self, figure_id):
        """
        Deletes a figure given its id
        """
        if self.__exists(figure_id):
            if self.__is_open(figure_id):
                self.close(figure_id)
            # remove the binary file
            binary_path = self.__figures[figure_id]["binary"]
            os.remove(binary_path)
            message = binary_path
            # remove graphic file
            if "graphic" in self.__figures[figure_id].keys():
                graphic_path = self.__figures[figure_id]["graphic"]
                os.remove(graphic_path)
                message += " and " + graphic_path
            # remove from the dict
            self.__figures.pop(figure_id)
            log.info("Deleted {} (figure id = {})".format(message, figure_id))
        else:
            log.warning("Request to delete figure {} is ignored".format(figure_id))

    def display(self, figure_id):
        """
        Displays a figure given its id
        """
        if self.__exists(figure_id) and not self.__is_open(figure_id):
            log.info("Display {}".format(self.__figures[figure_id]["binary"]))
            self.__load(figure_id)
            plt.show()
        else:
            log.warning("Request to display figure {} is ignored".format(figure_id))

    def close(self, figure_id):
        """
        Closes a figure given its id
        """
        log.info("Close {}".format(self.__figures[figure_id]["binary"]))
        if self.__exists(figure_id) and self.__is_open(figure_id):
            plt.close(figure_id)
        else:
            log.warning("Request to close figure {} is ignored".format(figure_id))
