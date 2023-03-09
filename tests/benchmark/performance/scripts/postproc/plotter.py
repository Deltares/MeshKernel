import logging
import os
import pickle as pkl
from enum import IntEnum, unique

import matplotlib.pyplot as plt
import numpy as np
from json_reader import JSONReader
from matplotlib.ticker import FormatStrFormatter

log = logging.getLogger("Benchmark")

"""
plots, saves and/or displays json results
"""


class Plotter:
    def __init__(self, json_reader, work_dir):
        self.__json_ids = json_reader.ids()
        self.__json_attributes = json_reader.attributes()
        self.__json_families = json_reader.families()
        self.__json_measurements = json_reader.measurements()
        self.__figures = dict()
        self.__figures_dir = self.__create_directory(work_dir, "figures")
        self.__pickles_dir = self.__create_directory(work_dir, "bin")
        self.__figure_id = -1

    @staticmethod
    def __create_directory(path, dir_name):
        dir = os.path.join(path, dir_name)
        if not os.path.isdir(dir):
            log.info("Creating the work directory: %s", dir)
            os.mkdir(dir)
        return dir

    def __exists(self, figure_id):
        """
        checks if fig exists in dict
        """
        return True if figure_id in self.__figures.keys() else False

    @staticmethod
    def __is_open(figure_id):
        """
        checks if the figure with a given id is open
        """
        return plt.fignum_exists(figure_id)

    def __dump(self, figure_id, figure_handle, file_name):
        """
        dumps a figure to a binary file
        """
        self.__figures[figure_id] = dict()
        bin_path = os.path.join(self.__pickles_dir, (file_name + ".bin"))
        self.__figures[figure_id]["binary"] = bin_path
        # open binary stream in write mode and dump the figure handle
        with open(bin_path, "wb") as binary_stream:
            pkl.dump(figure_handle, binary_stream)
        log.info("Dumped %s", bin_path)

    def __load(self, figure_id):
        """
        loads a figure from a binary file
        """
        bin_path = self.__figures[figure_id]["binary"]
        # open binary stream in read mode and load the figure handle
        with open(bin_path, "rb") as binary_stream:
            binary = pkl.load(binary_stream)
        log.info("Loaded %s", bin_path)
        return binary

    @staticmethod
    def __grid_dim(n_figures):
        """
        given a number of sub-figures, returns the nearest integer
        (dimension) required to fit the sub-figures on a square grid
        """
        return int(np.ceil(np.sqrt(n_figures)))

    """
    private enumeration for accessing label size
    """

    @unique
    class __LabelSize(IntEnum):
        AXIS = 0
        TICK = 1
        LEGEND = 2
        TITLE = 3

    __default_format = "png"

    def plot(self, family, attributes, scale="linear", font_size=(10, 8, 10, 16)):
        """
        Plots a figure given a family of benchmarks and set of benchmark attributes.
        Optional:
        - The scale of the y-axis can be specified (see matplotlib.axes.Axes.set_yscale).
          For ex. scale="log".
        - The font size of the axes labels, the legend and tick labels can be set using the parameter
          font_size = (axes_label_font_size tick_labels_font_size, legend_font_size, title_font_size)
        Returns a unique figure id.
        """
        self.__figure_id += 1
        assert self.__exists(self.__figure_id) == False, "Figure ID must be unique"
        n_subplots = len(attributes)
        assert n_subplots > 0
        figure_handle = plt.figure(self.__figure_id)
        grid_dim = self.__grid_dim(n_subplots)
        grid = figure_handle.add_gridspec(grid_dim, grid_dim)
        axes = np.ravel(grid.subplots())

        # delete unused axes
        for axis in axes[n_subplots:]:
            figure_handle.delaxes(axis)
        axes = np.delete(axes, np.s_[n_subplots:])

        # x data
        x_data = [id.pretty_name for id in self.__json_ids]

        for i_attribute, attribute in enumerate(attributes):
            # get axis of subplot
            axis = axes[i_attribute]

            # plot
            for arg in self.__json_families[family]:
                experiment = JSONReader.join_family(family, arg)
                y_data = self.__json_measurements[experiment][attribute]
                axis.plot(x_data, y_data, label=arg, marker="o")

            # set y-axis label abd its font size
            y_label = (
                self.__json_attributes[attribute]["pretty_name"]
                + " ["
                + self.__json_attributes[attribute]["unit"]
                + "]"
            )
            axis.set_ylabel(y_label, fontsize=font_size[self.__LabelSize.AXIS])

            # set y-axis scale
            axis.set_yscale(scale)

            # set the font size of x-axis and y-axis labels
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(font_size[self.__LabelSize.TICK])

            # activate grid
            axis.grid()

        # set title and its font size
        figure_handle.suptitle(family, fontsize=font_size[self.__LabelSize.TITLE])

        # set legend: all figures have the same legends, get those of the last axis
        handles, labels = axis.get_legend_handles_labels()
        figure_handle.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1, 0.75),
            fontsize=font_size[self.__LabelSize.LEGEND],
        )

        # dump to binary file and close
        self.__dump(self.__figure_id, figure_handle, family)

        # all done
        plt.close(figure_handle)

        return self.__figure_id

    def save(self, figure_id, file_name, fmt=__default_format, res="figure"):
        """
        saves a figure given an id, file name, and optionally a format (default is png)
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
            log.info("Saved %s.", path)

    def delete(self, figure_id):
        """
        deletes a figure given its id
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
            log.info("Deleted %s (figure id = %d)", message, figure_id)

    def display(self, figure_id):
        """
        displays a figure given its id
        """
        if self.__exists(figure_id) and not self.__is_open(figure_id):
            log.info("Display %s.", self.__figures[figure_id]["binary"])
            self.__load(figure_id)
            plt.show()

    def close(self, figure_id):
        """
        closes a figure given its id
        """
        log.info("Close %s.", self.__figures[figure_id]["binary"])
        if self.__exists(figure_id) and self.__is_open(figure_id):
            plt.close(figure_id)
