import os
from enum import IntEnum, unique

import matplotlib.pyplot as plt
import numpy as np
from json_reader import JSONReader
from matplotlib.ticker import FormatStrFormatter

"""
plots, saves and/or displays json results
"""


class Plotter:
    def __init__(self, json_reader, work_dir):
        self.__work_dir = work_dir
        self.__json_ids = json_reader.ids()
        self.__json_attributes = json_reader.attributes()
        self.__json_families = json_reader.families()
        self.__json_measurements = json_reader.measurements()
        self.__figures = dict()

    def __exists(self, fig_id):
        """
        checks if fig exists in dict
        """
        return True if fig_id in self.__figures.keys() else False

    def __register(self, fig_id, fig, family, attributes):
        """
        registers a figure
        """
        self.__figures[fig_id] = {
            "FIGURE": fig,
            "FAMILY": family,
            "ATTRIBUTES": attributes,
        }

    @staticmethod
    def __grid_dim(n_figs):
        """
        given a number of figures, returns the nearest integer
        (dimension) required to fit the figures on a square grid
        """
        return int(np.ceil(np.sqrt(n_figs)))

    """
    private enumeration for accessing label size
    """

    @unique
    class __LabelSize(IntEnum):
        AXIS = 0
        LEGEND = 1
        TICK = 2

    __default_format = "png"

    def plot(self, fig_id, family, attributes, scale="linear", font_size=(10, 6, 8)):
        """
        Plots a figure given an id, a family of benchmarks and set of benchmark attributes.
        Optional:
        - The scale of the y-axis can be specified (see matplotlib.axes.Axes.set_yscale).
          For ex. scale="log".
        - The font size of the axes labels, the legend and tick labels can be set using the parameter
          font_size = (axes_label_font_size, legend_font_size, tick labels_font_size)
        """
        assert self.__exists(fig_id) == False, "Figure ID must be unique"
        n_subplots = len(attributes)
        assert n_subplots > 0
        fig = plt.figure(fig_id)
        # fig.set_size_inches(12, 7)
        grid_dim = self.__grid_dim(n_subplots)
        grid = fig.add_gridspec(grid_dim, grid_dim)
        axes = np.ravel(grid.subplots())
        # delete unused axes
        for axis in axes[n_subplots:]:
            fig.delaxes(axis)
        axes = np.delete(axes, np.s_[n_subplots:])
        # x data
        x = [id.pretty_name for id in self.__json_ids]

        label_font_size = font_size[self.__LabelSize.AXIS]

        for i_attribute, attribute in enumerate(attributes):
            axis = axes[i_attribute]

            for arg in self.__json_families[family]:
                experiment = JSONReader.join_family(family, arg)
                y = self.__json_measurements[experiment][attribute]
                axis.plot(x, y, label=arg, marker="o")

            axis.set_xlabel("Benchmark", fontsize=label_font_size)
            y_label = (
                self.__json_attributes[attribute]["pretty_name"]
                + " ["
                + self.__json_attributes[attribute]["unit"]
                + "]"
            )
            axis.set_ylabel(y_label, fontsize=label_font_size)
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(font_size[self.__LabelSize.TICK])
            axis.set_yscale(scale)
            axis.grid()
        fig.suptitle(family)
        # all figures have the same legends, get those of the axis
        handles, labels = axis.get_legend_handles_labels()
        # common legend
        fig.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1, 0.75),
        )
        # plt.subplots_adjust(right=0.85)
        # fig.tight_layout()
        self.__register(fig_id, fig, family, attributes)
        plt.close(fig)

    """
    saves a figure given an id, file name, and optionally a format (default is png)
    """

    def save(self, fig_id, file_name, fmt=__default_format, res="figure"):
        if self.__exists(fig_id):
            path = os.path.join(self.__work_dir, file_name + "." + fmt)
            self.__figures[fig_id]["FIGURE"].savefig(
                path, format=fmt, dpi=res, bbox_inches="tight"
            )

    """
    deletes a figure given its id
    """

    def delete(self, fig_id):
        if self.__exists(fig_id):
            if self.__is_open(fig_id):
                self.close(fig_id)
            self.__figures.pop(fig_id)

    @staticmethod
    def __is_open(fig_id):
        """
        checks if the figure with a given id is open
        """
        return plt.fignum_exists(fig_id)

    def display(self, fig_id):
        """
        displays a figure given an id
        create a dummy figure and use its manager to display a figure with id fig_id: expensive!
        """
        if self.__exists(fig_id) and not self.__is_open(fig_id):
            manager = plt.figure(fig_id).canvas.manager
            manager.canvas.figure = self.__figures[fig_id]["FIG"]
            self.__figures[fig_id]["FIG"].set_canvas(manager.canvas)

    def display_and_save(self, fig_id, file_name, fmt=__default_format):
        """
        convenience method that displays a figure and saves it given an id, file name, and optionally a format
        """
        self.display(fig_id)
        self.save(fig_id, file_name, fmt)

    def close(self, fig_id):
        """
        closes a figure given an id
        """
        if self.__exists(fig_id) and self.__is_open(fig_id):
            plt.close(self.__figures[fig_id]["FIGURE"])
