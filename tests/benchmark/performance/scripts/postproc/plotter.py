import os
from enum import IntEnum, unique

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

"""
plots, saves and/or dsiplays json results
"""


class Plotter:
    def __init__(self, result, work_dir):
        self.__work_dir = work_dir
        self.__result = result
        self.__figures = dict()

    def __exists(self, fig_id):
        """
        checks if fig exists in dict
        """
        return True if fig_id in self.__figures.keys() else False

    def __register(self, fig_id, fig, x, y_list):
        """
        registers a figure
        """
        self.__figures[fig_id] = {"FIG": fig, "XMETA": x, "YMETA": y_list}

    @staticmethod
    def __grid_dim(n_figs):
        """
        given a number of figures, returns the nearest integer
        (dimension) required to fit the figures on a square grid
        """
        return int(np.ceil(np.sqrt(n_figs)))

    """
    private enumeration for accessing curve meta bu index
    """

    @unique
    class __CurveMeta(IntEnum):
        GROUP = 0
        SET = 1
        SCALE = 2

    """
    private enumeration for accesssing label size
    """

    @unique
    class __LabelSize(IntEnum):
        AXIS = 0
        LEGEND = 1
        TICK = 2

    __defaut_format = "png"

    def plot(self, fig_id, x, y_list, font_size=(10, 6, 8)):
        """
        plots a figure given an id, and x and y meta
        optional parameter: fonst_size = (axes label, legend, ticks)
        """
        assert self.__exists(fig_id) == False, "Figure ID must be unique"
        n_subplots = len(y_list)
        assert len(x) == 3 and n_subplots > 0
        fig = plt.figure(fig_id)
        grid_dim = self.__grid_dim(n_subplots)
        grid = fig.add_gridspec(grid_dim, grid_dim)
        axes = np.ravel(grid.subplots())
        # delete unused axes
        for axis in axes[n_subplots:]:
            fig.delaxes(axis)
        axes = np.delete(axes, np.s_[n_subplots:])
        # x data
        x_attrs, x_data = self.__result.get_set_from_group(
            x[self.__CurveMeta.GROUP], x[self.__CurveMeta.SET]
        )
        for i_y, y in enumerate(y_list):
            assert len(y) == 3
            prop_unit_are_same = True
            last_prop_unit = ()
            axis = axes[i_y]
            for i_set, y_set in enumerate(y[self.__CurveMeta.SET]):
                y_attrs, y_data = self.__result.get_set_from_group(
                    y[self.__CurveMeta.GROUP], y_set
                )
                axis.plot(
                    x_data,
                    y_data,
                    label=y_attrs["Notation"] + " (" + y_attrs["Unit"] + ")",
                )
                prop_unit = (y_attrs["Property"], y_attrs["Unit"])
                if i_set > 0 and last_prop_unit != prop_unit:
                    prop_unit_are_same = False
                else:
                    last_prop_unit = prop_unit
            label_font_size = font_size[self.__LabelSize.AXIS]
            # set axis properties
            axis.set_xlabel(
                x_attrs["Property"] + " (" + x_attrs["Unit"] + ")",
                fontsize=label_font_size,
            )
            if len(y[self.__CurveMeta.SET]) == 1:
                axis.set_ylabel(
                    y_attrs["Name"] + " (" + y_attrs["Unit"] + ")",
                    fontsize=label_font_size,
                )
            elif prop_unit_are_same:
                axis.set_ylabel(
                    last_prop_unit[0] + " (" + last_prop_unit[1] + ")",
                    fontsize=label_font_size,
                )
            axis.legend(fontsize=font_size[self.__LabelSize.LEGEND])
            axis.xaxis.set_major_formatter(FormatStrFormatter("%1.0e"))
            for label in axis.get_xticklabels() + axis.get_yticklabels():
                label.set_fontsize(font_size[self.__LabelSize.TICK])
            axis.set_xscale(x[self.__CurveMeta.SCALE])
            axis.set_yscale(y[self.__CurveMeta.SCALE])
            axis.grid()
        fig.tight_layout()
        self.__register(fig_id, fig, x, y_list)
        plt.close(fig)

    """
    saves a figure given an id, file name, and optionally a format (default is png)
    """

    def save(self, fig_id, file_name, fmt=__defaut_format):
        if self.__exists(fig_id):
            path = os.path.join(self.__work_dir, file_name + "." + fmt)
            self.__figures[fig_id]["FIG"].savefig(path, format=fmt)

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

    def display_and_save(self, fig_id, file_name, fmt=__defaut_format):
        """
        convenience mthod that displays a figure and saves it given an id, file name, and optionally a format
        """
        self.display(fig_id)
        self.save(fig_id, file_name, fmt)

    def close(self, fig_id):
        """
        closes a figure given an id
        """
        if self.__exists(fig_id) and self.__is_open(fig_id):
            plt.close(self.__figures[fig_id]["FIG"])

    def __print_register(self):
        for key, value in self.__figures.items():
            print(
                "ID: {}\nHandle: {}\nx-meta: {}\ny-meta: {}\n".format(
                    key, value["FIG"], value["XMETA"], value["YMETA"], "\n"
                )
            )
