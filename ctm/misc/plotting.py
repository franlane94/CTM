import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

from matplotlib import animation
from IPython.display import HTML

from matplotlib import rc
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.transforms import blended_transform_factory

mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
mpl.rcParams["text.usetex"] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cm'
mpl.rcParams["lines.linewidth"] = 2.2
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.labelsize"] = 14.
mpl.rcParams["xtick.top"] = True
mpl.rcParams["xtick.labelsize"] = 14.
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.right"] = True
mpl.rcParams["ytick.labelsize"] = 14.
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.minor.bottom"] = False
mpl.rcParams["xtick.minor.top"] = False
mpl.rcParams["ytick.minor.left"] = False
mpl.rcParams["ytick.minor.right"] = False

colour_1 = "dodgerblue"
colour_2 = "#EE442F"
colour_3 = "forestgreen"
colour_4 = "purple"
colour_5 = "deeppink"

colours = ["black", colour_1, colour_2, colour_3, colour_4, colour_5]
linestyles = ["-", "--", "-.", ":"]

class Plotting(object):

        def single_panel_figure(self, x_data, y_data, line_labels, x_label, y_label, x_min, x_max, y_min, y_max, type, legend=True, loc="upper left", fontsize=12., save=False, filename="power_spec.pdf"):

            fig = plt.figure(figsize=(6, 5))
            gs = gridspec.GridSpec(1, 1)

            gs.update(wspace=0.05, hspace=0.1)

            ax1 = fig.add_subplot(gs[0])

            if type == "loglog":

                for i in range(len(x_data)):

                    ax1.loglog(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "semilogx":

                for i in range(len(x_data)):

                    ax1.semilogx(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "semilogy":

                for i in range(len(x_data)):

                    ax1.semilogy(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "plot":

                for i in range(len(x_data)):

                    ax1.plot(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if legend == True:

                ax1.legend(loc=loc, frameon=False, fontsize=fontsize)

            if save == True:

                plt.savefig(filename, bbox_inches="tight")

        def two_panel_figure(self, x_data_1, y_data_1, x_data_1, y_data_1, line_labels, x_label, y_label, x_min, x_max, y_min, y_max, type, legend=True, loc="upper left", fontsize=12., save=False, filename="power_spec.pdf"):

            fig = plt.figure(figsize=(6, 5))
            gs = gridspec.GridSpec(1, 1)

            gs.update(wspace=0.05, hspace=0.1)

            ax1 = fig.add_subplot(gs[0])

            if type == "loglog":

                for i in range(len(x_data)):

                    ax1.loglog(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "semilogx":

                for i in range(len(x_data)):

                    ax1.semilogx(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "semilogy":

                for i in range(len(x_data)):

                    ax1.semilogy(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if type == "plot":

                for i in range(len(x_data)):

                    ax1.plot(x_data[i], y_data[i], label=line_labels[i], color=colours[i], linestyle=linestyles[i])
                    ax1.set_xlim([x_min, x_max])
                    ax1.set_ylim([y_min, y_max])
                    ax1.set_xlabel(x_label)
                    ax1.set_ylabel(y_label)

            if legend == True:

                ax1.legend(loc=loc, frameon=False, fontsize=fontsize)

            if save == True:

                plt.savefig(filename, bbox_inches="tight")
