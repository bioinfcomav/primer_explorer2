from os.path import splitext

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

FIGURE_SIZE = (15.0, 11.0)  # inche
COLORMAPS = ['Blues', 'Reds', 'Greens', 'Accent', 'Winter', 'Bone', 'Binary']
COLORS = ['blue', 'red', 'green', 'gray']
MARKERGLHYPS = ['^', 's', 'o', '.', 'D', 'h', 's']
MARKER_SIZE2 = 100
MARKER_SIZE3 = 5.0

BAR = 'bar'
LINE = 'line'


def _guess_output_for_matplotlib(fhand):
    'Given an fhand it guesses if we need png or svg'
    output = None
    if fhand is not None:
        output = splitext(fhand.name)[-1].strip('.')
    if not output:
        output = 'png'
    return output


def get_fig_and_canvas(num_rows=1, num_cols=1, figsize=None):
    if figsize is None:
        height = 5.5 * num_rows
        width = 7.5 * num_cols
        if height > 320.0:
            height = 320.0
        figsize = (width, height)

    fig = Figure(figsize=figsize)
    canvas = FigureCanvas(fig)

    return fig, canvas


class HistogramPlotter(object):

    def __init__(self, counters, plots_per_chart=1, kind=LINE, num_cols=1,
                 xlabel=None, ylabel=None, ylog_scale=False, ylimits=None,
                 distrib_labels=None, titles=None, xmax=None, xmin=None,
                 linestyles=None, figsize=None, xtickslabel_rotation=None,
                 num_bins=None, vlines=None):
        if plots_per_chart > 1 and kind == BAR:
            error_msg = 'if kind is BAR only one plot per chart is allowed'
            raise ValueError(error_msg)
        self.kind = kind
        self.counters = counters
        self.num_cols = num_cols
        self.plots_per_chart = plots_per_chart
        self.ylog_scale = ylog_scale
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.num_bins = num_bins
        self.num_plots, self.num_rows = self._get_plot_dimensions()
        fig, canvas = get_fig_and_canvas(num_rows=self.num_rows,
                                         num_cols=self.num_cols,
                                         figsize=figsize)
        self.figure = fig
        self.canvas = canvas
        axes = self._draw_plot(distrib_labels=distrib_labels, ylimits=ylimits,
                               titles=titles, xmax=xmax, xmin=xmin,
                               linestyles=linestyles, num_bins=num_bins,
                               xtickslabel_rotation=xtickslabel_rotation,
                               vlines=vlines)

        self.axes = axes

    def _get_plot_dimensions(self):
        num_plots, mod = divmod(len(self.counters), self.plots_per_chart)
        if mod != 0:
            num_plots += 1

        num_rows, mod = divmod(num_plots, self.num_cols)
        if mod != 0:
            num_rows += 1
        return num_plots, num_rows

    def write_figure(self, fhand):
        plot_format = _guess_output_for_matplotlib(fhand)
        self.canvas.print_figure(fhand, format=plot_format)
        fhand.flush()

    def _draw_histogram_in_axe(self, counter, axe, xmax=None, xmin=None,
                               title=None, distrib_label=None, linestyle=None,
                               ylimits=None, xtickslabel_rotation=None,
                               num_bins=None):

        try:
            distrib = counter.calculate_distribution(bins=num_bins, max_=xmax,
                                                     min_=xmin)
        except RuntimeError:
            axe.set_title(title + ' (NO DATA)')
            return axe
        except AttributeError as error:
            # if distributions is None
            err_msg = "'NoneType' object has no attribute "
            err_msg += "'calculate_distribution'"
            if err_msg in error:

                axe.set_title(title + ' (NO DATA)')
                return axe
            raise
        if distrib is None:
            axe.set_title(title + ' (NO DATA)')
            return axe

        counts = distrib['counts']
        bin_limits = distrib['bin_limits']
        if self.ylog_scale:
            axe.set_yscale('log')

        if self.xlabel:
            axe.set_xlabel(self.xlabel)
        if self.ylabel:
            axe.set_ylabel(self.ylabel)
        if title:
            axe.set_title(title)

        if self.kind == BAR:
            xvalues = range(len(counts))
            axe.bar(xvalues, counts)

            # the x axis label
            xticks_pos = [value + 0.5 for value in xvalues]

            left_val = None
            right_val = None
            xticks_labels = []
            for value in bin_limits:
                right_val = value
                if left_val:
                    xticks_label = (left_val + right_val) / 2.0
                    if xticks_label >= 10:
                        fmt = '%d'
                    elif xticks_label >= 0.1 and xticks_label < 10:
                        fmt = '%.1f'
                    elif xticks_label < 0.1:
                        fmt = '%.1e'
                    xticks_label = fmt % xticks_label
                    xticks_labels.append(xticks_label)
                left_val = right_val

            # we don't want to clutter the plot
            num_of_xlabels = 15
            step = int(len(counts) / float(num_of_xlabels))
            step = 1 if step == 0 else step
            xticks_pos = xticks_pos[::step]
            xticks_labels = xticks_labels[::step]
            axe.set_xticks(xticks_pos)
            axe.set_xticklabels(xticks_labels, rotation=xtickslabel_rotation)
        elif self.kind == LINE:
            kwargs = {}
            if distrib_label is not None:
                kwargs['label'] = distrib_label
            if linestyle is not None:
                kwargs['linestyle'] = linestyle

            x_values = []
            for index, i in enumerate(bin_limits):
                try:
                    i2 = bin_limits[index + 1]
                except IndexError:
                    break
                x_values.append((i + i2) / 2.0)
            y_values = counts
            axe.plot(x_values, y_values, **kwargs)

        if ylimits is not None:
            axe.set_ylim(ylimits)

        return axe

    def _draw_plot(self, distrib_labels=None, titles=None, xmax=None,
                   xmin=None, linestyles=None, ylimits=None,
                   xtickslabel_rotation=None, num_bins=None, vlines=None):
        counter_index = 0
        axes = []
        for plot_num in range(1, self.num_plots + 1):
            # print num_rows, num_cols, plot_num
            axe = self.figure.add_subplot(self.num_rows, self.num_cols,
                                          plot_num)
            for _ in range(self.plots_per_chart):
                try:
                    counter = self.counters[counter_index]
                    if distrib_labels is None:
                        distrib_label = None
                    else:
                        distrib_label = distrib_labels[counter_index]
                    if linestyles is None:
                        linestyle = None
                    else:
                        linestyle = linestyles[counter_index]
                except IndexError:
                    break
                title = titles[counter_index] if titles else None
                self._draw_histogram_in_axe(counter, axe=axe, xmin=xmin,
                                            xmax=xmax, title=title,
                                            distrib_label=distrib_label,
                                            linestyle=linestyle,
                                            ylimits=ylimits, num_bins=num_bins,
                                            xtickslabel_rotation=xtickslabel_rotation)
                counter_index += 1

            if distrib_labels is not None:
                axe.legend()
            # if vlines:
            #     axe.vlines(vlines, ymin=0, ymax=axe.get_ylim()[1], color='r')

            axes.append(axe)

        return axes


def draw_histogram_in_fhand(counter, fhand, title=None, xlabel=None, xmin=None,
                            xmax=None, ylabel=None, kind=BAR, ylimits=None,
                            ylog_scale=False, figsize=None, num_bins=None,
                            xtickslabel_rotation=None):
    'It draws an histogram and if the fhand is given it saves it'
    plot_hist = HistogramPlotter([counter], xlabel=xlabel, ylabel=ylabel,
                                 xmax=xmax, xmin=xmin, titles=[title],
                                 ylimits=ylimits, kind=kind, figsize=figsize,
                                 xtickslabel_rotation=xtickslabel_rotation,
                                 ylog_scale=ylog_scale, num_bins=num_bins)
    plot_hist.write_figure(fhand)


def draw_histograms(counters, fhand, distrib_labels=None, num_cols=2,
                    plots_per_chart=3, xlabel=None, ylabel=None, titles=None,
                    kind=LINE, xmax=None, xmin=None, linestyles=None,
                    ylimits=None, ylog_scale=False, num_bins=None,
                    xtickslabel_rotation=None, vlines=None):

    plot_hist = HistogramPlotter(counters, xlabel=xlabel, ylabel=ylabel,
                                 xmax=xmax, xmin=xmin, titles=titles,
                                 distrib_labels=distrib_labels, kind=kind,
                                 linestyles=linestyles, ylimits=ylimits,
                                 num_cols=num_cols, ylog_scale=ylog_scale,
                                 plots_per_chart=plots_per_chart,
                                 num_bins=num_bins, vlines=vlines,
                                 xtickslabel_rotation=xtickslabel_rotation)
    plot_hist.write_figure(fhand)
