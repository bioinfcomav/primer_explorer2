from os.path import splitext
from collections import Counter

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
LABELS = {'title': 'histogram', 'xlabel': 'values',
          'ylabel': 'count', 'minimum': 'minimum',
          'maximum': 'maximum', 'average': 'average',
          'variance': 'variance', 'sum': 'sum',
          'items': 'items', 'quartiles': 'quartiles'}


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
        height = 5.0 * num_rows
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
                 num_bins=None):
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
                               xtickslabel_rotation=xtickslabel_rotation)
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
                   xtickslabel_rotation=None, num_bins=None):
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
                    ylimits=None, ylog_scale=False, num_bins=None):

    plot_hist = HistogramPlotter(counters, xlabel=xlabel, ylabel=ylabel,
                                 xmax=xmax, xmin=xmin, titles=titles,
                                 distrib_labels=distrib_labels, kind=kind,
                                 linestyles=linestyles, ylimits=ylimits,
                                 num_cols=num_cols, ylog_scale=ylog_scale,
                                 plots_per_chart=plots_per_chart,
                                 num_bins=num_bins)
    plot_hist.write_figure(fhand)


class IntCounter(Counter):
    '''This is a subclass  of the Counter on python collections.

    It adds some statistical functionalities to the class'''

    def __init__(self, *args, **kwargs):
        'It inialices a '
        super(IntCounter, self).__init__(*args, **kwargs)
        self.labels = LABELS.copy()

    @property
    def min(self):
        'Get the minimun value'
        return min(self.keys())

    @property
    def max(self):
        'Get the maximun value'
        return max(self.keys())

    @property
    def count(self):
        'It returns the count of the values stored in the array'
        return sum(self.values())

    @property
    def median(self):
        'It calculates the median of the values appended'
        quotient, remainder = divmod(self.count, 2)
        if remainder == 0:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            return (val1 + val2) / 2
        else:
            return self._get_value_for_index(quotient)

    @property
    def sum(self):
        'It gets the sum of the values'
        sum_ = 0
        for index, value in self.items():
            sum_ += (index * value)
        return int(sum_)

    @property
    def average(self):
        'It calculates the average'
        count = self.count
        sum_ = self.sum
        return sum_ / count

    @property
    def variance(self):
        'It gets the variance of the values'
        mean = self.average
        sum_ = 0
        for index, counts in self.items():
            sum_ += ((index - mean) ** 2) * counts
        return sum_ / self.count

    @property
    def quartiles(self):
        'It returns the quartiles'
        num_items = self.count
        if num_items < 4:
            msg = 'At least 4 values are required to calculate the quartiles'
            raise RuntimeError(msg)
        # quartile 1
        quotient, remainder = divmod(num_items + 1, 4)
        if not remainder:
            quartile1 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile1 = (val1 + val2) / 2
        # quartile 3
        quotient, remainder = divmod((num_items + 1) * 3, 4)
        if not remainder:
            quartile3 = self._get_value_for_index(quotient - 1)
        else:
            val1 = self._get_value_for_index(quotient - 1)
            val2 = self._get_value_for_index(quotient)
            quartile3 = (val1 + val2) / 2
        return quartile1, self.median, quartile3

    @property
    def irq(self):
        'It gets the interquartile range'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]
        return quart3 - quart1

    @property
    def outlier_limits(self):
        'It returns the intercuartile'
        # pylint: disable=W0612
        quart1 = self.quartiles[0]
        quart3 = self.quartiles[2]

        iqr = self.irq
        limit_distance = round(iqr * 1.5)

        start = int(quart1 - limit_distance)
        end = int(quart3 + limit_distance)
        return (start, end)

    def _get_value_for_index(self, position):
        '''It takes a position and it returns the value for the given index'''
        cum_count = 0
        for index in sorted(self.keys()):
            count = self[index]
            cum_count += count
            if position <= cum_count - 1:
                return index
        else:
            if position >= cum_count:
                raise IndexError('You asked for an index beyond the scope')
            return index

    def _calculate_dist_range(self, min_, max_, outlier_threshold):
        'it calculates the range for the histogram'
        if ((min_ is not None or max_ is not None) and
                outlier_threshold is not None):
            msg = 'You can not pass max, min and outlier_threslhosld to '
            msg += 'calculate distribution range'
            raise ValueError(msg)

        if min_ is None:
            min_ = self.min
        if max_ is None:
            max_ = self.max
        if outlier_threshold:
            left_limit = self.count * outlier_threshold / 100
            rigth_limit = self.count - left_limit
            left_value = self._get_value_for_index(left_limit)
            rigth_value = self._get_value_for_index(rigth_limit)

            if min_ < left_value:
                min_ = left_value
            if max_ > rigth_value:
                max_ = rigth_value
        return min_, max_

    def calculate_bin_edges(self, min_, max_, n_bins=None):
        'It calculates the bin_edges'
        min_bins = 10
        max_bins = 20
        if n_bins is None:
            num_values = int(max_ - min_)
            if num_values == 0:
                n_bins = 1
            elif num_values < min_bins:
                n_bins = num_values
            else:
                n_bins = int(self.count / 10)
                if n_bins < min_bins:
                    n_bins = min_bins
                if n_bins > max_bins:
                    n_bins = max_bins
                if n_bins > num_values:
                    n_bins = num_values

        # now we can calculate the bin edges
        distrib_span = max_ - min_ if max_ != min_ else 1

        if distrib_span % n_bins:
            distrib_span = distrib_span + n_bins - (distrib_span % n_bins)
        bin_span = distrib_span // n_bins
        bin_edges = [min_ + bin_ * bin_span for bin_ in range(n_bins + 1)]
        return bin_edges

    def calculate_distribution(self, bins=None, min_=None, max_=None,
                               outlier_threshold=None):
        'It returns an histogram with the given range and bin'
        if self.count == 0:
            raise RuntimeError('No items in IntCounter')
        distrib = []
        min_, max_ = self._calculate_dist_range(min_, max_, outlier_threshold)
        if min_ is None or max_ is None:
            return None
        bin_edges = self.calculate_bin_edges(min_, max_, bins)
        for bin_index, left_edge in enumerate(bin_edges):
            try:
                rigth_edge = bin_edges[bin_index + 1]
            except IndexError:
                break
            sum_values = 0

            for index2 in sorted(self.keys()):
                value = self[index2]
                if index2 > rigth_edge:
                    break

                elif (left_edge <= index2 and index2 < rigth_edge or
                      left_edge <= index2 and index2 == max_):
                    sum_values += value

            distrib.append(sum_values)
        return {'counts': distrib, 'bin_limits': bin_edges}

    def update_labels(self, labels):
        'It prepares the labels for output files'
        self.labels.update(labels)

    def count_relative_to_value(self, value, comparison):
        'It counts the ints greater, equal, etc, relative to the given value.'
        return sum([counts for val, counts in self.items()
                    if comparison(val, value)])

    def __add__(self, other):
        'Add counts from two counters.'
        counter_python = super(IntCounter, self).__add__(other)
        return self.__class__(counter_python)

    def __str__(self):
        'It writes some basic stats of the values'
        if self.count != 0:
            labels = self.labels

            # now we write some basic stats
            def format_num(x):
                return '{:,d}'.format(x) if isinstance(x, int) else '%.2f' % x

            text = '{}: {}\n'.format(labels['minimum'], format_num(self.min))
            text += '{}: {}\n'.format(labels['maximum'], format_num(self.max))
            text += '{}: {}\n'.format(labels['average'],
                                      format_num(self.average))

            if labels['variance'] is not None:
                text += '{}: {}\n'.format(labels['variance'],
                                          format_num(self.variance))
            if labels['sum'] is not None:
                text += '{}: {}\n'.format(labels['sum'],
                                          format_num(self.sum))
            if labels['items'] is not None:
                text += '{}: {}\n'.format(labels['items'], self.count)
            text += '\n'
            return text
        return ''
