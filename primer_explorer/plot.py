from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

MAX_SECUENCIABLE_LEN = 700


def plot_length_distribs(lengths_by_pair, out_fhand):
    num_rows = len(lengths_by_pair)
    num_cols = 2

    fig = Figure((10, 5 * num_rows))
    FigureCanvas(fig)  # Don't remove it or savefig will fail later

    axes_index = 1
    axess0, axess1 = [], []
    for pair, lengths in lengths_by_pair.items():
        pair = [p.decode() for p in pair]
        if axess0:
            axes = fig.add_subplot(num_rows, num_cols, axes_index, sharex=axess0[0])
        else:
            axes = fig.add_subplot(num_rows, num_cols, axes_index)
        axess0.append(axes)
        axes_index += 1

        title = "PCR length products for pair {}, {} (all)"
        axes.hist(lengths, range=(0, 10000), bins=40)
        axes.set_title(title.format(*pair), fontdict={'fontsize': 10})
        axes.set_ylabel("Number of products")

        title = "PCR length products for pair {}, {} (max {} lenth)"
        lengths = [len_ for len_ in lengths if len_ <= MAX_SECUENCIABLE_LEN]
        if axess1:
            axes2 = fig.add_subplot(num_rows, num_cols, axes_index, sharex=axess1[0])
        else:
            axes2 = fig.add_subplot(num_rows, num_cols, axes_index)

        axess1.append(axes2)
        axes_index += 1
        axes2.hist(lengths, range=(0, MAX_SECUENCIABLE_LEN), bins=28)
        ylims = axes2.get_ylim()

        axes2.vlines(100, 0, ylims[1], colors='r')

        axes2.set_title(title.format(pair[0], pair[1], MAX_SECUENCIABLE_LEN),
                        fontdict={'fontsize': 10})
        axes2.set_ylabel("Number of products")

    fig.tight_layout()
    fig.savefig(out_fhand.name)
