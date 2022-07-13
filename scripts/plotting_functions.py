import pickle
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

def tuple_to_graph(tup):
    edges = [(tup[i-1], tup[i]) for i in range(1, len(tup))]
    G = nx.DiGraph()
    for e in edges:
        if not G.has_edge(e[0], e[1]):
            G.add_edge(e[0], e[1], weight=1)
        else:
            G[e[0]][e[1]]['weight'] += 1
    return G


def plot_motif_row(plot_data, normalized_y=False,
                   log_y=True, ylabel_xcord=-0.08, keys=['h','rw'],
                   rescale_yticks=False, output_file='', xticklabels=[],
                  plot_motifs=True, title='', mapping=None, plot_hypa=False):
    '''
    plot_data will be a dictionary keyed by motifs with the
    following entries:
        frequency: Observed count of motif in the data.
        avg_hypa/std_hypa: Average and standard deviation of
            motif count when sampled from Hypa ensemble.
        avg_rw/std_rw: Average and standard deviation of
            motif count when sampled from higher-order line
            graph/random walk ensemble.
        hypa_over/hypa_under: Number of edges corresponding to
            the motif that were *-represented based on hypa.
    '''
    # initialization settings
    node_size = 100
    labelsize = 15
    ticklabelsize = 14
    textsize = 16
    if mapping is None:
        # Map motifs to rows in the figure
        mapping = {tup: idx for idx, tup in enumerate(sorted(plot_data.keys()))}
    fig, axs = plt.subplots(nrows=len(mapping), ncols=1, figsize=(8, 8))
    plt.subplots_adjust(wspace=0.25, hspace=0.25)
    gs = matplotlib.gridspec.GridSpec(nrows=len(mapping.keys()), ncols=2,
                                      width_ratios=[.5, 3])

    node_colors = {seq: list(['black'] + ['white']*(len(np.unique(seq))-1))
                   for seq in mapping.keys()}
    foc_m_cols = {'solo': '#777777', 0: '#0B9A85', 1: "#76C1CA", 2:'orange'}
    foc_d_cols = {'solo': '#222222', 0: '#105C51', 1: "#208896", 2:'orange'}
    if plot_hypa:
        x = list(range(0, len(plot_data[next(iter(plot_data))])+4, 2))
    else:
        x = list(range(0, len(plot_data[next(iter(plot_data))]), 2))

    # for each motif
    for row, motif in enumerate(mapping):
        ax = plt.subplot(gs[row, 1])
        bar_heights = [plot_data[motif]['frequency']] +\
                       [plot_data[motif]['avg_' + key]
                        for key in keys]
        num_avg_bars = len(bar_heights)-1
        bar_errs = [0] + [plot_data[motif]['std_' + key] for key in keys]
        if plot_hypa:
            bar_heights.append(plot_data[motif]['hypa_over'] +
                               plot_data[motif]['hypa_under'])
            bar_heights.append(plot_data[motif]['hypa_over'])
            bar_heights.append(plot_data[motif]['hypa_under'])
            bar_errs += [0]*(len(bar_heights)-len(bar_errs))
        nudge = 0.005
        color = [foc_m_cols['solo']] + [foc_m_cols[2]]*2 + [foc_m_cols[0]] + [foc_m_cols[1]]*2
        if plot_hypa:
            color += [foc_m_cols[1]]*3

        ax.bar(x, bar_heights, yerr=bar_errs, width=0.4,
               bottom=[-nudge]*len(bar_heights),
                color=color,
                ec='0.1', lw=1.2)

        if not log_y:
            ax.set_ylim(bottom=0)
            if normalized_y:
                ax.set_ylim(top=1.1)

            if rescale_yticks:
                yticks = ax.get_yticks()
                ax.set_yticks(yticks)
                ax.set_yticklabels([str(round(t/10**3, 1))
                                    if round(int(t/10**3)-t/10**3, 1) != 0
                                    else str(int(t/10**3))
                                    for t in yticks], size=ticklabelsize)
                # ax.set_ylabel(r'Freq. ($10^3$)', size=labelsize)
            else:
                ax.set_ylabel('Frequency', size=labelsize)
        else:
            ax.set_yscale('symlog')
            ax.set_ylabel('Frequency (log)', size=labelsize)
            if normalized_y:
                ax.set_ylim(top=1.01)
                ax.set_yticks([0.01, 0.1, 1.0])
                ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        ax.set_axisbelow(True)
        ax.set_xticks([])
        ax.grid(which='major', axis='y')
        ax.set_xlim(-1.0, max(x)+1)
        ax.yaxis.set_label_coords(ylabel_xcord, 0.5)
        ymax = ax.get_ylim()[1]
        ax.fill_between([-1, 1], 0, ymax, alpha=0.2,
                        color=foc_m_cols['solo'], lw=1)
        ax.fill_between([1, 5], 0, ymax, alpha=0.2,
                        color=foc_m_cols[2], lw=1)
        ax.fill_between([5, 7], 0, ymax, alpha=0.2,
                        color=foc_m_cols[0], lw=1)
        ax.fill_between([7, 11], 0, ymax, alpha=0.2,
                        color=foc_m_cols[1], lw=1)


        if plot_hypa:
            ax.fill_between([num_avg_bars*2+1, num_avg_bars*2+8], 0, ymax, alpha=0.2,
                        color=foc_m_cols[1], lw=1)
        if row == 0:
            if title != '':
                ax.text(5.0, ymax*1.05, title, ha='center', va='bottom',
                        fontsize=textsize+2, weight='bold')
            if plot_hypa:
                ax.text(3.8, ymax*1.02, 'Ensemble View', ha='center', va='bottom',
                    color=foc_d_cols[0], fontsize=textsize)
                ax.text(3.5+num_avg_bars*2, ymax*1.02, 'HYPA View', ha='center', va='bottom',
                        color=foc_d_cols[1], fontsize=textsize)
        if row == 4:
            if plot_motifs:
                ax.text(-2.7, ymax / 2, fr"Frequency ($10^3$)", ha='center', va='bottom',
                    fontsize=textsize, rotation=90)

            else:
                ax.text(-2.4, ymax / 2, fr"Frequency ($10^3$)", ha='center', va='bottom',
                    fontsize=textsize, rotation=90)


    ax.set_xticks(x)
    if len(xticklabels) == 0:
        ax.set_xticklabels(['Total', 'Hypa', r'$G^k_c\nUnweighted$', 'Significant', 'Over', 'Under'],
                               size=textsize)
    else:
        ax.set_xticklabels(xticklabels, size=textsize)

    if plot_motifs:
        # Y axis: Plot motifs
        for q, seq in enumerate(mapping.keys()):
            ax = plt.subplot(gs[q, 0])
            g = tuple_to_graph(seq)
            pos = nx.circular_layout(g)
            xpos_q = list(list(zip(*list(pos.values())))[0])
            ypos_q = list(list(zip(*list(pos.values())))[1])
            xdiff = max(xpos_q) - min(xpos_q)
            ydiff = max(ypos_q) - min(ypos_q)
            if seq in [('A', 'B', 'C', 'D'), ('A','B','C')]:
                pos = {k: j*1.15 for k, j in pos.items()}

            nx.draw_networkx_nodes(g, pos, node_size=node_size, node_color=node_colors[seq],
                                   linewidths=1.75, edgecolors='black',ax=ax)
            for u, v, edat in g.edges(data=True):
                for enum in range(edat['weight']):
                    ax.annotate("",
                        xy=pos[u], xycoords='data',
                        xytext=pos[v], textcoords='data',
                        arrowprops=dict(arrowstyle="<|-", color='black',
                                        alpha=0.9,
                                        shrinkA=6, shrinkB=6,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr',str(0.2+0.2*enum)),
                                       ),
                        )
            ax.set_xlim(min(xpos_q)-xdiff*1.5, max(xpos_q)+xdiff*1.5)
            ax.set_ylim(min(ypos_q)-ydiff*0.3, max(ypos_q)+ydiff*0.3)
            ax.set_axis_off()

    plt.tight_layout()

    if output_file != '':
        plt.savefig(output_file)


def read_julia_output(filename, keys=['lg-s'], kvals=[2,3]):
    plot_data = dict()
    for k in kvals:
        for ensemble in keys:
            with open(filename.format(k, ensemble), 'r') as fin:
                for line in fin:
                    arr = line.strip().split('|')
                    # Motif
                    m = tuple([u.upper() for u in arr[0].split(',')])
                    # Frequency
                    freq = float(arr[1])
                    # Samples
                    sampled_counts = list(map(float, arr[2].split(',')))
                    plot_data.setdefault(m, dict())
                    if sum(sampled_counts) == 0:
                        print(f"Warning: no observations of motif {m} in ensemble {ensemble}.")
                        continue
                    plot_data[m].update({'avg_'+ensemble: np.mean(sampled_counts),
                                          'std_'+ensemble: np.std(sampled_counts),
                                          'frequency': freq})

    return plot_data

def get_plot_data(filename, normalized=False, keys=['h', 'rw'], kvals=[2,3]):
    if not normalized:
        plot_data = dict()
        for k in kvals:
            for ensemble in keys:
                with open(filename.format(k, ensemble), 'rb') as fpickle:
                    sampled_counts, data = pickle.load(fpickle)
                for m, d in data.items():
                    plot_data.setdefault(m, dict())
                    if len(sampled_counts[m]) == 0:
                        print(f"Warning: no observations of motif {m} inensemble {ensemble}.")
                        continue
                    if 'hypa_over' in d:
                        plot_data[m].update({'avg_'+ensemble: np.mean(sampled_counts[m]),
                                                 'std_'+ensemble: np.std(sampled_counts[m]),
                                                 'frequency':d['frequency'],
                                                 'hypa_over':d['hypa_over'],
                                                 'hypa_under':d['hypa_under']})
                    else:
                        plot_data[m].update({'avg_'+ensemble: np.mean(sampled_counts[m]),
                                                 'std_'+ensemble: np.std(sampled_counts[m]),
                                                 'frequency':d['frequency'],})

        return plot_data
    else:
        normalized_plot_data = dict()
        for k in kvals:
            for ensemble in keys:
                with open(filename.format(k, ensemble), 'rb') as fpickle:
                    sampled_counts, data = pickle.load(fpickle)
                arr = []
                normed_arrs = []
                motifs = list(sampled_counts.keys())
                total_freq = sum([d['frequency'] for _, d in data.items()])
                for i, m in enumerate(motifs):
                    d = sampled_counts[m]
                    if i > 0:
                        if len(m) > len(motifs[i-1]):
                            arr = np.array(arr, dtype='float64')
                            arr /= arr.sum(axis=0)
                            normed_arrs.append(arr)
                            arr = []
                    arr.append(d)
                arr = np.array(arr, dtype='float64')
                arr /= arr.sum(axis=0)
                normed_arrs.append(arr)
                norm_arrs = np.array(arr)
                for i, m in enumerate(sampled_counts.keys()):
                    normalized_plot_data.setdefault(m, dict())
                    d = data[m]
                    if 'hypa_over' in d:
                        normalized_plot_data[m].update({'avg_'+ensemble: norm_arrs[i].mean(),
                                                 'std_'+ensemble: norm_arrs[i].std(),
                                                 'frequency':d['frequency'] / total_freq,
                                                 'hypa_over':d['hypa_over'] / total_freq,
                                                 'hypa_under':d['hypa_under'] / total_freq})
                    else:
                        normalized_plot_data[m].update({'avg_'+ensemble: norm_arrs[i].mean(),
                                                 'std_'+ensemble: norm_arrs[i].std(),
                                                 'frequency':d['frequency']})

                keysums = dict()
                for m in data.keys():
                    for key in data[m]:
                        if 'avg' not in key and 'std' not in key:
                            keysums.setdefault(key, 0)
                            keysums[key] += data[m][key]

                for m in data.keys():
                    for key in data[m]:
                        if 'avg' not in key and 'std' not in key:
                            normalized_plot_data[m][key] = data[m][key] / keysums[key]
        return normalized_plot_data
