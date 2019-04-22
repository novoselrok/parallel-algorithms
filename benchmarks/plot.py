import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.style.use('ggplot')

LABELS = {
    'sequential': 'Zaporedno',
    'kmeans': 'k-means',
    'sample-sort': 'Urejanje z vzorčenjem',
    'nbody-bh': 'Barnes-Hut simulacija n teles',
    'c': 'C',
    'julia': 'Julia',
    'chapel': 'Chapel'
}

COLORS = {
    'chapel': '#8dc63f',
    'julia': '#aa7dc0',
    'c': '#1b75b3'
}

languages = ['c', 'chapel', 'julia']

def plot_sequential(seq_df, xlabel, title, filename, use_xaxis_fmt=False):
    for lang in languages:
        lang_df = seq_df[seq_df['language'] == lang]
        sizes = lang_df['size']

        plt.plot(sizes, lang_df['time'], label=LABELS[lang], marker='o', color=COLORS[lang])

    plt.xticks(sizes)
    plt.xlabel(xlabel)
    plt.ylabel('Čas [sekunde]')

    if use_xaxis_fmt:
        plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))

    plt.title(title)
    plt.legend()
    plt.savefig('figures/' + filename + '.pdf', transparent=False, dpi=80)
    plt.clf()
    # plt.show()

data_sizes = {
    'kmeans': [64000, 128000, 256000, 512000],
    'sample-sort': [10000000, 20000000, 40000000, 80000000],
    'nbody-bh': [10000, 20000, 40000, 80000]
}

def plot_workers(workers_df, problem, filename, title, subplot_title_prefix, xlabel):
    fig, axes = plt.subplots(2, 2)
    
    for x in range(2):
        for y in range(2):
            idx = x * 2 + y
            size = data_sizes[problem][idx]
            ax = axes[x, y]

            for lang in languages:
                lang_size_df = workers_df[(workers_df['language'] == lang) & (workers_df['size'] == size)].sort_values(by=['nworkers'])
                nworkers = lang_size_df['nworkers']

                ax.plot(nworkers, lang_size_df['time'], label=LABELS[lang], marker='o', color=COLORS[lang])
                ax.set_xticks(nworkers)
            
            ax.set_ylabel('Čas [sekunde]', fontsize=8)
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_title(subplot_title_prefix + str(size), fontsize=9)
            ax.legend(fontsize=8)
    
    fig.suptitle(title)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    plt.savefig('figures/' + filename + '.pdf', transparent=False, dpi=80)
    plt.clf()
    # fig.show()
    # plt.show()



def main():
    df = pd.read_csv('benchmarks.csv')

    plot_sequential(df[(df['implementation'] == 'sequential') & (df['problem'] == 'kmeans')], 'Število točk', 'k-means - zaporedno', 'kmeans-seq')
    plot_sequential(df[(df['implementation'] == 'sequential') & (df['problem'] == 'sample-sort')], 'Število elementov', 'Urejanje z vzorčenjem - zaporedno', 'sample-sort-seq', use_xaxis_fmt=True)
    plot_sequential(df[(df['implementation'] == 'sequential') & (df['problem'] == 'nbody-bh')], 'Število teles', 'Barnes-Hut simulacija n teles - zaporedno', 'nbody-bh-seq')

    kmeans_sm_df = df[(df['implementation'] == 'sm') & (df['problem'] == 'kmeans')]
    sample_sort_sm_df = df[(df['implementation'] == 'sm') & (df['problem'] == 'sample-sort')]
    nbody_sm_df = df[(df['implementation'] == 'sm') & (df['problem'] == 'nbody-bh')]

    plot_workers(kmeans_sm_df, 'kmeans', 'kmeans-sm', 'k-means - skupni pomnilnik', 'Število točk: ', 'Število niti')
    plot_workers(sample_sort_sm_df, 'sample-sort', 'sample-sort-sm', 'Urejanje z vzorčenjem - skupni pomnilnik', 'Število elementov: ', 'Število niti')
    plot_workers(nbody_sm_df, 'nbody-bh', 'nbody-bh-sm', 'Barnes-Hut simulacija n teles - skupni pomnilnik', 'Število teles: ', 'Število niti')

    kmeans_dm_df = df[(df['implementation'] == 'dm') & (df['problem'] == 'kmeans')]
    sample_sort_dm_df = df[(df['implementation'] == 'dm') & (df['problem'] == 'sample-sort')]
    nbody_dm_df = df[(df['implementation'] == 'dm') & (df['problem'] == 'nbody-bh')]

    plot_workers(kmeans_dm_df, 'kmeans', 'kmeans-dm', 'k-means - porazdeljeni pomnilnik', 'Število točk: ', 'Število procesov')
    plot_workers(sample_sort_dm_df, 'sample-sort', 'sample-sort-dm', 'Urejanje z vzorčenjem - porazdeljeni pomnilnik', 'Število elementov: ', 'Število procesov')
    plot_workers(nbody_dm_df, 'nbody-bh', 'nbody-bh-dm', 'Barnes-Hut simulacija n teles - porazdeljeni pomnilnik', 'Število teles: ', 'Število procesov')


if __name__ == "__main__":
    main()
