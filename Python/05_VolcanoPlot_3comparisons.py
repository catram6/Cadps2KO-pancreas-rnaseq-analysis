import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Volcano plot
# Three comparisons
# ============================================================

plt.rcParams['font.family'] = 'DejaVu Sans'


# ------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------

configs = [
    {
        'file': '6month_deg.xlsx',
        'xlabel': 'log2FC (KO vs WT, 6m)',
        'out': 'volcano_6m_KO_vs_WT.png',
        'fc_thresh': 1,
    },
    {
        'file': '14month_deg.xlsx',
        'xlabel': 'log2FC (KO vs WT, 14m)',
        'out': 'volcano_14m_KO_vs_WT.png',
        'fc_thresh': 1,
    },
    {
        'file': 'WT_deg.xlsx',
        'xlabel': 'log2FC (WT 14m vs 6m)',
        'out': 'volcano_WT_14m_vs_6m.png',
        'fc_thresh': 2,
    },
]


Y_CAP = 60
PADJ_SIG = 0.05


# ------------------------------------------------------------
# Generate volcano plots
# ------------------------------------------------------------

for cfg in configs:

    # Load DESeq2 results
    df = pd.read_excel(
        cfg['file']
    )

    # Remove rows with missing values
    df = df.dropna(
        subset=['log2FoldChange', 'padj']
    ).copy()

    # Remove padj = 0 because log10(0) is undefined
    df = df[df['padj'] > 0]

    # Calculate -log10(adjusted p-value)
    neglog = -np.log10(df['padj'])

    # Limit the maximum y-axis value
    neglog_capped = neglog.clip(
        upper=Y_CAP
    )

    fc = df['log2FoldChange']
    padj = df['padj']
    thresh = cfg['fc_thresh']


    # --------------------------------------------------------
    # Classify genes
    # --------------------------------------------------------

    is_up = (
        (padj < PADJ_SIG) &
        (fc > thresh)
    )

    is_down = (
        (padj < PADJ_SIG) &
        (fc < -thresh)
    )

    is_other = ~(is_up | is_down)


    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 7.5)
    )

    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')


    # Non-significant genes
    ax.scatter(
        fc[is_other],
        neglog_capped[is_other],
        s=14,
        color='#b0b0b0',
        alpha=0.5,
        linewidths=0,
        zorder=2
    )


    # Downregulated genes
    ax.scatter(
        fc[is_down],
        neglog_capped[is_down],
        s=16,
        color='#3b6fd6',
        alpha=0.75,
        linewidths=0,
        zorder=3,
        label=f'DOWN (log2FC<-{thresh}, padj<{PADJ_SIG})'
    )


    # Upregulated genes
    ax.scatter(
        fc[is_up],
        neglog_capped[is_up],
        s=16,
        color='#d63b3b',
        alpha=0.75,
        linewidths=0,
        zorder=3,
        label=f'UP (log2FC>{thresh}, padj<{PADJ_SIG})'
    )


    # --------------------------------------------------------
    # Figure formatting
    # --------------------------------------------------------

    ax.grid(
        True,
        color='#e0e0e0',
        linewidth=1,
        zorder=0
    )

    ax.set_axisbelow(True)


    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color('#888888')


    ax.set_xlabel(
        cfg['xlabel'],
        fontsize=18
    )

    ax.set_ylabel(
        '-log10(padj)',
        fontsize=18
    )

    ax.tick_params(
        axis='both',
        labelsize=14
    )

    ax.set_ylim(
        0,
        Y_CAP + 2
    )


    ax.legend(
        loc='upper left',
        fontsize=11,
        frameon=True
    )


    plt.tight_layout()


    # Save figure
    plt.savefig(
        cfg['out'],
        dpi=150
    )

    plt.close()
