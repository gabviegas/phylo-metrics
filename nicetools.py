import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# ------------ CARREGAR FASTAS ------------
def load_fasta(filename: str) -> dict[str, str]:
    """Carrega sequências de um arquivo FASTA

    Args:
        filename (str): caminho para o arquivo FASTA a ser lido

    Returns:
        dict[str, str]: dicionário onde as chaves são IDs das
        sequências e os valores são as sequências correspondentes.
    """

    sequences = {}
    current_id = None
    current_seq = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line.split()[0][1:]
                current_seq = []
            elif line and current_id:
                current_seq.append(line.upper())

        if current_id:
            sequences[current_id] = "".join(current_seq)

    return sequences

# ------------ DENDOGRAMAS ------------

def plotar_arvore(distance_matrix, labels, title):
    """
    Plota um dendrograma (árvore) UPGMA a partir de uma matriz de distância.

    Argumentos:
    - distance_matrix (np.array): A matriz de distância N x N (quadrada).
    - labels (list): A lista de N nomes (IDs) das sequências.
    - title (str): O título que você quer dar ao gráfico.
    """

    matriz_total = distance_matrix[np.triu_indices(len(labels), k=1)]
    Z = linkage(matriz_total, method="average")

    sns.set_theme(style="white")
    plt.figure(figsize=(12, 8))

    max_dist = np.max(Z[:, 2])

    dn = dendrogram(
        Z,
        labels=labels,
        orientation="left",
        leaf_font_size=11,
        color_threshold=max_dist * 0.4,
        above_threshold_color="gray",
    )

    plt.title(title, fontsize=16, pad=20)
    plt.xlabel("Distância", fontsize=12)
    plt.ylabel("Sequências", fontsize=12)

    ax = plt.gca()

    ax.xaxis.grid(True, linestyle="--", color="lightgray", zorder=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    plt.tight_layout()
    plt.show()