import numpy as np
from scipy.spatial.distance import euclidean
from Levenshtein import distance as levenshtein_distance

def voss_vectors(seqs_dict: dict[str, str]) -> tuple[np.ndarray, list[str]]:
    """Vetoriza uma sequência nucleotídica utilizando o mapeamento VOSS

    Args:
        seqs_dict (dict[str, str]): dicionário em que as chaves são IDs e os
            valores são sequências (strings)

    Returns:
        tuple[np.ndarray, list[str]]: Matriz de distâncias e listas de sequências.
    """
    
    ids = list(seqs_dict.keys())
    
    max_len = 0
    for seq in seqs_dict.values():
        if len(seq) > max_len:
            max_len = len(seq)

    vetores = {}
    for seq_id in ids:
        seq = seqs_dict[seq_id]

        u_A = np.zeros(max_len)
        u_G = np.zeros(max_len)
        u_C = np.zeros(max_len)
        u_T = np.zeros(max_len)

        for i, base in enumerate(seq):
            if i >= max_len:
                break
            if base == "A":
                u_A[i] = 1
            elif base == "G":
                u_G[i] = 1
            elif base == "C":
                u_C[i] = 1
            elif base == "T":
                u_T[i] = 1

        final_vector = np.concatenate((u_A, u_G, u_C, u_T))
        vetores[seq_id] = final_vector
    
    return vetores

# ------------ Matriz Euclidiana ---------------

def euclidean_matrix(seqs_dict: dict[str, str]) -> tuple[np.ndarray, list[str]]:
    """Calcula a matriz Euclidiana usando o mapeamento VOSS

    Args:
        seqs_dict (dict[str, str]): dicionário em que as chaves são IDs e os
            valores são sequências (strings)

    Returns:
        tuple[np.ndarray, list[str]]: _description_
    """
    ids = list(seqs_dict.keys())
    n = len(ids)

    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = seqs_dict[ids[i]]
            seq2 = seqs_dict[ids[j]]

            dist = euclidean(seq1, seq2)
            matrix[i, j] = dist
            matrix[j, i] = dist

    return matrix, ids

# --------------- DISTÂNCIA E MATRIZ HAMMING ------------------

def hamming_distance(seq1: str, seq2: str) -> int:
    """Função auxiliar que contabiliza a distância Hamming para duas strings

    Args:
        seq1 (str): sequência 1 a ser avaliada
        seq2 (str): sequência 2 a ser avaliada

    Returns:
        int: distância hamming
    """

    distance = 0

    for i, j in zip(seq1, seq2):
        if i != j:
            distance += 1

    return distance

def hamming_matrix(seqs_dict: dict[str, str]) -> tuple[np.ndarray, list[str]]:
    """Calcula uma matriz das distâncias hamming para as sequências nucleotídicas informadas.

    Args:
        seqs_dict (dict[str, str]): Dicionário com todas as sequências a serem analisadas

    Raises:
        ValueError: Alguma das sequências não é do mesmo tamanho.

    Returns:
        tuple[np.ndarray, list[str]]: Matriz de distâncias e listas de sequências.
    """

    ids = list(seqs_dict.keys())
    n = len(ids)
    matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = seqs_dict[ids[i]]
            seq2 = seqs_dict[ids[j]]

            try:
                assert len(seq1) == len(seq2)

                dist = hamming_distance(seq1, seq2)

                matrix[i, j] = dist
                matrix[j, i] = dist

            except AssertionError as e:
                raise ValueError("As sequências não tem mesmo tamanho")

    return matrix, ids

# --- Matriz Levenshtein ----

def levenshtein_matrix(seqs_dict: dict[str, str]) -> tuple[np.ndarray, list[str]]:
    """Calcula a matriz de distância de Levenshtein

    Args:
        seqs_dict (dict[str, str]): dicionário em que as chaves são IDs e os
            valores são sequências (strings)

    Returns:
        tuple[np.ndarray, list[str]]: Matriz de distâncias e listas de sequências.
    """

    ids = list(seqs_dict.keys())
    n = len(ids)
    matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = seqs_dict[ids[i]]
            seq2 = seqs_dict[ids[j]]

            dist = levenshtein_distance(seq1, seq2)

            matrix[i, j] = dist
            matrix[j, i] = dist

    return matrix, ids

