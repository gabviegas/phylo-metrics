import subprocess

from Bio import SeqIO 
from io import StringIO

def aligning_with_mafft(file_path='sequences.fasta', filename='aligned_seqs.fasta') -> None:
    """Alinhamento de sequências a partir do MAFFT. 

    Args:
        file_path (str, optional): Nome do arquivo de entrada. Defaults to 'sequences.fasta'.
        filename (str, optional): Nome do arquivo de saída. Defaults to 'aligned_seqs.fasta'.
    """

    try:
        with open(file_path, "r") as file:
            content = file.read()
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}") 


    # this mafft code is credited to Sanbomics on YouTube. 
    # link: https://www.youtube.com/watch?v=_cTbADrGLCQ

    print("Alinhando...")

    child = subprocess.Popen(["mafft", "--quiet", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child.stdin.write(content.encode()) # encodando as sequencias 
    child_out = child.communicate()[0].decode("utf8") # decodando as sequencias

    seq_ali = list(SeqIO.parse(StringIO(child_out), "fasta"))
    child.stdin.close()

    print("Fim do alinhamento!")

    # writing a .txt to keep the alignment safe

    with open(filename, "w+") as file:
        for seq in seq_ali:
            file.write(">"+seq.description+"\n")
            file.write(str(seq.seq)+"\n")

    return