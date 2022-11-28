"""
This module contains various utilities used during in silico PCR.
"""


def reverse_complement(dna_string: str) -> str:
    """Returns the reverse complement of a DNA string.

    Inputs
    ------
    dna_string: str
        A string representing a DNA sequence. Supported bases are A, C, G, and T.

    Outputs
    -------
    The reverse complement of dna_string.

    Raises
    ------
    KeyError
        Raised if there is a base in dna_string that is not one of ACGT.
    """

    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    try:
        rev_seq = "".join([complements[s] for s in dna_string[::-1]])
    except KeyError as e:
        print(f"Base {e.args} not supported.")
        raise
    return rev_seq
