def reverse_complement(nuc_sequence: str) -> str:
    """
    Returns the reverse complement of a nucleotide sequence.
    >>> reverse_complement('ACGT')
    'ACGT'
    >>> reverse_complement('ATCGTGCTGCTGTCGTCAAGAC')
    'GTCTTGACGACAGCAGCACGAT'
    >>> reverse_complement('TGCTAGCATCGAGTCGATCGATATATTTAGCATCAGCATT')
    'AATGCTGATGCTAAATATATCGATCGACTCGATGCTAGCA'
    """
    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rev_seq = "".join([complements[s] for s in nuc_sequence[::-1]])
    return rev_seq
