"""
This module contains various utilities used during in silico PCR.
"""

from typing import Iterator, List, TextIO

from ispcr.FastaSequence import FastaSequence


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

    Example
    -------
    >>> reverse_complement('GCTGA')
    'TCAGC'
    """

    complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    try:
        rev_seq = "".join([complements[s] for s in dna_string[::-1]])
    except KeyError as e:
        print(f"Base {e.args} not supported.")
        raise
    return rev_seq


def read_fasta(fasta_file: TextIO) -> Iterator[FastaSequence]:
    """An iterator for fasta files.

    Inputs
    ------
    fasta_file: TextIO
        An open file for reading

    Outputs
    -------
    An iterator yielding the the sequence names and sequences from a fasta file

    Example
    -------
    input_file = 'tests/test_data/sequences/met_r.fa.fasta'
    with open(input_file) as fin:
    for name, seq in read_fasta(fin):
        print(f'{name}\n{seq}')

    """
    name = None
    seq: List[str] = []
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                header = name
                sequence = "".join(seq)
                yield FastaSequence(header, sequence)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name:
        header = name
        sequence = "".join(seq)
        yield FastaSequence(header, sequence)


def read_sequences_from_file(primer_file: str) -> List[FastaSequence]:
    primers = []
    with open(primer_file) as fin:
        for fasta_sequence in read_fasta(fin):
            primers.append(fasta_sequence)

    return primers
