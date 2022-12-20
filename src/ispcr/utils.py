"""
This module contains various utilities used during in silico PCR.
"""

from typing import Iterator, List, TextIO, Union

from ispcr.FastaSequence import FastaSequence

COLUMN_HEADERS = {
    "fpri": 0,
    "rpri": 1,
    "start": 2,
    "end": 3,
    "length": 4,
    "pname": 5,
    "pseq": 6,
}

BASE_HEADER = (
    "forward_primer\treverse_primer\tstart\tend\tlength\tproduct_name\tproduct_sequence"
)


def desired_product_size(
    potential_product_length: int,
    min_product_length: Union[int, None] = None,
    max_product_length: Union[int, None] = None,
) -> bool:
    """Determines if a potential product's size is in the user's desired product range.

    Inputs
    potential_product_length - int
        The length of the potential product
    min_product_length: int | None
        The minimum product size the user will accept. If None, there is no lower limit.
    max_product_length: int | None
        The maximum product size the user will accept. If None, there is no lower limit.

    Outputs
    A boolean for whether the product length is between the min and max product length.

    Example
    desired_product_size(100, 75, 125)
    True
    desired_product_size(200, 75, 125)
    False
    """

    if min_product_length is None:
        if max_product_length is None:
            return True
        else:
            return potential_product_length <= max_product_length
    else:
        if max_product_length is None:
            return min_product_length <= potential_product_length
        else:
            if max_product_length < min_product_length:
                raise ValueError(
                    "min_product_length cannot be larger than max_product_length"
                )
            return min_product_length <= potential_product_length <= max_product_length


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
    sequences = []
    with open(primer_file) as fin:
        for fasta_sequence in read_fasta(fin):
            sequences.append(fasta_sequence)

    return sequences


def is_valid_cols_string(header_string: str) -> bool:
    """
    Internal helper to check if a header string is valid.
    """
    if not header_string:
        return False

    for col_name in header_string.split():
        if col_name not in COLUMN_HEADERS:
            return False

    return True


def get_column_indices(header_string: str) -> List[int]:
    """
    Returns column indices based on a header string.
    """
    return [COLUMN_HEADERS[col] for col in header_string.split()]


def filter_output_line(output_line: str, column_indices: List[int]) -> str:
    """
    Filters a single line of isPCR results based on selected column indices.
    """

    if not output_line.split():
        return ""
    elif column_indices == list(range(7)):
        return output_line
    else:
        columns = output_line.split()
        return "\t".join([columns[i] for i in column_indices])


def parse_selected_cols(cols: str) -> List[int]:
    if cols != "all":
        if not is_valid_cols_string(cols):
            raise InvalidColumnSelectionError("Invalid header string.")
        else:
            selected_column_indices = get_column_indices(cols)
    else:
        selected_column_indices = list(range(len(BASE_HEADER.split())))
    return selected_column_indices


class InvalidColumnSelectionError(Exception):
    pass
