import re
from typing import Union

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import read_sequences_from_file, reverse_complement


def get_pcr_product(
    sequence: FastaSequence,
    forward_primer: FastaSequence,
    reverse_primer: FastaSequence,
    min_product_length: Union[int, None] = None,
) -> str:
    """Returns the products amplified by a pair of primers against a single sequence.

    This function is meant to be used by find_pcr_product, but if you'd like to use this separately,
    you must supply FastaSequence objects. A FastaSequence is just a convenient object to couple a header
    with a DNA string. For example,
        >>> forward_primer = FastaSequence("test_forward_primer", "ACTG")
        >>> reverse_primer = FastaSequence("test_reverse_primer", "ATTA")
        >>> target_sequence = FastaSequence("target_sequence", "ATGCTGATGCATGCTA")

    Inputs
    ------
    sequence: FastaSequence
        The fasta sequence to test for amplification.

    forward_primer: FastaSequence
        The forward primer to use.

    reverse_primer: FastaSequence
        The forward primer to use. Note that you should supply this

    min_product_length: None | int
        If provided, only return those products whose length are greater than or equal to this number.
        Defaults to None, which returns all products found.


    Outputs
    -------
    A tab-separated string containing all of the products amplified by the primers contained in the primer file.
    The fields are:
        1. Forward primer name
        2. Reverse primer name
        3. Start position of the product in the target sequence
        4. End position of the product in the target sequence
        5. Product length
        6. The product

    """
    forward_matches = [
        match.start()
        for match in re.finditer(forward_primer.sequence, sequence.sequence)
    ]
    products = []

    for forward_match in forward_matches:
        tempseq = sequence[forward_match:]
        reverse_matches = [
            match.start()
            for match in re.finditer(
                reverse_complement(reverse_primer.sequence), tempseq
            )
        ]

        for reverse_match in reverse_matches:
            product = tempseq[: reverse_match + len(reverse_primer)]
            start = forward_match
            end = (
                forward_match + reverse_match + len(reverse_primer)
            )  # This is the end in the original sequence
            product_length = len(product)
            if min_product_length is not None and product_length < min_product_length:
                continue

            product_line = f"{forward_primer.header}\t{reverse_primer.header}\t{start}\t{end}\t{product_length}\t{product}"
            products.append(product_line)
    return "\n".join(products)


def find_pcr_product(
    primer_file: str, sequence_file: str, min_product_length: Union[int, None] = None
) -> str:
    """Returns all the products amplified by a set of primers in all sequences in a fasta file.

    Inputs
    ------
    primer_file: str
        The path to the fasta file containing the primers to be tested. Currently, this primer
        file is expected to only contain two sequences, with the forward sequence appearing first.
        For an example:
            >test_1.f
            AGTCA
            >test_2.r
            TTATGC

    sequence_file: str
        The path to the fasta file containing the sequences to test the primers against.

    min_product_length: None | int
        If provided, only return those products whose length are greater than or equal to this number.
        Defaults to None, which returns all products found.


    Outputs
    -------
    A tab-separated string containing all of the products amplified by the primers contained in the primer file.
    The fields are:
        1. Forward primer name
        2. Reverse primer name
        3. Start position of the product in the target sequence
        4. End position of the product in the target sequence
        5. Product length
        6. The product

    """
    primers = read_sequences_from_file(primer_file)
    forward_primer, reverse_primer = primers

    sequences = read_sequences_from_file(sequence_file)
    products = []

    for sequence in sequences:
        new_products = get_pcr_product(
            sequence=sequence,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=min_product_length,
        )
        products.append(new_products)
    return "\n".join(products)
