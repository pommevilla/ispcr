import re
from typing import Union

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import (
    desired_product_size,
    filter_output_line,
    parse_selected_cols,
    read_sequences_from_file,
    reverse_complement,
)

BASE_HEADER = (
    "forward_primer\treverse_primer\tstart\tend\tlength\tproduct_name\tproduct_sequence"
)


def calculate_pcr_product(
    sequence: FastaSequence,
    forward_primer: FastaSequence,
    reverse_primer: FastaSequence,
    min_product_length: Union[int, None] = None,
    max_product_length: Union[int, None] = None,
    header: bool = True,
    cols: str = "all",
    output_file: Union[bool, str] = False,
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

    max_product_length: None | int
        If provided, only return those products whose length are less than or equal to this number.
        Defaults to None, which returns all products found.

    header: bool | str
        Whether or not to print a header on the results. Defaults to True. False will not print out the header.

    cols: str
        Which columns to print out. Defaults to "all," which prints out all the columns. A string can be supplied to
        only output the strings of interest. For example, cols="fpri rpri pname" will only output the names of the forward
        primer, reverse primer, and the target sequence when a target is found.
        Available options are:
            fpri - the name of the forward primer
            rpri - the name of the reverse primer
            start - the start location of the product in the target sequence
            end - the end location of the product in the target sequence
            length - the length of the product
            pname - the name of the sequence in which the target was found
            pseq - the nucleotide sequnce of the amplified product

    output_file: bool | str
        The file to write the results out to. Defaults to False, which will not print anything out. Providing a string
        will create that input file at that location. If set to True without providing a string, the output file
        will be of the form <DD-MM-YYYY_HH:MM:SS>.txt



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

    products = []

    # Check cols string
    selected_column_indices = parse_selected_cols(cols)

    if header is True:
        products.append(filter_output_line(BASE_HEADER, selected_column_indices))

    forward_matches = [
        match.start()
        for match in re.finditer(forward_primer.sequence, sequence.sequence)
    ]

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
            # if min_product_length is not None and product_length < min_product_length:
            #     continue
            if not desired_product_size(
                product_length, min_product_length, max_product_length
            ):
                continue

            product_line = f"{forward_primer.header}\t{reverse_primer.header}\t{start}\t{end}\t{product_length}\t{sequence.header}\t{product}"
            products.append(filter_output_line(product_line, selected_column_indices))

    results = "\n".join(products)

    if isinstance(output_file, str) is True:
        with open(output_file, "w") as fout:
            fout.write(results)

    return results


def get_pcr_products(
    primer_file: str,
    sequence_file: str,
    min_product_length: Union[int, None] = None,
    max_product_length: Union[int, None] = None,
    header: Union[bool, str] = True,
    cols: str = "all",
    output_file: Union[bool, str] = False,
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

    max_product_length: None | int
        If provided, only return those products whose length are less than or equal to this number.
        Defaults to None, which returns all products found.

    header: bool
        Whether or not to print a header on the results. Defaults to True. False will not print out the header.

    cols: str
        Which columns to print out. Defaults to "all," which prints out all the columns. A string can be supplied to
        only output the strings of interest. For example, cols="fpri rpri pname" will only output the names of the forward
        primer, reverse primer, and the target sequence when a target is found.
        Available options are:
            fpri - the name of the forward primer
            rpri - the name of the reverse primer
            start - the start location of the product in the target sequence
            end - the end location of the product in the target sequence
            length - the length of the product
            pname - the name of the sequence in which the target was found
            pseq - the nucleotide sequnce of the amplified product

    output_file: bool | str
        The file to write the results out to. Defaults to False, which will not print anything out. Providing a string
        will create that input file at that location. If set to True without providing a string, the output file
        will be of the form <DD-MM-YYYY_HH:MM:SS>.txt


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

    products = []

    # If anything gets passed for the header, it gets handled here instead of in
    # calculate_pcr_product.

    selected_column_indices = parse_selected_cols(cols)

    if header is True:
        products.append(filter_output_line(BASE_HEADER, selected_column_indices))

    primers = read_sequences_from_file(primer_file)
    forward_primer, reverse_primer = primers

    sequences = read_sequences_from_file(sequence_file)

    for sequence in sequences:
        new_products = calculate_pcr_product(
            sequence=sequence,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=min_product_length,
            max_product_length=max_product_length,
            header=False,
            cols=cols,
        )
        if new_products:
            products.append(new_products)

    results = "\n".join(products)

    if isinstance(output_file, str) is True:
        with open(output_file, "w") as fout:
            fout.write(results)

    return results
