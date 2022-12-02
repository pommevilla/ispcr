import re

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import read_sequences_from_file, reverse_complement


def get_pcr_product(
    sequence: FastaSequence,
    forward_primer: FastaSequence,
    reverse_primer: FastaSequence,
) -> str:
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

            product_line = f"{forward_primer.header}\t{reverse_primer.header}\t{start}\t{end}\t{product_length}\t{product}"
            products.append(product_line)
    return "\n".join(products)


def find_pcr_product(primer_file: str, sequence_file: str) -> str:
    primers = read_sequences_from_file(primer_file)
    forward_primer, reverse_primer = primers

    sequences = read_sequences_from_file(sequence_file)
    products = []

    for sequence in sequences:
        new_products = get_pcr_product(
            sequence=sequence,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
        )
        products.append(new_products)
    return "\n".join(products)
