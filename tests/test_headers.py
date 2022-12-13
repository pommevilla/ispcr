from typing import Iterator, List

import pytest

from ispcr import find_pcr_product, get_pcr_product
from ispcr.FastaSequence import FastaSequence


class TestHeaders:
    @pytest.fixture(scope="class")
    def primers(self) -> Iterator[List[FastaSequence]]:
        forward_primer = FastaSequence("test_forward", "GGAG")
        reverse_primer = FastaSequence("test_reverse", "TAAT")
        yield [forward_primer, reverse_primer]

    @pytest.fixture(scope="class")
    def small_sequence_1(self) -> Iterator[FastaSequence]:
        sequence = FastaSequence("test_sequence", "GGAGCATGCTATGTCGTAGCTGATGCAATTA")
        yield sequence

    def test_no_header_get_pcr_product(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = ""

        forward_primer, reverse_primer = primers
        pcr_results = get_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            header=False,
        ).splitlines()
        actual_results = "".join(
            [line for line in pcr_results if "forward_primer" in line]
        )

        assert expected_results == actual_results

    def test_basic_header_get_pcr_product(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = (
            "forward_primer\treverse_primer\tstart\tend\tlength\tsequence"
        )

        forward_primer, reverse_primer = primers
        pcr_results = get_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            header=True,
        ).splitlines()
        actual_results = "".join(
            [line for line in pcr_results if line == expected_results]
        )

        assert expected_results == actual_results

    def test_no_header_find_pcr_product(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = ""

        pcr_results = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            min_product_length=100,
            header=False,
        ).splitlines()
        actual_results = "".join(
            [line for line in pcr_results if "forward_primer" in line]
        )

        assert expected_results == actual_results

    def test_basic_header_find_pcr_product(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = (
            "forward_primer\treverse_primer\tstart\tend\tlength\tsequence"
        )

        pcr_results = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            min_product_length=100,
            header=True,
        ).splitlines()
        actual_results = "".join(
            [line for line in pcr_results if line == expected_results]
        )

        assert expected_results == actual_results
