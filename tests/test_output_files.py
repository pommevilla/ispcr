from os import listdir, remove
from os.path import exists
from typing import Iterator, List

import pytest

from ispcr import calculate_pcr_product, get_pcr_products
from ispcr.FastaSequence import FastaSequence


class TestCalculatePCRPRoduct:
    @pytest.fixture(scope="class")
    def primers(self) -> Iterator[List[FastaSequence]]:
        forward_primer = FastaSequence("test_forward", "GGAG")
        reverse_primer = FastaSequence("test_reverse", "TAAT")
        yield [forward_primer, reverse_primer]

    @pytest.fixture(scope="class")
    def small_sequence_1(self) -> Iterator[FastaSequence]:
        sequence = FastaSequence("test_sequence", "GGAGCATGCTATGTCGTAGCTGATGCAATTA")
        yield sequence

    def test_no_output_file_when_false(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        cwd_before = listdir()
        forward_primer, reverse_primer = primers
        _ = calculate_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            header=False,
            output_file=False,
        )
        cwd_after = listdir()

        assert cwd_before == cwd_after

    def test_output_file_with_name(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        output_file_name = "test.txt"
        forward_primer, reverse_primer = primers
        _ = calculate_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            header=False,
            output_file=output_file_name,
        )
        assert exists(output_file_name)
        remove(output_file_name)


class TestGetPCRProducts:
    def test_no_output_file_when_false(self) -> None:
        cwd_before = listdir()
        _ = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            output_file=False,
        )
        cwd_after = listdir()

        assert cwd_before == cwd_after

    def test_output_file_with_name(self) -> None:
        output_file_name = "test.txt"
        _ = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            output_file=output_file_name,
        )
        assert exists(output_file_name)
        remove(output_file_name)
