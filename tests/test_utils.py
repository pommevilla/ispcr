from typing import Iterator, List

import pytest

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import read_fasta, reverse_complement


class TestReverseComplement:
    def test_reverse_complement(self) -> None:
        test_string = "ACGT"
        expect_reverse_complement = "ACGT"
        actual_reverse_complement = reverse_complement(test_string)

        assert actual_reverse_complement == expect_reverse_complement

    def test_incorrect_bases(self) -> None:
        test_string = "AVGT"

        with pytest.raises(KeyError):
            reverse_complement(test_string)


class TestReaderUtils:
    @pytest.fixture(scope="class")
    def single_test_sequence(self) -> Iterator[List[FastaSequence]]:
        input_file = "tests/test_data/sequences/single_test.fa"
        fasta_sequences = []
        with open(input_file) as fin:
            for fasta_sequence in read_fasta(fin):
                fasta_sequences.append(fasta_sequence)
        yield fasta_sequences

    @pytest.fixture(scope="class")
    def multiple_test_sequences(self) -> Iterator[List[FastaSequence]]:
        input_file = "tests/test_data/sequences/met_r.fa"
        fasta_sequences = []
        with open(input_file) as fin:
            for fasta_sequence in read_fasta(fin):
                fasta_sequences.append(fasta_sequence)
        yield fasta_sequences

    def test_correct_number_of_entries(
        self, multiple_test_sequences: List[FastaSequence]
    ) -> None:
        assert len(multiple_test_sequences) == 20

    def test_read_single_fasta(self, single_test_sequence: FastaSequence) -> None:
        assert len(single_test_sequence) == 1

    def test_correct_single_header(self, single_test_sequence: FastaSequence) -> None:
        expected_single_header = "single_test_sequence"
        single_test_sequence = FastaSequence("single_test_sequence", "GGG")
        actual_single_header = single_test_sequence.header

        assert expected_single_header == actual_single_header
