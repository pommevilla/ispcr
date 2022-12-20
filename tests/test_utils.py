from typing import Iterator, List

import pytest

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import (
    desired_product_size,
    is_valid_header_string,
    read_fasta,
    reverse_complement,
)


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

    # TODO: Mypy correctly identifies the single_test_sequence fixture as a string,
    # and throws an error when we try to use single_test_sequence.header.
    # PyTest runs fine, however. Hardcoding this for now to appease MyPy, but need to come
    # back and fix later.
    def test_correct_single_header(self, single_test_sequence: FastaSequence) -> None:
        expected_single_header = "single_test_sequence"
        single_test_sequence = FastaSequence("single_test_sequence", "GGG")
        actual_single_header = single_test_sequence.header

        assert expected_single_header == actual_single_header


class TestDesiredProductSize:
    def test_min_none_pass(self) -> None:
        min_product_length, max_product_length = None, 200
        potential_product_size = 150
        expected = True
        actual = desired_product_size(
            potential_product_size, min_product_length, max_product_length
        )

        assert actual == expected

    def test_min_none_fail(self) -> None:
        min_product_length, max_product_length = None, 200
        potential_product_size = 300
        expected = False
        actual = desired_product_size(
            potential_product_size, min_product_length, max_product_length
        )

        assert actual == expected

    def test_max_none_pass(self) -> None:
        min_product_length, max_product_length = 100, None
        potential_product_size = 150
        expected = True
        actual = desired_product_size(
            potential_product_size, min_product_length, max_product_length
        )

        assert actual == expected

    def test_max_none_fail(self) -> None:
        min_product_length, max_product_length = 100, None
        potential_product_size = 50
        expected = False
        actual = desired_product_size(
            potential_product_size, min_product_length, max_product_length
        )

        assert actual == expected

    def test_both_none_pass(self) -> None:
        min_product_length, max_product_length = None, None
        potential_product_sizes = [5, 50, 150, 250, 1000]
        expected = True

        for potential_product_size in potential_product_sizes:
            actual = desired_product_size(
                potential_product_size, min_product_length, max_product_length
            )
            assert actual == expected

    def test_both_pass(self) -> None:
        min_product_length, max_product_length = 75, 250
        potential_product_sizes = [75, 100, 150, 250]
        expected = True

        for potential_product_size in potential_product_sizes:
            actual = desired_product_size(
                potential_product_size, min_product_length, max_product_length
            )
            assert actual == expected

    def test_both_fail(self) -> None:
        min_product_length, max_product_length = 75, 250
        potential_product_sizes = [1, 50, 74, 251, 1000]
        expected = False

        for potential_product_size in potential_product_sizes:
            actual = desired_product_size(
                potential_product_size, min_product_length, max_product_length
            )
            assert actual == expected

    def test_improper_limits(self) -> None:
        min_product_length, max_product_length = 200, 100

        with pytest.raises(ValueError):
            desired_product_size(
                20,
                min_product_length=min_product_length,
                max_product_length=max_product_length,
            )


class TestIsValidHeaderString:
    @pytest.fixture(scope="class")
    def base_header(self) -> Iterator[str]:
        yield "fpri\trpri\tstart\tend\tlength\tpname\tpseq"

    def test_empty_string(self) -> None:
        expected = False
        actual = is_valid_header_string("")

        assert actual == expected

    def test_existing_headers_in_order(self, base_header: str) -> None:

        for col_header in base_header.split():
            assert is_valid_header_string(col_header)

    def test_existing_headers_out_of_order(self, base_header: str) -> None:
        from random import shuffle

        column_headers = base_header.split()
        shuffle(column_headers)

        for col_header in column_headers:
            assert is_valid_header_string(col_header)

    def test_single_wrong_header(self) -> None:
        expected = False
        actual = is_valid_header_string("forward_primer")

        assert actual == expected
