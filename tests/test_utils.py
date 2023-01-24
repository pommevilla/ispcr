from typing import Dict, Iterator, List

import pytest

from ispcr.FastaSequence import FastaSequence
from ispcr.utils import (
    InvalidColumnSelectionError,
    desired_product_size,
    filter_output_line,
    get_column_indices,
    get_invalid_bases,
    is_valid_cols_string,
    parse_selected_cols,
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
        actual = is_valid_cols_string("")

        assert actual == expected

    def test_existing_headers_in_order(self, base_header: str) -> None:

        for col_header in base_header.split():
            assert is_valid_cols_string(col_header)

    def test_existing_headers_out_of_order(self, base_header: str) -> None:
        from random import shuffle

        column_headers = base_header.split()
        shuffle(column_headers)

        for col_header in column_headers:
            assert is_valid_cols_string(col_header)

    def test_single_wrong_header(self) -> None:
        expected = False
        actual = is_valid_cols_string("forward_primer")

        assert actual == expected


class TestGetColumnIndicies:
    @pytest.fixture(scope="class")
    def base_header(self) -> Iterator[str]:
        yield "fpri\trpri\tstart\tend\tlength\tpname\tpseq"

    def test_empty_string(self) -> None:
        expected: List[int] = []
        actual = get_column_indices("")

        assert expected == actual

    def test_product_columns(self) -> None:
        expected = [5, 6]
        actual = get_column_indices("pname pseq")

        assert expected == actual

    def test_out_of_order(self) -> None:
        expected = [1, 0, 4, 3]
        actual = get_column_indices("rpri fpri length end")

        assert expected == actual

    def test_default_columns(self, base_header: str) -> None:
        expected = list(range(len(base_header.split())))
        actual = get_column_indices("fpri\trpri\tstart\tend\tlength\tpname\tpseq")

        assert actual == expected


class TestFilterOutputLine:
    @pytest.fixture(scope="class")
    def base_header(self) -> Iterator[str]:
        yield "fpri\trpri\tstart\tend\tlength\tpname\tpseq"

    @pytest.fixture(scope="class")
    def small_product_line(self) -> Iterator[str]:
        yield "forward_primer.f\treverse_primer.r\t0\t31\t31\ttest_sequence\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"

    def test_no_output_columns_small_products(self, small_product_line: str) -> None:
        expected = ""
        col_indices: List[int] = []
        actual = filter_output_line(small_product_line, column_indices=col_indices)

        assert expected == actual

    def test_checks_for_default_header(self, base_header: str) -> None:
        expected = base_header
        actual = filter_output_line(base_header, list(range(7)))

        assert expected == actual

    def test_no_early_return_when_line_header_string_out_of_order(
        self, base_header: str
    ) -> None:
        expected = "pname\trpri\tend\tlength\tstart\tpseq\tfpri"
        actual = filter_output_line(base_header, [5, 1, 3, 4, 2, 6, 0])

        assert expected == actual

    def test_product_info_only(self, small_product_line: str) -> None:
        expected = "test_sequence\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"
        actual = filter_output_line(small_product_line, [5, 6])

        assert expected == actual

    def test_returns_empty_string(self) -> None:
        expected = ""
        actual = filter_output_line("", list(range(7)))

        assert expected == actual


class TestParseSelectedCols:
    @pytest.fixture(scope="class")
    def base_header(self) -> Iterator[str]:
        yield "fpri\trpri\tstart\tend\tlength\tpname\tpseq"

    def test_exits_on_invalid_cols(self) -> None:
        with pytest.raises(InvalidColumnSelectionError):
            parse_selected_cols("invalid col string")

    def test_outputs_default_string(self, base_header: str) -> None:
        expected = list(range(len(base_header.split())))
        actual = parse_selected_cols(base_header)

        assert expected == actual


class TestGetInvalidBases:
    def test_simple_invalid_base_count(self) -> None:
        test_string = "ACGTV"
        expected = {"V": [4]}
        actual = get_invalid_bases(test_string)
        assert expected == actual

    def test_no_invalid_bases_found(self) -> None:
        test_string = "ACGT"
        expected: Dict[str, List[int]] = {}
        actual = get_invalid_bases(test_string)
        assert expected == actual

    def test_multiple_single_invalid_bases(self) -> None:
        test_string = "TZACFGH"
        expected = {"Z": [1], "F": [4], "H": [6]}
        actual = get_invalid_bases(test_string)
        assert expected == actual

    def test_multiple_times(self) -> None:
        test_string = "TZACFGH"
        expected = {"Z": [1], "F": [4], "H": [6]}
        actual = get_invalid_bases(test_string)
        assert expected == actual
