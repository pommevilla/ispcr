from typing import Iterator

import pytest

from ispcr.FastaSequence import FastaSequence


class TestFastaSequence:
    @pytest.fixture(scope="class")
    def test_sequence_1(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test1", "ATTGAC")

    @pytest.fixture(scope="class")
    def test_sequence_2(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test2", "GTCTAGTACAC")

    def test_contains(self, test_sequence_1: FastaSequence) -> None:
        assert "TT" in test_sequence_1

    def test_does_not_contain(self, test_sequence_1: FastaSequence) -> None:
        assert "CC" not in test_sequence_1

    def test_length(self, test_sequence_1: FastaSequence) -> None:
        expected_length = len(test_sequence_1.sequence)
        actual_length = len(test_sequence_1)

        assert actual_length == expected_length

    def test_iterator(self, test_sequence_1: FastaSequence) -> None:
        base_count = 0
        for _base in test_sequence_1:
            base_count += 1

        assert base_count == len(test_sequence_1)

    def test_indexing(self, test_sequence_1: FastaSequence) -> None:
        expected_base = "G"
        actual_base = test_sequence_1[3]

        assert expected_base == actual_base

    def test_slicing(self, test_sequence_2: FastaSequence) -> None:
        expected_slice = "CAC"
        actual_slice = test_sequence_2[-3:]

        assert expected_slice == actual_slice

    def test_str(self, test_sequence_1: FastaSequence) -> None:
        expected_str = ">test1\nATTGAC"
        actual_str = str(test_sequence_1)

        assert expected_str == actual_str

    def test_repr(self, test_sequence_2: FastaSequence) -> None:
        expected_repr = "FastaSeq(header=test2..., seq=GTCTA...)"
        actual_repr = repr(test_sequence_2)

        assert expected_repr == actual_repr
