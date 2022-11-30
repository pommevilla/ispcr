from typing import Iterator

import pytest

from ispcr.FastaSequence import FastaSequence


class TestFastaSequence:
    @pytest.fixture(scope="class")
    def test_sequence_1(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test", "ATTGAC")

    @pytest.fixture(scope="class")
    def test_sequence_2(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test", "GTCTAGTACAC")

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
