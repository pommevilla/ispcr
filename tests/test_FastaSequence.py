from typing import Iterator

import pytest

from ispcr.FastaSequence import FastaSequence


class TestFastSequenceDunders:
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
        expected_str = ">test1\nATTGAC\n"
        actual_str = str(test_sequence_1)

        assert expected_str == actual_str

    def test_repr(self, test_sequence_2: FastaSequence) -> None:
        expected_repr = "FastaSeq(header=test2..., seq=GTCTA...)"
        actual_repr = repr(test_sequence_2)

        assert expected_repr == actual_repr


class TestFastaSequenceDerivedFeatures:
    @pytest.fixture(scope="class")
    def test_sequence_1(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test1", "ATTGAC")

    @pytest.fixture(scope="class")
    def test_sequence_2(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test2", "GTCTAGTACAC")

    @pytest.fixture(scope="class")
    def test_sequence_3(self) -> Iterator[FastaSequence]:
        yield FastaSequence("test3", "GAGCATGCTATGTCG")

    def test_simple_gc_content(self) -> None:
        test_fasta_sequence = FastaSequence("inner_sequence", "GCAA")
        expected_gc_content = 0.5
        actual_gc_content = test_fasta_sequence.gc_content

        assert expected_gc_content == actual_gc_content

    def test_simple_base_count(self) -> None:
        test_fasta_sequence = FastaSequence("inner_sequence", "GCAT")
        expected_base_counts = {k: 1 for k in ["A", "C", "G", "T"]}
        actual_base_counts = test_fasta_sequence.base_counts

        assert expected_base_counts == actual_base_counts

    def test_simple_base_count_1(self, test_sequence_1: FastaSequence) -> None:
        expected_base_counts = {"A": 2, "C": 1, "G": 1, "T": 2}
        actual_base_counts = test_sequence_1.base_counts

        assert expected_base_counts == actual_base_counts

    def test_calc_melting_temp_short_sequences(
        self, test_sequence_1: FastaSequence, test_sequence_2: FastaSequence
    ) -> None:
        expected = [16, 32]
        actual = [seq.melting_temp for seq in [test_sequence_1, test_sequence_2]]

        assert expected == actual

    def test_calc_melting_temp_long_sequence(
        self, test_sequence_3: FastaSequence
    ) -> None:
        expected = 41.9
        actual = test_sequence_3.melting_temp

        assert expected == actual
