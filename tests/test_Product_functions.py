from typing import Iterator, List

import pytest

from ispcr import find_pcr_product, get_pcr_product
from ispcr.FastaSequence import FastaSequence


class TestGetPCRPRoduct:
    @pytest.fixture(scope="class")
    def primers(self) -> Iterator[List[FastaSequence]]:
        forward_primer = FastaSequence("test_forward", "GGAG")
        reverse_primer = FastaSequence("test_reverse", "TAAT")
        yield [forward_primer, reverse_primer]

    @pytest.fixture(scope="class")
    def small_sequence_1(self) -> Iterator[FastaSequence]:
        sequence = FastaSequence("test_sequence", "GGAGCATGCTATGTCGTAGCTGATGCAATTA")
        yield sequence

    def test_simple_sequence(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = (
            "test_forward\ttest_reverse\t0\t31\t31\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"
        )

        forward_primer, reverse_primer = primers
        actual_results = get_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
        )

        assert expected_results == actual_results

    def test_minimum_product_length(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = ""

        forward_primer, reverse_primer = primers
        actual_results = get_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=100,
        )

        assert expected_results == actual_results


class TestFindProduct:
    def test_simple_sequence(self) -> None:
        expected_result = "forward_primer.f\treverse_primer.r\t0\t31\t31\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"
        actual_result = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
        )

        assert expected_result == actual_result

    def test_min_product_length(self) -> None:
        expected_result = "forward_primer.f\treverse_primer.r\t152\t586\t434\tGGAGCTTCATCGAGTCAAAACTGTTGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTATTGCCGACGCTGGTGACCATGGGATTTTATCTTTCAAAGTTTCTGACCTGAAAGAAGATATACCCTCGCAAATATCGACGGCTAAGGAAGAGTCCTTCAGTGGTAATGAAGAAGAAGAAGAAGAAGAAGGTGATGACGATGATAAGATAACCCTTCAGGATTTTGTTTGTAATGAAAAGAACCAAAAAGAAATGGGTGAACAAAGAAATGACGTAAGCTCGTCTTCTTGGGTACAAACTGAGCTGTTGTTTCTTCTCCTAAAGGGAAGTATCGGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA\nforward_primer.f\treverse_primer.r\t177\t586\t409\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTATTGCCGACGCTGGTGACCATGGGATTTTATCTTTCAAAGTTTCTGACCTGAAAGAAGATATACCCTCGCAAATATCGACGGCTAAGGAAGAGTCCTTCAGTGGTAATGAAGAAGAAGAAGAAGAAGAAGGTGATGACGATGATAAGATAACCCTTCAGGATTTTGTTTGTAATGAAAAGAACCAAAAAGAAATGGGTGAACAAAGAAATGACGTAAGCTCGTCTTCTTGGGTACAAACTGAGCTGTTGTTTCTTCTCCTAAAGGGAAGTATCGGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        actual_result = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/single_test.fa",
            min_product_length=150,
        )

        assert expected_result == actual_result
