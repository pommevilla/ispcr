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

    @pytest.fixture(scope="class")
    def medium_sequence_1(self) -> Iterator[FastaSequence]:
        seq_string = "AAAACCCACCCAAAAAATATTCTTTTGCATCCACTGTCAACTTTTCACAGAAACCCATTAAGTCAGGATCCTTAAGAGTTTCCGAGTGTTCATCTGCTGATATTCCAACAACAAACTCTACCGAGTGTCTGAATTTGCTGCTTGAAAAGAGAGGAGCTTCATCGAGTCAAAACTGTTGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTATTGCCGACGCTGGTGACCATGGGATTTTATCTTTCAAAGTTTCTGACCTGAAAGAAGATATACCCTCGCAAATATCGACGGCTAAGGAAGAGTCCTTCAGTGGTAATGAAGAAGAAGAAGAAGAAGAAGGTGATGACGATGATAAGATAACCCTTCAGGATTTTGTTTGTAATGAAAAGAACCAAAAAGAAATGGGTGAACAAAGAAATGACGTAAGCTCGTCTTCTTGGGTACAAACTGAGCTGTTGTTTCTTCTCCTAAAGGGAAGTATCGGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        sequence = FastaSequence("medium_test_sequence", seq_string)
        yield sequence

    def test_simple_sequence(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = "test_forward\ttest_reverse\t0\t31\t31\ttest_sequence\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"

        forward_primer, reverse_primer = primers
        actual_results = calculate_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            header=False,
        )

        assert expected_results == actual_results

    def test_minimum_product_length(
        self, primers: List[FastaSequence], small_sequence_1: FastaSequence
    ) -> None:
        expected_results = ""

        forward_primer, reverse_primer = primers
        actual_results = calculate_pcr_product(
            sequence=small_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=100,
            header=False,
        )

        assert expected_results == actual_results

    def test_maximum_product_length(
        self, primers: List[FastaSequence], medium_sequence_1: FastaSequence
    ) -> None:
        expected_results = "test_forward\ttest_reverse\t177\t255\t78\tmedium_test_sequence\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA\ntest_forward\ttest_reverse\t528\t586\t58\tmedium_test_sequence\tGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        forward_primer, reverse_primer = primers
        actual_results = calculate_pcr_product(
            sequence=medium_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            max_product_length=100,
            header=False,
        )

        assert expected_results == actual_results

    def test_product_length_range(
        self, primers: List[FastaSequence], medium_sequence_1: FastaSequence
    ) -> None:
        expected_results = "test_forward\ttest_reverse\t177\t255\t78\tmedium_test_sequence\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA"
        forward_primer, reverse_primer = primers
        actual_results = calculate_pcr_product(
            sequence=medium_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=75,
            max_product_length=100,
            header=False,
        )

        assert expected_results == actual_results


class TestGetPCRProducts:
    def test_simple_sequence(self) -> None:
        expected_result = "forward_primer.f\treverse_primer.r\t0\t31\t31\tsmall_1\tGGAGCATGCTATGTCGTAGCTGATGCAATTA"
        actual_result = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            header=False,
        )

        assert expected_result == actual_result

    def test_min_product_length(self) -> None:
        expected_result = ""
        actual_result = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            min_product_length=100,
            header=False,
        )
        assert expected_result == actual_result

    def test_maximum_product_length(self) -> None:
        expected_results = "forward_primer.f\treverse_primer.r\t177\t255\t78\tsingle_test_sequence\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA\nforward_primer.f\treverse_primer.r\t528\t586\t58\tsingle_test_sequence\tGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        actual_result = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/single_test.fa",
            max_product_length=100,
            header=False,
        )

        assert expected_results == actual_result

    def test_product_length_range(self) -> None:
        expected_results = "forward_primer.f\treverse_primer.r\t177\t255\t78\tsingle_test_sequence\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA"
        actual_results = get_pcr_products(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/single_test.fa",
            min_product_length=75,
            max_product_length=100,
            header=False,
        )

        assert expected_results == actual_results
