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

    @pytest.fixture(scope="class")
    def medium_sequence_1(self) -> Iterator[FastaSequence]:
        seq_string = "AAAACCCACCCAAAAAATATTCTTTTGCATCCACTGTCAACTTTTCACAGAAACCCATTAAGTCAGGATCCTTAAGAGTTTCCGAGTGTTCATCTGCTGATATTCCAACAACAAACTCTACCGAGTGTCTGAATTTGCTGCTTGAAAAGAGAGGAGCTTCATCGAGTCAAAACTGTTGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTATTGCCGACGCTGGTGACCATGGGATTTTATCTTTCAAAGTTTCTGACCTGAAAGAAGATATACCCTCGCAAATATCGACGGCTAAGGAAGAGTCCTTCAGTGGTAATGAAGAAGAAGAAGAAGAAGAAGGTGATGACGATGATAAGATAACCCTTCAGGATTTTGTTTGTAATGAAAAGAACCAAAAAGAAATGGGTGAACAAAGAAATGACGTAAGCTCGTCTTCTTGGGTACAAACTGAGCTGTTGTTTCTTCTCCTAAAGGGAAGTATCGGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        sequence = FastaSequence("medium_test_sequence", seq_string)
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

    def test_maximum_product_length(
        self, primers: List[FastaSequence], medium_sequence_1: FastaSequence
    ) -> None:
        expected_results = "test_forward\ttest_reverse\t177\t255\t78\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA\ntest_forward\ttest_reverse\t528\t586\t58\tGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        forward_primer, reverse_primer = primers
        actual_results = get_pcr_product(
            sequence=medium_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            max_product_length=100,
        )

        assert expected_results == actual_results

    def test_product_length_range(
        self, primers: List[FastaSequence], medium_sequence_1: FastaSequence
    ) -> None:
        expected_results = "test_forward\ttest_reverse\t177\t255\t78\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA"
        forward_primer, reverse_primer = primers
        actual_results = get_pcr_product(
            sequence=medium_sequence_1,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            min_product_length=75,
            max_product_length=100,
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
        expected_result = ""
        actual_result = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/small_sequence.fa",
            min_product_length=100,
        )
        assert expected_result == actual_result

    def test_maximum_product_length(self) -> None:
        expected_results = "forward_primer.f\treverse_primer.r\t177\t255\t78\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA\nforward_primer.f\treverse_primer.r\t528\t586\t58\tGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        actual_result = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/single_test.fa",
            max_product_length=100,
        )

        assert expected_results == actual_result

    def test_product_length_range(self) -> None:
        expected_results = "forward_primer.f\treverse_primer.r\t177\t255\t78\tGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTA"
        actual_results = find_pcr_product(
            primer_file="tests/test_data/primers/test_primers_1.fa",
            sequence_file="tests/test_data/sequences/single_test.fa",
            min_product_length=75,
            max_product_length=100,
        )

        assert expected_results == actual_results


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

    @pytest.fixture(scope="class")
    def medium_sequence_1(self) -> Iterator[FastaSequence]:
        seq_string = "AAAACCCACCCAAAAAATATTCTTTTGCATCCACTGTCAACTTTTCACAGAAACCCATTAAGTCAGGATCCTTAAGAGTTTCCGAGTGTTCATCTGCTGATATTCCAACAACAAACTCTACCGAGTGTCTGAATTTGCTGCTTGAAAAGAGAGGAGCTTCATCGAGTCAAAACTGTTGGAGAAAGATTTCTCTTGAAGATCTTTTCTGTTCCACTTCAAACCTTCCTTCCCCTACTAAAGGGAATCTCCCAATTATTGCCGACGCTGGTGACCATGGGATTTTATCTTTCAAAGTTTCTGACCTGAAAGAAGATATACCCTCGCAAATATCGACGGCTAAGGAAGAGTCCTTCAGTGGTAATGAAGAAGAAGAAGAAGAAGAAGGTGATGACGATGATAAGATAACCCTTCAGGATTTTGTTTGTAATGAAAAGAACCAAAAAGAAATGGGTGAACAAAGAAATGACGTAAGCTCGTCTTCTTGGGTACAAACTGAGCTGTTGTTTCTTCTCCTAAAGGGAAGTATCGGGAGTAATGATACTCAAACAACACTGAGAAAACCAACCCTGTTTCTGATTCCACATTA"
        sequence = FastaSequence("medium_test_sequence", seq_string)
        yield sequence

    def test_no_header(
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

    def test_basic_header(
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
