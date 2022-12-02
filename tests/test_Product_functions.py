from typing import Iterator, List

import pytest

from ispcr import get_pcr_product
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
