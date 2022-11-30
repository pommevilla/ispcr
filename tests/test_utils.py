from pytest import raises

from ispcr.utils import read_fasta, reverse_complement


class TestReverseComplement:
    def test_reverse_complement(self) -> None:
        test_string = "ACGT"
        expect_reverse_complement = "ACGT"
        actual_reverse_complement = reverse_complement(test_string)

        assert actual_reverse_complement == expect_reverse_complement

    def test_incorrect_bases(self) -> None:
        test_string = "AVGT"

        with raises(KeyError):
            reverse_complement(test_string)


class TestReadFasta:
    def test_correct_number_of_entries(self) -> None:
        expected_num_lines = 20
        input_file = "tests/test_data/sequences/met_r.fa.fasta"

        actual_num_lines = 0
        with open(input_file) as fin:
            for _name, _seq in read_fasta(fin):
                actual_num_lines = actual_num_lines + 1

        assert expected_num_lines == actual_num_lines
