from pytest import raises

from ispcr.utils import reverse_complement


def test_reverse_complement() -> None:
    test_string = "ACGT"
    expect_reverse_complement = "ACGT"
    actual_reverse_complement = reverse_complement(test_string)

    assert actual_reverse_complement == expect_reverse_complement


def test_incorrect_bases() -> None:
    test_string = "AVGT"

    with raises(KeyError):
        reverse_complement(test_string)
