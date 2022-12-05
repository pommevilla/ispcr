# ispcr

[![PyPI](https://img.shields.io/pypi/v/ispcr?style=flat-square)](https://pypi.python.org/pypi/ispcr/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ispcr?style=flat-square)](https://pypi.python.org/pypi/ispcr/)
[![PyPI - License](https://img.shields.io/pypi/l/ispcr?style=flat-square)](https://pypi.python.org/pypi/ispcr/)
[![Coookiecutter - Wolt](https://img.shields.io/badge/cookiecutter-Wolt-00c2e8?style=flat-square&logo=cookiecutter&logoColor=D4AA00&link=https://github.com/woltapp/wolt-python-package-cookiecutter)](https://github.com/woltapp/wolt-python-package-cookiecutter)


---

**Documentation**: [https://pommevilla.github.io/ispcr](https://pommevilla.github.io/ispcr)

**Source Code**: [https://github.com/pommevilla/ispcr](https://github.com/pommevilla/ispcr)

**PyPI**: [https://pypi.org/project/ispcr/](https://pypi.org/project/ispcr/)

---

A simple, light-weight package written in base Python to perform *in silico* PCR to determine primer performance.

**Currently in development**

## Installation

```sh
pip install ispcr
```
## Demonstration

The main function to use in this package is `find_pcr_product`, which takes as input two file paths:
  * `primer_file` - the path to fasta file containing your primers
    * This is currently limited to a fasta file containing two sequences, with the forward primer coming first and the reverse primer coming second
  * `sequence_file` the path to the fasta file containing the sequences to test your primers against

`find_pcr_product` will then iterate through the sequences in `sequence` file and find all products amplified by the forward and reverse primer.

![](imgs/find_pcr_product.png)

`find_pcr_product` also takes a `minimum_product_length` argument:

![](imgs/find_pcr_product_min_length.png)
