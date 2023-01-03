"""
Generic fasta sequence class.
"""
from dataclasses import dataclass
from typing import Iterator, Union


@dataclass(frozen=True)
class FastaSequence:
    header: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)

    def __iter__(self) -> Iterator[str]:
        yield from self.sequence

    def __str__(self) -> str:
        return f">{self.header}\n{self.sequence}"

    def __repr__(self) -> str:
        return f"FastaSeq(header={self.header[:5]}..., seq={self.sequence[:5]}...)"

    def __contains__(self, x: str) -> bool:
        return x in self.sequence

    def __getitem__(self, i: Union[int, slice]) -> str:
        return self.sequence[i]
