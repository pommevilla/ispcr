"""
Generic fasta sequence class.
"""
from dataclasses import dataclass, field
from typing import Iterator, Union

VALID_BASES = set("ACGT")


@dataclass(frozen=True)
class FastaSequence:
    header: str
    sequence: str
    base_counts: dict[str, int] = field(init=False)
    gc_content: float = field(init=False)
    melting_temp: Union[float, int] = field(init=False)

    def __post_init__(self) -> None:
        from collections import Counter

        base_counts = Counter(self.sequence)
        object.__setattr__(self, "base_counts", base_counts)

        gc_count = self.base_counts["G"] + self.base_counts["C"]
        gc_content = gc_count / len(self.sequence)
        object.__setattr__(self, "gc_content", gc_content)

        melting_temp = self._get_melting_temp()
        object.__setattr__(self, "melting_temp", melting_temp)

    def __len__(self) -> int:
        return len(self.sequence)

    def __iter__(self) -> Iterator[str]:
        yield from self.sequence

    def __str__(self) -> str:
        return f">{self.header}\n{self.sequence}\n"

    def __repr__(self) -> str:
        return f"FastaSeq(header={self.header[:5]}..., seq={self.sequence[:5]}...)"

    def __contains__(self, x: str) -> bool:
        return x in self.sequence

    def __getitem__(self, i: Union[int, slice]) -> str:
        return self.sequence[i]

    def _get_melting_temp(self) -> Union[float, int]:

        melting_temp: Union[float, int]

        if len(self.sequence) < 14:
            melting_temp = 2 * (self.base_counts["A"] + self.base_counts["T"]) + 4 * (
                self.base_counts["G"] + self.base_counts["C"]
            )
        else:
            melting_temp = 64.9 + 41 * (
                self.base_counts["G"] + self.base_counts["C"] - 16.4
            ) / len(self.sequence)

        melting_temp = round(melting_temp, 1)

        return melting_temp
