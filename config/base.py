from dataclasses import dataclass


@dataclass
class IntRange:
    left: int
    right: int


@dataclass
class FloatRange:
    left: int
    right: int