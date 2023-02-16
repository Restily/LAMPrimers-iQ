from dataclasses import dataclass


@dataclass
class Primer:
    GC: float
    Tm: float
    seq: str
    end_idx: int
    length: int

    def convert_to_old(self) -> list:
        return [self.seq, [self.GC, self.Tm, self.end_idx, self.length]]


class PrimersSet:

    def __init__(self) -> None:
        self.primers = []

        # Iterating
        self.__cur_idx = 0

    def append(self, primer: Primer) -> None:
        self.primers.append(primer)
    
    def remove_last(self) -> None:
        self.primers.pop()
        self.__cur_idx -= 1

    def first(self) -> Primer:
        return self.primers[0]
    
    def last(self) -> Primer:
        return self.primers[-1]

    def __iter__(self) -> 'Primer':
        return self

    def __next__(self) -> Primer:
        if self.__cur_idx == len(self.primers):
            raise StopIteration
        
        primer = self.primers[self.__cur_idx]
        self.__cur_idx += 1

        return primer
    
    def __len__(self) -> int:
        return len(self.primers)