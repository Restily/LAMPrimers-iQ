from dataclasses import dataclass

from config.base import IntRange
from config.sorting import *


@dataclass
class Step:
    idx: int
    primer_idx: int
    distance: IntRange
    next_idx: int
    main: bool


class StepGenerator:

    def __init__(self) -> None:
        self.__steps = []
        self.__cur_step = 0

    def iter_next(self):
        self.__cur_step += 1

    def iter_back(self):
        self.__cur_step -= 1

    def current(self) -> Step:
        return self.__steps[self.__cur_step]

    def previous(self) -> Step:
        return self.__steps[self.__cur_step - 1]

    def add_step(self, step: Step) -> None:
        self.__steps.append(step)

    def reset(self) -> None:
        self.__cur_step = 0

    def check_last(self) -> bool:
        return self.__cur_step == len(self.__steps) - 1

    @classmethod
    def build(cls) -> 'StepGenerator':
        generator = StepGenerator()
        distances = [F3_F2, F2_F1C, F1C_B1C, B2_B1C, B3_B2]
        next_idxs = [0, 1, 0, 2, 1, 2]
        mains = [True, True, False, True, False, False]

        for idx, in range(len(distances)):
            generator.add_step(Step(
                idx=idx, 
                distance=distances[idx], 
                next_idx=next_idxs[idx],
                main=mains[idx]
            ))

        return generator
