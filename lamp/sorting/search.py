"""
primers = [primer, [%GC, Tm, ind, primer_length]]

primer_set = [F3, F2, F1c, B1c, B2, B3]
"""
from config.sorting import *
from .utils import StepGenerator
from lamp.struct import Primer, PrimersSet

from lamp.thermodynamics.dimers import Dimers
from lamp.thermodynamics.thermodynamics import Temperature


class Sorting:

    @classmethod
    def _check_amplicon_length(cls, first_primer: Primer,
                               cur_primer: Primer, last_primer_flag: bool) -> bool:
        """
        Checking that the set is in the valid amplicon length range

        :param first_primer_params: params of first primers (primers[0])
        :param cur_primer_params: params of current primer for checking

        :return: True, if distance in range of valid range of length amplicon, else False
        """
        amplicon_length = cur_primer.end_idx - \
            first_primer.end_idx + first_primer.length

        if last_primer_flag:
            if amplicon_length < MIN_LENGTH_AMPLICON:
                return False

        if amplicon_length > MAX_LENGTH_AMPLICON:
            return False

        return True

    @classmethod
    def _check_undesigned_primers(cls, primers: PrimersSet,
                                  cur_primer: Primer,
                                  distance: IntRange) -> int:
        """
        Checking parameters for characteristics: Tm (temperature), %GC, length of primers

        :param primers: primers
        :param cur_primer: primer for checking

        :return: int: 
            -1 - if the primer doesn't fit the max distance between last primer
                or primer doesn't fit the parameters
            0 - if primer reached the max distance between last primer
            1 - if primer suitable
        """
        # Flag if current primer is last
        last_primer_flag = len(primers) == 5

        # Check amplicon length range
        if not cls._check_amplicon_length(primers.first(), cur_primer, last_primer_flag):
            return 0

        last_primer = primers.last()
        primers_distance = cur_primer.end_idx - cur_primer.length - last_primer.end_idx

        # Check max primers distance
        if primers_distance > distance.right:
            return 0

        # Check min primers distance
        if primers_distance > distance.left:

            # Checking for diffrence lengths
            if abs(last_primer.length - cur_primer.length) <= MAX_DIFFRENCE_LENGTHS_PRIMERS:

                # Check temperature diffrence of primers
                if Temperature.check_temperature_diffrence(primers, cur_primer):

                    # Check dimers of primers
                    if not Dimers.check_primers_dimers(primers, cur_primer):
                        return 1
        else:
            return 2

        return -1

    @classmethod
    def search_sets(cls, main_primers: list[Primer],
                    compl_primers: list[Primer]) -> list[PrimersSet]:
        """
        Main function for build sets of primers

        :param primers: primers from main sequence
        :param compl_primers: primers from complementary sequence

        :return: sets of primers
        """
        steps_generator = StepGenerator.build()

        # Sets of primers
        result = []

        # у зонда температура на 5 градусов меньше, чем в наборе
        len_compl_primers = len(compl_primers)
        len_main_primers = len(main_primers)

        while True:
            step = steps_generator.current()

            if step.main:
                if step.next_idx >= len_main_primers:
                    break
            else:
                if step.next_idx >= len_compl_primers:
                    break

            if step.idx == 0:
                primers = PrimersSet()
                primers.append(main_primers[step.idx])

                steps_generator.iter_next()
                continue

            primer_idx = max(
                steps_generator.previous().primer_idx + 1,
                step.next_idx
            )

            if step.main:
                cur_primer = main_primers[primer_idx]
            else:
                cur_primer = compl_primers[primer_idx]

            check_primer = cls._check_undesigned_primers(primers,
                                                         cur_primer,
                                                         step.distance)

            if check_primer == 0:
                steps_generator.iter_back()
                primers.remove_last()
                continue
            elif check_primer == 1:
                primers.append(cur_primer)

                if steps_generator.check_last():
                    result.append(primers)
                    steps_generator.reset()
                else:
                    steps_generator.iter_next()
            elif check_primer == 2:
                step.next_idx = primer_idx

            step.primer_idx += 1

        return result
