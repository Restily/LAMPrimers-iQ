# file: dimers.py
from lamp.struct import PrimersSet, Primer


class Dimers:
    """
    Class for checking of Dimers primers
    """

    def __init__(self):
        super()

    @staticmethod
    def check_complementary_nucls(fnucl: str, snucl: str) -> bool:
        """
        Check nucleotides are they complementary

        :param fnucl: nucleotide
        :param snucl: nucleotide

        :return: True if complementary, else False
        """
        compl_pairs = {'A': 'T', 'T': 'A', 'G': 'C',
                       'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

        if fnucl in compl_pairs:
            if snucl == compl_pairs[fnucl]:  # Check nucls complementary
                return True

        return False

    @staticmethod
    def compare_primers(first_primer: str, second_primer: str) -> list[str, int]:
        """
        Compares two primers and checks how much complementary nucls they have

        :param first_primer: first primer (from 5' to 3') (or opposite)
        :param second_primer: second primer (from 3' to 5') (or opposite)
        
        :return: list[String: matches, Integer: count of complementary pairs]
        """
        pairs = 0

        for i in range(len(first_primer)):  # Get count of complementary nucls
            if Dimers.check_complementary_nucls(first_primer[i], second_primer[i]):
                pairs += 1

        return pairs

    @staticmethod
    def check_dimert(first_primer: str, second_primer: str) -> bool:
        """
        Function for identifying dimers of two primers

        :param first_primer: first sequence (from 5' to 3')
        :param second_primer: second sequence (from 5' to 3')

        :return: True, if sequences have dimer, else False
        """

        second_primer = second_primer[::-1]

        if Dimers.compare_primers(first_primer[-3:], second_primer[:3]) == 3:
            return True

        for i in range(1, max(len(first_primer), len(second_primer))):
            if 3 + i <= len(second_primer):
                if Dimers.compare_primers(first_primer[-3:], second_primer[i: 3 + i]) == 3:
                    return [0, i]

            if 3 + i <= len(first_primer):
                if Dimers.compare_primers(second_primer[:3], first_primer[-3 - i: -i]) == 3:
                    return [1, i]

        return False
    
    @staticmethod
    def check_dimer(first_primer: str, second_primer: str) -> bool:
        """
        Function for identifying dimers of two primers

        :param first_primer: first sequence (from 5' to 3')
        :param second_primer: second sequence (from 5' to 3')

        :return: True, if sequences have dimer, else False
        """

        second_primer = second_primer[::-1]

        if Dimers.compare_primers(first_primer[-3:], second_primer[:3]) == 3:
            return True

        for i in range(1, max(len(first_primer), len(second_primer))):
            if 3 + i <= len(second_primer):
                if Dimers.compare_primers(first_primer[-3:], second_primer[i: 3 + i]) == 3:
                    return True

            # if 3 + i <= len(first_primer):
            #     if Dimers.compare_primers(second_primer[:3], first_primer[-3 - i: -i]) == 3:
            #         return True

        return False

    @staticmethod
    def check_homodimer(primer: str) -> bool:
        """
        Check dimers in primer

        :param primer: primer (from 5' to 3')

        :return: True if have dimers else False
        """
        if primer[-2:] == 'GC':
            return True

        if 'AAAA' in primer or 'GGGG' in primer or 'CCCC' in primer or 'TTTT' in primer:
            return True

        return Dimers.check_dimer(primer, primer)

    @staticmethod
    def check_primers_dimers(primers: PrimersSet, cur_primer: Primer) -> bool:
        """
        Checks if primer can form at least one dimer with some primer from
        set of primers without dimers

        :param primers: set of primers without dimers
        :param primer: primer

        :return: True if set have dimer with primer else False
        """
        for primer in primers:
            if Dimers.check_dimer(primer.seq, cur_primer.seq):
                return True

        return False
