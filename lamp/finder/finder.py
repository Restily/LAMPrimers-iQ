"""
Main file for design primers
"""
from Bio import Seq, SeqRecord

from lamp.struct import Primer
from ..sorting.search import Sorting
from config.finder import *

from lamp.thermodynamics.dimers import Dimers
from lamp.thermodynamics.thermodynamics import Thermodynamics


class PrimerFinder(Thermodynamics):
    """ 
    Class for design primers for LAMP
    
    Methods:
    >>> start_design_primers (main function): start design primers 
    >>> get_complementary_seq: Building a complementary DNA
    """

    def __init__(self):
        Thermodynamics.__init__(self)

        # Good primers
        self._primers = []
        self._compl_primers = []

        # Range of primers lengths
        self._lengths_primers = PRIMERS_LENGTH_RANGE

    def find_primers(self, seq: str, complementary_seq: str, seq_length: list[int]):
        """        
        Gets primer from sequence

        :param seq: sequence
        :param complementary_seq: complementary sequence
        :param seq_length: indices start/end searching
        """

        # Work with first primer
        primer = seq[seq_length[0]:seq_length[0] + self._lengths_primers[0] - 1]
        self._get_first_primer_gc_count(primer)

        # Brute force (Linear search)
        for ind in range(seq_length[0], seq_length[1] - self._lengths_primers[1]):
            for primer_length in range(self._lengths_primers[0], self._lengths_primers[1] + 1):

                # Get last nucleotide in primer
                nucl = seq[ind + primer_length - 1]

                # Checking %GC and Tm primer conditions
                GC, Tm = self._check_and_get_primer_params(nucl, primer_length)

                # If conditions suitable,
                # we get primers and checking their homodimers
                if GC and Tm:

                    # Get primer on main sequence
                    primer = seq[ind:ind + primer_length]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(primer):
                        # Add params
                        cur_primer = Primer(
                            GC=GC,
                            Tm=Tm,
                            seq=primer,
                            end_idx=ind + primer_length,
                            length=primer_length
                        )

                        self._primers.append(cur_primer)

                    # Get primer on complementary sequence
                    compl_primer = complementary_seq[ind:ind +
                                                     primer_length][::-1]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(compl_primer):
                        if not cur_primer:
                            cur_primer = Primer(
                                GC=GC,
                                Tm=Tm,
                                seq=compl_primer,
                                end_idx=ind + primer_length,
                                length=primer_length
                            )

                        self._compl_primers.append(cur_primer)

            # Check first and last nucleotides for calculation dynamic %GC
            self._check_gc(seq[ind], seq[ind + self._lengths_primers[0] - 1])

            # TODO: comment
            cur_primer = None
        
        return self._primers, self._compl_primers