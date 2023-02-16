"""
Main file for design primers
"""
from Bio import Seq, SeqRecord

from lamp.struct import Primer
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
        primer = seq[seq_length.left:seq_length.left + self._lengths_primers.left - 1]
        self._get_first_primer_gc_count(primer)

        # Brute force (Linear search)
        for idx in range(seq_length.left, seq_length.right - self._lengths_primers.right):
            for primer_length in range(self._lengths_primers.left, self._lengths_primers.right + 1):

                # Get last nucleotide in primer
                nucl = seq[idx + primer_length - 1]

                # Checking %GC and Tm primer conditions
                params = self._check_and_get_primer_params(nucl, primer_length)

                # If conditions suitable,
                # we get primers and checking their homodimers
                if params:
                    GC, Tm = params
                    end_idx = idx + primer_length

                    # Get primer on main sequence
                    primer = seq[idx:end_idx]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(primer):
                        # Add params
                        cur_primer = Primer(
                            GC=GC,
                            Tm=Tm,
                            seq=primer,
                            end_idx=idx + primer_length,
                            length=primer_length
                        )

                        self._primers.append(cur_primer)

                    # Get primer on complementary sequence
                    compl_primer = complementary_seq[idx:end_idx][::-1]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(compl_primer):
                        if not cur_primer:
                            cur_primer = Primer(
                                GC=GC,
                                Tm=Tm,
                                seq=compl_primer,
                                end_idx=idx + primer_length,
                                length=primer_length
                            )

                        self._compl_primers.append(cur_primer)

            # Check first and last nucleotides for calculation dynamic %GC
            self._check_gc(seq[idx], seq[idx + self._lengths_primers.left - 1])

            # TODO: comment
            cur_primer = None
        
        return self._primers, self._compl_primers