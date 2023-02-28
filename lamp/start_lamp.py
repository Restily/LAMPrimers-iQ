"""
Main file for design primers
"""
from Bio import Seq, SeqRecord

from design.search_primers import Design
from lamp.dimers import Dimers
from lamp.thermodynamics import Thermodynamics
from lamp.config import LAMPConfig
from utils.history import History


class LAMP(Thermodynamics):
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
        self._lengths_primers = LAMPConfig.PRIMERS_LENGTH_RANGE

    def _get_and_check_primers(self, seq: str, complementary_seq: str, seq_length: list[int]):
        """
        Private method 
        
        Gets primer from genome

        :param seq: sequence
        :param complementary_seq: complementary sequence
        :param seq_length: indices start/end searching
        """

        # Work with first primer
        primer = seq[seq_length[0]:seq_length[0] +
                     self._lengths_primers[0] - 1]
        self._get_first_primer_gc_count(primer)

        # Brute force (Linear search)
        for ind in range(seq_length[0], seq_length[1] - self._lengths_primers[1]):
            for primer_length in range(self._lengths_primers[0], self._lengths_primers[1] + 1):

                # Get last nucleotide in primer
                nucl = seq[ind + primer_length - 1]

                # Checking %GC and Tm primer conditions
                self._primer_params = self._check_and_get_primer_params(
                    nucl, primer_length)

                # If conditions suitable,
                # we get primers and checking their homodimers
                if self._primer_params:

                    # Get primer on main sequence
                    primer = seq[ind:ind + primer_length]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(primer):
                        # Add params
                        self._primer_params.extend([ind + 1, primer_length])

                        self._primers.append([primer, self._primer_params])

                    # Get primer on complementary sequence
                    compl_primer = complementary_seq[ind:ind +
                                                     primer_length][::-1]

                    # Check primer on homodimer
                    if not Dimers.check_homodimer(compl_primer):
                        if len(self._primer_params) != 4:
                            self._primer_params.extend(
                                [ind + 1, primer_length])

                        self._compl_primers.append(
                            [compl_primer, self._primer_params])

            # Check first and last nucleotides for calculation dynamic %GC
            self._check_gc(seq[ind], seq[ind + self._lengths_primers[0] - 1])

    def get_complementary_seq(self, seq: Seq.Seq) -> str:
        """
        Building a complementary DNA
        
        :param seq: DNA sequence (Seq.Seq from BioPython)
        
        :return: complementary DNA sequence (str)
        """
        return str(seq.complement())

    # seq_length: list[int, str]) -> list:
    def start_design_primers(self, record: SeqRecord) -> list:
        """
        Design primers (main function)

        First, we get all the suitable primers from the genome with function get_primers,
        and then we sort by sets with function set_search.

        seq ---> 5' - 3'
        complementary_seq ---> 3' - 5'
        
        :param seq_file: path to sequence file
        :param seq_length: indices of start/end sequence searching

        :return: primers_sets: selected sets of primers for LAMP
        """
        seq_length = [1, 'all']

        # Genome sequence reading
        seq = str(record.seq)

        # Indexes of sequence ends for search
        seq_length = [seq_length[0] - 1,
                      len(seq) if seq_length[1] == 'all' else int(seq_length[1])]
        # Исправить индексы

        # Get complementary sequence
        complementary_seq = self.get_complementary_seq(record.seq)

        # Get best primers from genom
        self._get_and_check_primers(seq, complementary_seq, seq_length)

        # Create design class for primer sets
        design = Design()

        print(len(self._primers))

        # Get primer sets
        sets_primers = design.search_sets(self._primers, self._compl_primers)

        # Sort primer sets
        sets_primers = self.sort_primer_sets(sets_primers)

        History.write_logs(sets_primers)

        for primer_set in sets_primers:
            print(primer_set, design.design_loop_primers(primer_set,
                                                         self._primers,
                                                         self._compl_primers), '\n')

        # Return primer sets
        return sets_primers
