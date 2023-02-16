from time import time
from Bio import Seq, SeqRecord

from . import utils
from config.finder import SEARCH_RANGE

from .finder.finder import PrimerFinder
from .design.search import Design


class LAMP:
        
    @classmethod
    def design_primers(cls, record: SeqRecord):
        # Genome sequence reading
        main_seq = str(record.seq)

        # Get complementary sequence
        complementary_seq = utils.get_complementary_seq(record.seq)

        if SEARCH_RANGE.right == -1:
            SEARCH_RANGE.right = len(main_seq)

        ftime = time()

        # Find primers
        main_primers, compl_primers = PrimerFinder().find_primers(
            main_seq, complementary_seq, SEARCH_RANGE
        )

        print('Search time:', time() - ftime)
        print(len(main_primers), len(compl_primers))

        # Get primer sets
        sets_primers = Design.search_sets(main_primers, compl_primers)

        old_sets = []
        # for primer_set in sets_primers:
        #     old_set = []
        #     for primer in primer_set:
        #         old_set.append(primer.convert_to_old())
            
        #     old_sets.append(old_set)

        # Return primer sets
        return old_sets