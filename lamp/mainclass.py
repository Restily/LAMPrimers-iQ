from Bio import Seq, SeqRecord

from . import utils
from config.finder import SEARCH_RANGE

from .finder.finder import PrimerFinder
from .sorting.search import Sorting


class LAMP:
        
    @classmethod
    def design_primers(cls, record: SeqRecord):
        # Genome sequence reading
        main_seq = str(record.seq)

        # Get complementary sequence
        complementary_seq = utils.get_complementary_seq(record.seq)

        # Find primers
        main_primers, compl_primers = PrimerFinder().find_primers(
            main_seq, complementary_seq, SEARCH_RANGE
        )

        # Get primer sets
        sets_primers = Sorting.search_sets(main_primers, compl_primers)

        old_sets = []
        for primer_set in sets_primers:
            old_set = []
            for primer in primer_set:
                old_set.append(primer.convert_to_old())
            
            old_sets.append(old_set)

        # Return primer sets
        return old_sets