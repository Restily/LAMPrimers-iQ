from Bio import Seq


def get_complementary_seq(seq: Seq.Seq) -> str:
    """
    Building a complementary DNA
    
    :param seq: DNA sequence (Seq.Seq from BioPython)
    
    :return: complementary DNA sequence (str)
    """
    return str(seq.complement())