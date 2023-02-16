from Bio import Seq, SeqRecord

from lamp import LAMP


def start_tests():
    with open('tests/genome.txt') as f:
        seq = f.read().replace('\n', '')
        
    record = SeqRecord.SeqRecord(
        seq=Seq.Seq(seq),
        id='LAMPrimers-IQ'
    )

    LAMP.design_primers(record)