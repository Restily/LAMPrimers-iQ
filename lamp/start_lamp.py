from ctypes import *

from Bio import Seq, SeqRecord

from .config import LAMPConfig, DesignConfig


class Sets(Structure):
    _fields_ = [
        ("stringSize", c_int),
        ("setsString", c_char_p)
    ]


class LAMP:

    @classmethod
    def parse_results(cls, results: str) -> list:
        result = []
        set_primers = []

        if not results:
            return result

        lines = results.split('\n')
        for line in lines:
            if not line:
                continue

            if 'set' in line:
                if set_primers != []:
                    result.append(set_primers)
                set_primers = []
            else:
                _, primer, idx, _, primer_length, GC, Tm = line.split()
                res = [primer, [
                    round(float(GC), 2),
                    round(float(Tm), 2),
                    int(idx),
                    int(primer_length)
                ]]

                set_primers.append(res)
        
        result.append(set_primers)
        return result

    @classmethod
    def get_config_params(cls) -> list:
        DISTANCE_RANGE = [
            *DesignConfig.F3_F2,
            *DesignConfig.F2_F1C,
            *DesignConfig.F1C_B1C,
            *DesignConfig.B2_B1C,
            *DesignConfig.B3_B2,
        ]

        params = [
            c_uint8(LAMPConfig.PRIMERS_LENGTH_RANGE[0]),
            c_uint8(LAMPConfig.PRIMERS_LENGTH_RANGE[1]),
            c_float(LAMPConfig.NA_PLUS),
            (c_float * len(LAMPConfig.TM_RANGE))(*LAMPConfig.TM_RANGE),
            (c_float * len(LAMPConfig.GC_RANGE))(*LAMPConfig.GC_RANGE),
            c_uint8(DesignConfig.MAX_DIFFERENCE_LENGTHS_PRIMERS),
            c_uint8(LAMPConfig.TM_MAX_DIFFERENCE),
            (c_int * len(DesignConfig.AMPLICON_RANGE))(*DesignConfig.AMPLICON_RANGE),
            (c_int * len(DISTANCE_RANGE))(*DISTANCE_RANGE)
        ]
        
        return params

    @classmethod
    def config_dll(cls) -> CDLL:
        lib = CDLL('./LAMP.dll')

        lib.lamp.argtypes = (
            POINTER(c_char), c_uint32, c_uint8,
            c_uint8, c_float, POINTER(c_float),
            POINTER(c_float), c_uint8, c_uint8,
            POINTER(c_int), POINTER(c_int)
        )

        lib.lamp.restype = Sets

        return lib

    @classmethod
    def start_design_primers(cls, record: SeqRecord) -> list:
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
        seq = str(record.seq).encode('utf-8')
        seq_size = c_uint32(len(seq))

        # Indexes of sequence ends for search
        seq_length = [seq_length[0] - 1,
                      len(seq) if seq_length[1] == 'all' else int(seq_length[1])]
        
        lib = cls.config_dll()

        config_params = cls.get_config_params()
        params = [seq, seq_size] + config_params

        lib_results = lib.lamp(*params)
        results = lib_results.setsString[:lib_results.stringSize].decode('utf-8')

        return cls.parse_results(results)
