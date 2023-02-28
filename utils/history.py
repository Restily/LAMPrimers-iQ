import os
import json
from datetime import datetime


class History:
    LOGS_PATH = '/history.txt'

    @classmethod
    def check_exists(cls):
        if not os.path.exists(cls.LOGS_PATH):
            with open(cls.LOGS_PATH, 'w') as f:
                data = []
                f.write(json.dumps(data))

    @classmethod
    def write_logs(cls, primers_sets: list[list]):
        with open(cls.LOGS_PATH) as f:
            history_data = json.load(f)

        date = datetime.now()

        data = {
            'date': date,
            'results': []
        }
        # primers = [primer, [%GC, Tm, ind, primer_length]]

        for primer_set in primers_sets:
            primers = []
            for primer in primer_set:
                primer_dict = {
                    'seq': primer[0],
                    'GC': primer[1][0],
                    'Tm': primer[1][1],
                    'index': primer[1][2],
                    'length': primer[1][3]
                }
                primers.append(primer_dict)
            data['results'].append(primers)

        history_data.append(data)

        with open(cls.LOGS_PATH, 'w') as f:
            f.write(json.dumps(history_data))