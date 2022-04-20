"""
primers = [primer, [%GC, Tm, ind, primer_length]]
"""
from xml.etree.ElementInclude import include
from lamp.dimers import Dimers
from lamp.thermodynamics import Temperature
from design.config import DesignConfig


class Design:

    def __init__(self):
        super()
        

    def _check_amplicon_length(self, first_primer_params: list[int, float], 
                              cur_primer_params: list[int, float], last_primer_flag: bool) -> bool: 
        """
        Checking that the set is in the valid amplicon length range

        :param first_primer_params: params of first primers (primers[0])
        :param cur_primer_params: params of current primer for checking

        :return: True, if distance in range of valid range of length amplicon, else False
        """
        amplicon_length = cur_primer_params[2] - first_primer_params[2] + first_primer_params[3]

        if last_primer_flag:
            if amplicon_length < DesignConfig.MIN_LENGTH_AMPLICON:
                return False

        if amplicon_length > DesignConfig.MAX_LENGTH_AMPLICON:
            return False

        return True  
    

    def _check_primers(self, primers: list[list[str, list[int, float]]], 
                      cur_primer: list[str, list[int, float]],
                      distance: list[int]):
        """
        Checking parameters for characteristics: Tm (temperature), %GC, length of primers

        :param primers: primers
        :param cur_primer: primer for checking

        :return: int: 
            -1, if the primer does not fit the max distance or amplicon length not suitable
            0, if the primer does not fit the min distance
            1, if the primer does not fit the parameters
            2, if primer suitable
        """
        # Flag if current primer is last
        last_primer_flag = len(primers) == 5

        # Check amplicon length range
        if not self._check_amplicon_length(primers[0][1], cur_primer[1], last_primer_flag):
            return -1

        primers_distance = cur_primer[1][2] - cur_primer[1][3] - primers[-1][1][2]

        # Check max primers distance
        if primers_distance > distance[1]:
            return -1

        # Check min primers distance
        if primers_distance > distance[0]:

            # Checking equality of current primer and last primer from set
            if primers[-1][0] != cur_primer[0]:

                # Checking for diffrence lengths
                if abs(primers[-1][1][3] - cur_primer[1][3]) <= DesignConfig.MAX_DIFFRENCE_LENGTHS_PRIMERS:

                    # Check temperature diffrence of primers
                    if Temperature.check_temperature_diffrence(primers, cur_primer):

                        # Check dimers of primers
                        if not Dimers.check_primers_dimers(primers, cur_primer[0]):
                            return 2

            return 1
        
        return 0


    def search_sets(self, main_primers: list[list[str, list[float, float, int, int]]], 
                    compl_primers: list[list[str, list[int, float, int, int]]]) -> tuple: #primers - 5'-3', compl_primers - 3'-5'
        """
        Main function for build sets of primers

        :param primers: primers from main sequence (5' - 3')
        :param compl_primers: primers from complementary sequence (3' - 5')

        :return: sets of primers
        """
        distances = [
            DesignConfig.F3_F2,
            DesignConfig.F2_F1C,
            DesignConfig.F1C_B1C,
            DesignConfig.B2_B1C,
            DesignConfig.B3_B2
        ]

        # F3, F2, F1c, B1c, B2, B3
        indices = [0, 1, 0, 2, 1, 2]
        flags = []
        cur_ind = 0

        len_compl_primers = len(compl_primers)
        len_main_primers = len(main_primers)    

        sets_primers = []
        primers = []

        while indices[0] != len_main_primers - 3:
            if cur_ind == 0:
                flags = [True, True, True, True, True]
                primers.append(main_primers[indices[0]])

                cur_ind += 1
                indices[0] += 1

                new_indices = indices[:]

            if (cur_ind % 5) % 2 == 1:
                cur_primers = main_primers
                cur_len_primers = len_main_primers
            else:
                cur_primers = compl_primers
                cur_len_primers = len_compl_primers

            if cur_ind == 3 or cur_ind == 4:
                indices[cur_ind] = max(indices[cur_ind], indices[cur_ind - 2] + 1)
            elif cur_ind != 0 and cur_ind != 2:
                indices[cur_ind] = max(indices[cur_ind], indices[cur_ind - 1] + 1)

            for primer_index in range(indices[cur_ind], cur_len_primers):                
                check_current_primers = self._check_primers(primers, 
                                        cur_primers[primer_index], distances[cur_ind - 1])

                print(cur_ind, indices[cur_ind], primer_index, check_current_primers)
                
                if check_current_primers == -1:
                    if len(primers) == cur_ind + 1:
                        del primers[-1]
                        
                    cur_ind -= 1

                    break
                elif check_current_primers == 0:
                    indices[cur_ind] = primer_index
                elif check_current_primers >= 1 and flags[cur_ind - 1]:
                    indices[cur_ind] = primer_index
                    flags[cur_ind - 1] = False
                else:
                    primers.append(cur_primers[primer_index])

                    new_indices[cur_ind] = primer_index + 1
                    cur_ind += 1

                    break
            
            if cur_ind == 6:
                sets_primers.append(primers)

                indices[0] = indices[3] + 1
                indices[2] = indices[5] + 1

                primers = []
                cur_ind = 0
        
        return sets_primers



