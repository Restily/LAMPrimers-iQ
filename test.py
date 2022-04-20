"""
primers = [primer, [%GC, Tm, ind, primer_length]]

primer_set = [F3, F2, F1c, B1c, B2, B3]
"""
import itertools
from turtle import right

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
    

    def _check_undesigned_primers(self, primers: list[list[str, list[int, float]]], 
                      cur_primer: list[str, list[int, float]],
                      distance: list[int]):
        """
        Checking parameters for characteristics: Tm (temperature), %GC, length of primers

        :param primers: primers
        :param cur_primer: primer for checking

        :return: int: 
            -1, if the primer does not fit the min distance
            0, if the primer does not fit the max distance
            1, if the primer does not fit the parameters
            2, if primer suitable
        """
        # Flag if current primer is last
        last_primer_flag = len(primers) == 5

        # Check amplicon length range
        if not self._check_amplicon_length(primers[0][1], cur_primer[1], last_primer_flag):
            return False

        primers_distance = cur_primer[1][2] - cur_primer[1][3] - primers[-1][1][2]

        # Check max primers distance
        if primers_distance > distance[1]:
            return False

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
                            return True
        
        return 'n'


    def search_sets(self, main_primers: list[list[str, list[int, float]]], 
                    compl_primers: list[list[str, list[int, float]]]): 
        """
        Main function for build sets of primers

        :param primers: primers from main sequence
        :param compl_primers: primers from complementary sequence

        :return: sets of primers
        """
        # Count of sets primers
        count = 0

        # Sets of primers
        sets_primers = []

        # Start indexes of searching
        F3_next_ind = 0
        F2_next_ind = 1
        F1C_next_ind = 2
        B1C_next_ind = 0
        B2_next_ind = 1
        B3_next_ind = 2

        # у зонда температура на 5 градусов меньше, чем в наборе
        len_compl_primers = len(compl_primers)
        len_main_primers = len(main_primers)    
        
        # Starting seatching sets
        # Step 0 - F3
        for F3_ind in range(F3_next_ind, len_main_primers):
            print(F3_ind, len_main_primers)
            unique = True

            # Add F3 primer
            primers = [main_primers[F3_ind]]
            
            # Get F2 primer
            F2_new_ind = max(F2_next_ind, F3_ind + 1)

            for F2_ind in range(F2_new_ind, len_main_primers): #step 1 - F2
                if not unique:
                    break

                # Check primers      
                check_F3_F2 = self._check_undesigned_primers(primers, main_primers[F2_ind], DesignConfig.F3_F2)

                # If primer is not suitable, delete him 
                if not check_F3_F2:
                    if len(primers) == 1:
                        del primers[-1]

                    break

                # If we get good primer F2 for F3
                if check_F3_F2 != 'n':
                    
                    # Add F2 primer
                    primers.append(main_primers[F2_ind])
                    
                    for F1C_ind in range(F1C_next_ind, len_compl_primers): #step 2 - F1c
                        if not unique:
                            break

                        # Check primers
                        check_F2_F1C = self._check_undesigned_primers(primers, compl_primers[F1C_ind], DesignConfig.F2_F1C)
                        
                        # If primer is not suitable, delete him 
                        if not check_F2_F1C:
                            if len(primers) == 2:
                                del primers[-1]
                                
                            break
                        
                        if check_F2_F1C != 'n':
                            primers.append(compl_primers[F1C_ind])

                            for B1C_ind in range(B1C_next_ind, len_main_primers): #step 3 - B1c
                                if not unique:
                                    break

                                check_F1C_B1C = self._check_undesigned_primers(primers, main_primers[B1C_ind], DesignConfig.F1C_B1C)
                                
                                if not check_F1C_B1C:
                                    if len(primers) == 3:
                                        del primers[-1]

                                    break
                                
                                if check_F1C_B1C != 'n':
                                    primers.append(main_primers[B1C_ind])
                                    
                                    for B2_ind in range(B2_next_ind, len_compl_primers): #step 4 - B2
                                        if not unique:
                                            break

                                        check_B2_B1C = self._check_undesigned_primers(primers, compl_primers[B2_ind], DesignConfig.B2_B1C)
                                    
                                        if not check_B2_B1C:
                                            if len(primers) == 4:
                                                del primers[-1]

                                            break
                                        
                                        if check_B2_B1C != 'n':
                                            primers.append(compl_primers[B2_ind])

                                            B3_new_ind = max(B3_next_ind, B2_ind + 1)

                                            for B3_ind in range(B3_new_ind, len_compl_primers): #step 5 - B3
                                                if not unique:
                                                    break

                                                check_B3_B2 = self._check_undesigned_primers(primers, compl_primers[B3_ind], DesignConfig.B3_B2)
                                                
                                                if not check_B3_B2:
                                                    if len(primers) == 5:
                                                        del primers[-1]

                                                    break
                                                
                                                if check_B3_B2 != 'n':                                                    

                                                    # Add last B3 primer
                                                    primers.append(compl_primers[B3_ind])

                                                    # Add set of primers in list of sets
                                                    sets_primers.append(primers[:]) #'F3 ', 'F2 ', 'F1c', 'B1c', 'B2 ', 'B3 '
                                                    count += 1

                                                    # Delete last primer
                                                    primers = []

                                                    F3_next_ind = B1C_ind + 1
                                                    F1C_next_ind = B3_ind + 1

                                                    # Flag that set of primers is unique
                                                    unique = False

                                                    break
        
        return sets_primers
    


















    def _check_loop_primers(self, primer_set: list[list[list]],
                            cur_primer: list,
                            compl: bool = False):
        if compl:
            left_primers_distance = cur_primer[1][2] - cur_primer[1][3] - primer_set[1][1][2]
            right_primers_distance = primer_set[2][1][2] - (cur_primer[1][2] + cur_primer[1][3])

            left_distance = DesignConfig.F2_LF
            right_distance = DesignConfig.LF_F1C
        else:
            left_primers_distance = cur_primer[1][2] - cur_primer[1][3] - primer_set[-3][1][2]
            right_primers_distance = primer_set[-2][1][2] - (cur_primer[1][2] + cur_primer[1][3])

            left_distance = DesignConfig.BF_B2
            right_distance = DesignConfig.B2_B1C
                
         # Check max primers distance
        if right_primers_distance < right_distance[0]:
            return 0

        # Check min primers distance
        if left_primers_distance < left_distance[0] \
            or left_primers_distance > left_distance[1] \
            or right_primers_distance > right_distance[1]:
            return -1

        # Check every primer in primer_set
        for primer in primer_set:
            # Checking equality of current primer and last primer from set
            if primer[-1][0] == cur_primer[0]:
                return 1

            # Checking for diffrence lengths
            if abs(primer[1][3] - cur_primer[1][3]) > DesignConfig.MAX_DIFFRENCE_LENGTHS_PRIMERS:
                return 1

        # Check temperature diffrence of primers
        if Temperature.check_temperature_diffrence(primer_set, cur_primer):

            # Check dimers of primers
            if not Dimers.check_primers_dimers(primer_set, cur_primer[0]):
                return 2
    
        return 1


    def _search_loop_primers_ps(self, primer_set: list[list], 
                                primers: list[list], 
                                compl_primers: list[list]):
        """ 
        Design loop primers (LF and LB) for set of primers

        LoopF - from complementary sequence between F2c and F1c
        LoopB - from main sequence between B1c and B2c

        :primer_set: set of primers
        :param primers: primers from main sequence (5' - 3')
        :param compl_primers: primers from complementary sequence 

        :return: sets of primers
        """
        LoopF_primers = []
        LoopB_primers = []

        # Check LF primers
        for compl_primer_ind in range(len(compl_primers)):
            check_loop_primer = self._check_loop_primers(primer_set,
                                                         compl_primers[compl_primer_ind],
                                                         compl=True)

            if check_loop_primer == 0:
                break

            elif check_loop_primer == 2:
                LoopF_primers.append(primers[primer_ind])

        # Check BF primers
        for primer_ind in range(len(primers)):
            check_loop_primer = self._check_loop_primers(primer_set,
                                                         primers[primer_ind])
                                                         
            if check_loop_primer == 0:
                break

            elif check_loop_primer == 2:
                LoopB_primers.append(primers[primer_ind])

        print(LoopF_primers, LoopB_primers)

        loop_pairs = []
        best_pair = []

        # Build pairs and check on dimers
        for loop_pair in itertools.product(LoopF_primers, LoopB_primers):
            if not Dimers.check_primers_dimers(loop_pair[0], loop_pair[1]):
                loop_pairs.append(loop_pair)
        
        # Sort Loop primers by temperature
        if loop_pairs:
            best_pair = loop_pairs.sort(key=lambda pair: (pair[0][1][1] + pair[1][1][1]) / 2)[0]

        return best_pair


    def design_loop_primers(self, primer_sets: list[list[list]], 
                            primers: list[list], 
                            compl_primers: list[list]):
        loop_primer_sets = []

        for set_ind in range(len(primer_sets)):
            loop_primers = self._search_loop_primers_ps(primer_sets[set_ind],
                                                        primers,
                                                        compl_primers)

            if loop_primers != []:
                loop_primer_sets.append((
                    primer_sets[set_ind][:2] +
                    [loop_primers[0]] +
                    primer_sets[set_ind][2:4] +
                    [loop_primers[1]] +
                    primer_sets[set_ind][4:]
                ))

        return loop_primer_sets
            

    # def _check_hybridization_probe(self,
    #                                primer_set: list[list],
    #                                probe: )


    # def design_hybridization_probe(self,
    #                                primer_set: list[list],
    #                                sequence: str):
    #     F1c_end_ind = primer_set[2][1][2] + primer_set[2][1][3] - 1
    #     B1c_start_ind = primer_set[3][1][2]

    #     indices_pairs = itertools.product(DesignConfig.F1C_PROBE,
    #                                       DesignConfig.B1C_PROBE)

    #     inds_pairs = [[F1c_end_ind + indp[0], 
    #                    B1c_start_ind - indp[1]] for indp in indices_pairs]
        
        

# d = Design()
# d.design_hybridization_probe([[], [], [0, [0, 0, 190, 25]], [0, [0, 0, 255, 24]]], 'A')

d = Design()
primer_set = [
    ['TGGCTACTACCGAAGAGCTACCAGA', [52.0, 58.22, 28525, 25]], 
    ['CGAATTCGTGGTGGTGACGGTAA', [52.17, 56.29, 28550, 23]], 
    ['GGCCCAGTTCCTAGGTAGTAGAAA', [50.0, 56.44, 28600, 24]], 
    ['AGAAGCTGGACTTCCCTATGGTG', [52.17, 56.29, 28624, 23]], 
    ['GGTGTATTCAAGGCTCCCTCAGTT', [50.0, 56.44, 28675, 24]], 
    ['GGATTGCGGGTGCCAATGTGAT', [54.55, 56.13, 28704, 22]]
]
d._check_loop_primers(primer_set, ['AGCATGCATGTAGCTAGCTAC', [52, 58, 28575, 24]], compl=True)