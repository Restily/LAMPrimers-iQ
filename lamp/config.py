class DesignConfig:

    # Range of length amplicon
    AMPLICON_RANGE = [120, 220]

    # Max diffrence of primers lengths
    MAX_DIFFERENCE_LENGTHS_PRIMERS = 3

    # Distance range between primers
    F3_F2 = [1, 10]
    F2_F1C = [10, 25]
    F1C_B1C = [0, 30] 
    B2_B1C = [10, 25]
    B3_B2 = [1, 10]

    # Distance for loop primers
    F2_F1C_LOOP = [30, 45]
    B2_B1C_LOOP = [30, 45]
    F2_LF = [2, 5]
    LF_F1C = [2, 5]
    B1C_BF = [2, 5]
    BF_B2 = [2, 5]
 
    # Distance with hybridization probe
    F1C_B1C_PROBE = [35, 50]
    F1C_PROBE = [5, 8]
    B1C_PROBE = [5, 8]


class LAMPConfig:
    """
    Constants for Temprature class
    """
    # Na+ concentration in reaction
    NA_PLUS = 0.05

    # Min difference of temperature
    TM_MAX_DIFFERENCE = 2

    # Max diffrence of temperature with probe
    TM_PROBE_MAX_VALUE = 5

    # Temperature range
    TM_RANGE = [55, 60]

    # % GC range
    GC_RANGE = [45, 55]

    # Length of primer range
    PRIMERS_LENGTH_RANGE = [24, 25]


    
