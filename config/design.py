from .base import IntRange


# Range of length amplicon
MIN_LENGTH_AMPLICON = 120
MAX_LENGTH_AMPLICON = 220

# Max diffrence of primers lengths
MAX_DIFFRENCE_LENGTHS_PRIMERS = 3

# Distance range between primers
F3_F2 = IntRange(1, 10)
F2_F1C = IntRange(10, 25)
F1C_B1C = IntRange(0, 30) 
B2_B1C = IntRange(10, 25)
B3_B2 = IntRange(1, 10)

# Distance for loop primers
F2_F1C_LOOP = IntRange(30, 45)
B2_B1C_LOOP = IntRange(30, 45)
F2_LF = IntRange(2, 5)
LF_F1C = IntRange(2, 5)
B1C_BF = IntRange(2, 5)
BF_B2 = IntRange(2, 5)

# Distance with hybridization probe
F1C_B1C_PROBE = IntRange(35, 50)
F1C_PROBE = IntRange(5, 8)
B1C_PROBE = IntRange(5, 8)
