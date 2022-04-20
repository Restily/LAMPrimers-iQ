# file: config.py

"""
File for storing data for calculating thermodynamics
"""

class LAMPConfig:
    """
    Constants for Temprature class
    """
    # Na+ concentration in reaction
    NA_PLUS = 0.05

    # Min difference of temperature
    TM_MIN_VALUE = 2

    # Max diffrence of temperature with probe
    TM_PROBE_MAX_VALUE = 5

    # Temperature range
    TM_RANGE = [55, 60]

    # % GC range
    GC_RANGE = [45, 55]

    # Length of primer range
    PRIMERS_LENGTH_RANGE = [24, 25]

    

