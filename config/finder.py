from .base import IntRange, FloatRange


SEARCH_RANGE = IntRange(0, 0)

# Na+ concentration in reaction
NA_PLUS = 0.05

# Min difference of temperature
TM_MIN_VALUE = 2

# Max diffrence of temperature with probe
TM_PROBE_MAX_VALUE = 5

# Temperature range
TM_RANGE = FloatRange(55, 60)

# % GC range
GC_RANGE = FloatRange(45, 55)

# Length of primer range
PRIMERS_LENGTH_RANGE = (24, 25)