# file: therdmodynamics.py
# primers = [primer, [%GC, Tm, ind, primer_length]]
import math

from config.finder import *
from lamp.struct import Primer, PrimersSet


class Temperature:
    """ 
    Class for working with temperature in LAMP
    """
    def __init__(self):
        # Get ranges for params from config
        self._TM_RANGE = TM_RANGE


    def calc_tm(self, primer_length: int, gc_count: int) -> float:
        """
        Calculating Primer Melting Temperature

        :param primer: primer  
        
        :return: Tm = 81.5 + 16.6(log10([Na+])) + 0.41(%GC) - 575 / primer length
        """
        return round(81.5 + 16.6 * math.log10(NA_PLUS) + 0.41 * gc_count - 575 / primer_length, 2)


    @staticmethod
    def check_temperature_diffrence(primers: PrimersSet, cur_primer: Primer) -> bool:
        """
        Checking that all temperatures from the set doesnt differ more than const value

        :param primers: set of primers
        :param cur_primer: primer

        :return: True, if diffrenece < 1, else False
        """
        for primer in primers:
            if abs(primer.Tm - cur_primer.Tm) > TM_MIN_VALUE:
                return False

        return True


    @staticmethod
    def check_temperature_diffrence_probe(primers: list, probe: list) -> bool:
        """
        Checking that all temperatures from the set doesnt differ more than const value

        :param primers: set of primers
        :param cur_primer: primer

        :return: True, if diffrenece < 1, else False
        """
        for primer in primers:
            if abs(primer.Tm - probe.Tm) < TM_PROBE_MAX_VALUE:
                return False

        return True


    def sort_primer_sets(self, primer_sets: list[list]) -> list[list]:
        """ 
        Sort primer sets by average temperature

        :param primer_sets: designed sets of primers

        :return: sorted sets of primers
        """
        # Calc average temperature
        for ind in range(len(primer_sets)):
            primer_sets[ind].append(sum(
                [primer_sets[ind][k][1][1] for k in range(len(primer_sets[ind]))]
            ) / len(primer_sets[ind]))

        # Sort by average temperature
        primer_sets.sort(key=lambda x: x[-1])

        return [primer_sets[i][:-1] for i in range(len(primer_sets))]
            

class GC:
    """
    Class for calculating %GC
    """
    def __init__(self):
        # Get ranges for params from config
        self._GC_RANGE = GC_RANGE
        
        # Vars for dynamic calculating %GC
        self._gc_count = 0
        self._current_gc_count = 0

    
    def _get_first_primer_gc_count(self, primer: str):
        """
        Get first GC count for dynamic calc %GC

        :param primer: primer
        """
        self._gc_count = primer.count('G') + primer.count('C')
        
        self._current_gc_count = self._gc_count
        

    def _calc_gc(self, nucl: str, primer_length: int) -> float:
        """        
        Dynamic calc GC

        :param nucl: last nucleotide of primer
        :param primer_length: length of primer

        :return: primer %GC
        """
        if nucl in ['G', 'C']:
            self._current_gc_count += 1 
        
        return round(self._current_gc_count / primer_length * 100, 2)
        

    def _check_gc(self, fnucl: str, snucl: str) -> None:
        """        
        Check nucleotide for calc GC

        :param fnucl: deleted nucl
        :param snucl: adding nucl
        :param ind: index of el
        :param dec: if True, GC--, else GC++
        """
        if fnucl in ['G', 'C']:
            self._gc_count -= 1
        
        if snucl in ['G', 'C']:
            self._gc_count += 1

        self._current_gc_count = self._gc_count


class Thermodynamics(Temperature, GC):

    def __init__(self):
        Temperature.__init__(self)
        GC.__init__(self)
    

    def _check_and_get_primer_params(self, nucl: str, primer_length: int) -> list[float, float]:
        """        
        Checking primer for parameters: GC_count, Tm
        (private, because we have dynamic calc GC)

        :param primer: primer

        :return: params: [Gc, Tm], if primer corresponds to the parameters, else False
        """
        # Calculation current primer %GC
        primer_GC = self._calc_gc(nucl, primer_length)

        # Check, if primer temperature is in config %GC range
        if primer_GC >= self._GC_RANGE.left and primer_GC <= self._GC_RANGE.right:

            # Calculation current primer temperature
            primer_Tm = self.calc_tm(primer_length, primer_GC)  

            # Check, if primer temperature is in config Tm range
            if primer_Tm >= self._TM_RANGE.left and primer_Tm <= self._TM_RANGE.right:
                return [primer_GC, primer_Tm]

        return False