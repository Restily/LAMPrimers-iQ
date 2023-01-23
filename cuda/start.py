import numpy


from pycuda.compiler import SourceModule


class Primer:
    def __init__(self) -> None:
        pass

    def __str__(self):
        return str(cuda.from_device(self.data, self.shape, self.dtype))


with open('cuda/code.cpp') as f:
    mod = SourceModule(f.read())


start_lamp = mod.get_function('start_lamp')

