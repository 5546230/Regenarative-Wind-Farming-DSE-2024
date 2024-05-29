'''
========== CD calculator ============

6 cylinders of length: 35.88m, diameter: 1m
2 cylinders length: 28.32m, diameter: 1m
2 cylinders length: 26.7m, diameter: 1m
2 cylinders projected length: 44m, diameter: 1m
'''

class Drag():
    def __init__(self):

        self.L1 = 35.88
        self.L2 = 28.32
        self.L3 = 26.7
        self.L4 = 44
        self.D = 1

    def compute_Reynolds(self, V, rho, mu):

        return self.D * V * rho/mu





drag = Drag()
mu = 1.8e-5
rho = 1.225
Re = drag.compute_Reynolds(10, rho, mu)
print(Re)