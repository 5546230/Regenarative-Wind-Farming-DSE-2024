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
        self.D = 0.5
        self.A1 = self.L1/self.D
        self.A2 = self.L2/self.D
        self.A3 = self.L3/self.D
        self.A4 = self.L4/self.D   

        self.CD_inf = 0.4
        self.kappa1 = 0.9
        self.kappa2 = 0.88
        self.kappa3 = 0.87
        self.kappa4 = 0.93     



    def compute_Reynolds(self, V, rho, mu):

        return self.D * V * rho/mu
    
    def compute_CD(self, V, rho):

        CD1 = self.kappa1 * self.CD_inf
        CD2 = self.kappa2 * self.CD_inf
        CD3 = self.kappa3 * self.CD_inf
        CD4 = self.kappa4 * self.CD_inf

        D = 0.5 * rho* (V**2)* (CD1*self.L1*self.D + CD2*self.L2*self.D + CD3 * self.L3*self.D + CD4 * self.L4*self.D)
        return D








drag = Drag()
mu = 1.8e-5
rho = 1.225
Re = drag.compute_Reynolds(10, rho, mu)

D_element = drag.compute_CD(60, 1.225)
print(D_element)

print(Re)