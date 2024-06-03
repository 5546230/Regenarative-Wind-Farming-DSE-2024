'''
========== CD calculator ============

6 cylinders of length: 35.88m, diameter: 1m
2 cylinders length: 28.32m, diameter: 1m
2 cylinders length: 26.7m, diameter: 1m
2 cylinders projected length: 44m, diameter: 1m
'''

class Drag():
    def __init__(self, V, rho, D_truss):

        self.L1 = 24.34
        self.n1 = 324
        self.L2 = 42.12
        self.n2 = 33
        self.D = D_truss
        self.A1 = self.L1/self.D
        self.A2 = self.L2/self.D
        self.V = V
        self.rho = rho
     

        self.CD_inf = 0.4
        self.kappa1 = 0.78
        self.kappa2 = 0.85
  


    def compute_Reynolds(self, mu):

        return self.D * self.V * self.rho/mu
    
    def compute_CD(self):

        CD1 = self.kappa1 * self.CD_inf
        CD2 = self.kappa2 * self.CD_inf
        

        D = 0.5 * self.rho* (self.V**2)* self.D * (self.n1*self.L1*CD1+self.n2*self.L2*CD2)
        D_c = D/self.D
        return D, D_c

    def placeholder(self, d):
        _, D_c = self.compute_CD()
        D_c = D_c - 0.5 * self.rho* (self.V**2)*self.n2*self.L2*self.kappa2 * self.CD_inf
        D = D_c*d/self.n1
        return D



if __name__ == "__main__":
    drag = Drag(35, 1.225, 1)
    mu = 1.8e-5
    rho = 1.225
    Re = drag.compute_Reynolds(mu)
    d = drag.placeholder(1)
    print(d)

    D_grid, _ = drag.compute_CD()
    print(D_grid)

    print(Re)