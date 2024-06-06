'''
========== CD calculator ============

6 cylinders of length: 35.88m, diameter: 1m
2 cylinders length: 28.32m, diameter: 1m
2 cylinders length: 26.7m, diameter: 1m
2 cylinders projected length: 44m, diameter: 1m
'''
import numpy as np
class Drag():
    def __init__(self, V, rho, D_truss):

        self.L1 = 24.34
        self.n1 = 324
        self.L2 = 42.12
        self.n2 = 33
        self.L_vert = 30.38
        self.n_vert = 12*17
        self.L_hor1 = 23.31 
        self.n_hor1 = 16
        self.L_hor2 = 6
        self.n_hor2 = 13*6
        self.L_hor3 = 27.38 
        self.n_hor3 = 4
        self.L_diag1 = 19.146
        self.n_diag1 = 12*4*8
        self.L_diag2 = 20.449 
        self.n_diag2 = 12*4*2



        self.L_side1 = 50
        self.L_side2 = 30.38
        self.L_side3 = np.sqrt(self.L_side1**2+self.L_side2**2)
        self.D = D_truss
        self.A1 = self.L1/self.D
        self.A2 = self.L2/self.D
        self.V = V
        self.rho = rho
     

        self.CD_inf = 0.5
        self.kappa1 = 0.8 #long
        self.kappa2 = 0.75 #medium

        self.kappa3 = 0.65 #short
  


    def compute_Reynolds(self, mu):

        return self.D * self.V * self.rho/mu
    
    def compute_CD_front(self):

        CD1 = self.kappa1 * self.CD_inf
        CD2 = self.kappa2 * self.CD_inf
        CD3 = self.kappa3 * self.CD_inf
        

        D = 0.5 * self.rho* (self.V**2)* self.D * ((self.n_vert*self.L_vert+self.n_hor3*self.L_hor3)*CD1+
                                                   (self.n_hor1*self.L_hor1+self.L_diag1*self.n_diag1+self.L_diag2*self.n_diag2)*CD2+
                                                   (self.n_hor2*self.L_hor2)*CD3)
        D_c = D/self.D
        return D, D_c
    
    def compute_D_side(self):
        D_d = 0.5 * self.rho* (self.V**2)*  (0.87*self.CD_inf*(self.L_side3*12+self.L_side1*13)+self.kappa1*self.CD_inf*24*self.L_vert)
        return D_d

    def placeholder(self, l, d, type: str):
        if type == 'front':
            _, D_c = self.compute_CD_front()
            D = D_c/(self.n_diag1+self.n_diag2+self.n_hor1+self.n_hor2+self.n_hor3+self.n_vert)/np.average(self.L_diag1+self.L_diag2+self.L_hor1+self.L_hor2+self.L_hor3+self.L_vert)*l*d

        else:
            D_c = self.compute_CD_side()
            D = D_c/(2*12+13+12)/np.average(self.L_side1+self.L_side2+self.L_side3)*l*d

        return D



if __name__ == "__main__":
    drag = Drag(10.59, 1.225, 1)
    mu = 1.8e-5
    rho = 1.225
    Re = drag.compute_Reynolds(mu)
    d = drag.placeholder(1,1, type = 'front')
    print(d)

    D_grid, _ = drag.compute_CD_front()
    print(D_grid)
    print(drag.compute_D_side())

    print(Re)