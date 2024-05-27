import numpy as np
import matplotlib.pyplot as plt
from Structures_3 import Steel

class scaling_tool:
    def __init__(self, P, P_ref, t_ref, d_ref, d, z_ref):
        self.S = np.sqrt(P/P_ref)
        self.t_ref = t_ref
        self.d_ref = d_ref
        self.d = d
        self.z_ref = z_ref

    def t(self, DDs=1):
        return self.t_ref * self.S * DDs

    def z(self):
        return self.z_ref*self.S + (self.d-self.d_ref)

    def mass(self, mat):
        rho = mat.rho
        t = self.t()
        z = self.z()



S = scaling_tool(P_ref=15, P = 30, t_ref=55.341/1000, d_ref=30, d=60, z_ref = 15+30+45,  )
print(S.t(), S.z()-60-15)
print(82/85)
