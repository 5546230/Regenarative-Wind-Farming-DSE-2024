import matplotlib.pyplot as plt
import numpy as np  
import inputs_BEM_powerCurve as inps
import pandas as pd
class Airfoil:
    def __init__(self, csv):
        self.csv = csv
    def import_polar_data(self):
        data = pd.read_csv(self.csv, header=0, names=["alfa", "cl", "cd", "cm"], sep=',')
        return data['alfa'], data['cl'], data['cd']

class BEMsolver:
    def __init__(self, Uinf, radius, NBlades, TSR, RootLocation_R, TipLocation_R, chord_distribution, twist_distribution, polar_alpha_root, polar_alpha_mid, polar_alpha_tip,polar_cl_root, polar_cl_mid, polar_cl_tip, polar_cd_root, polar_cd_mid, polar_cd_tip, root_boundary_R, mid_boundary_R):
        self.Uinf = Uinf
        self.radius = radius
        self.NBlades = NBlades
        self.TSR = TSR
        self.omega = self.TSR*self.Uinf/self.radius
        self.RootLocation_R = RootLocation_R
        self.TipLocation_R = TipLocation_R
        self.chord_distribution = chord_distribution
        self.twist_distribution = twist_distribution
        self.polar_alpha_root = polar_alpha_root
        self.polar_cl_root = polar_cl_root
        self.polar_cd_root = polar_cd_root
        self.polar_alpha_mid = polar_alpha_mid
        self.polar_cl_mid = polar_cl_mid
        self.polar_cd_mid = polar_cd_mid
        self.polar_alpha_tip = polar_alpha_tip
        self.polar_cl_tip = polar_cl_tip
        self.polar_cd_tip = polar_cd_tip
        self.root_boundary_R = root_boundary_R
        self.mid_boundary_R = mid_boundary_R
    def ainduction(self,CT):
        """
        This function calculates the induction factor 'a' as a function of thrust coefficient CT 
        including Glauert's correction
        """
        a = np.zeros(np.shape(CT))
        CT1=1.816   
        CT2=2*np.sqrt(CT1)-CT1
        a[CT>=CT2] = 1 + (CT[CT>=CT2]-CT1)/(4*(np.sqrt(CT1)-1))
        a[CT<CT2] = 0.5-0.5*np.sqrt(1-CT[CT<CT2])
        return a
    def PrandtlTipRootCorrection(self, r_R, axial_induction):
        """
        This function calcualte steh combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
        given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial induction factor
        """       
        
        temp1 = -self.NBlades/2*(self.TipLocation_R-r_R)/r_R*np.sqrt( 1+ ((self.TSR*r_R)**2)/((1-np.clip(axial_induction, 0, 0.999))**2))
        Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
        Ftip[np.isnan(Ftip)] = 0
        temp1 = self.NBlades/2*(self.RootLocation_R-r_R)/r_R*np.sqrt( 1+ ((self.TSR*r_R)**2)/((1-np.clip(axial_induction, 0, 0.999))**2))
        Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
        Froot[np.isnan(Froot)] = 0
        #print(Froot*Ftip)
        return Froot*Ftip, Ftip, Froot
    def loadBladeElement(self, vnorm, vtan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd):
        """
        calculates the load in the blade element
        """
        vmag2 = vnorm**2 + vtan**2
        inflowangle = np.arctan2(vnorm,vtan)
        alpha = twist + inflowangle*180/np.pi
        cl = np.interp(alpha, polar_alpha, polar_cl)
        cd = np.interp(alpha, polar_alpha, polar_cd)
        lift = 0.5*vmag2*cl*chord
        drag = 0.5*vmag2*cd*chord
        fnorm = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
        ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
        gamma = 0.5*np.sqrt(vmag2)*cl*chord
        return fnorm , ftan, gamma
    def solveStreamtube(self, r1_R, r2_R, chord, twist, polar_alpha, polar_cl, polar_cd ):
        """
        solve balance of momentum between blade element load and loading in the streamtube
        input variables:
        Uinf - wind speed at infinity
        r1_R,r2_R - edges of blade element, in fraction of Radius ;
        rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
        Radius is the rotor radius
        Omega -rotational velocity
        NBlades - number of blades in rotor
        """
        Area = np.pi*((r2_R*self.radius)**2-(r1_R*self.radius)**2) #  area streamtube
        r_R = (r1_R+r2_R)/2 # centroide
        # initiatlize variables
        a = 0.0 # axial induction
        aline = 0.0 # tangential induction factor
        
        Niterations = 100
        Erroriterations =0.00001 # error limit for iteration rpocess, in absolute value of induction
        
        for i in range(Niterations):
            # ///////////////////////////////////////////////////////////////////////
            # // this is the block "Calculate velocity and loads at blade element"
            # ///////////////////////////////////////////////////////////////////////
            Urotor = self.Uinf*(1-a) # axial velocity at rotor
            Utan = (1+aline)*self.omega*r_R*self.radius # tangential velocity at rotor
            # calculate loads in blade segment in 2D (N/m)
            fnorm, ftan, gamma = self.loadBladeElement(Urotor, Utan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd)
            load3Daxial =fnorm*self.radius*(r2_R-r1_R)*self.NBlades # 3D force in axial direction
            # load3Dtan =loads[1]*Radius*(r2_R-r1_R)*NBlades # 3D force in azimuthal/tangential direction (not used here)
        
            # ///////////////////////////////////////////////////////////////////////
            # //the block "Calculate velocity and loads at blade element" is done
            # ///////////////////////////////////////////////////////////////////////

            # ///////////////////////////////////////////////////////////////////////
            # // this is the block "Calculate new estimate of axial and azimuthal induction"
            # ///////////////////////////////////////////////////////////////////////
            # // calculate thrust coefficient at the streamtube 
            CT = load3Daxial/(0.5*Area*self.Uinf**2)
            
            # calculate new axial induction, accounting for Glauert's correction
            anew =  self.ainduction(CT)
            
            # correct new axial induction with Prandtl's correction
            Prandtl, Prandtltip, Prandtlroot = self.PrandtlTipRootCorrection(r_R, anew)
            if (Prandtl < 0.0001): 
                Prandtl = 0.0001 # avoid divide by zero
            anew = anew/Prandtl # correct estimate of axial induction
            a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction

            # calculate aximuthal induction
            aline = ftan*self.NBlades/(2*np.pi*self.Uinf*(1-a)*self.omega*2*(r_R*self.radius)**2)
            aline =aline/Prandtl # correct estimate of azimuthal induction with Prandtl's correction
            # ///////////////////////////////////////////////////////////////////////////
            # // end of the block "Calculate new estimate of axial and azimuthal induction"
            # ///////////////////////////////////////////////////////////////////////
            
            #// test convergence of solution, by checking convergence of axial induction
            if (np.abs(a-anew) < Erroriterations): 
                # print("iterations")
                # print(i)
                break

        return [a , aline, r_R, fnorm , ftan, gamma]
    
    def solveBEM(self):
        delta_r_R = .01
        r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
        results =np.zeros([len(r_R)-1,6]) 

        for i in range(len(r_R)-1):
            r_avg = (r_R[i] + r_R[i+1]) / 2

            chord = np.interp((r_R[i]+r_R[i+1])/2, r_R, self.chord_distribution)
            twist = np.interp((r_R[i]+r_R[i+1])/2, r_R, self.twist_distribution)
            if r_avg <= self.root_boundary_R:
                polar_alpha = self.polar_alpha_root
                polar_cl = self.polar_cl_root
                polar_cd = self.polar_cd_root
            elif self.root_boundary_R < r_avg <= self.mid_boundary_R:
                polar_alpha = self.polar_alpha_mid
                polar_cl = self.polar_cl_mid
                polar_cd = self.polar_cd_mid
            else:
                polar_alpha = self.polar_alpha_tip
                polar_cl = self.polar_cl_tip
                polar_cd = self.polar_cd_tip
            
            results[i,:] = self.solveStreamtube(r_R[i], r_R[i+1], chord, twist, polar_alpha, polar_cl, polar_cd )
            return results
        

class Blade:
    def __init__(self, delta_r_R ,pitch, tip_chord, root_chord,root_twist):
        self.delta_r_R = delta_r_R
        self.pitch = pitch
        self.tip_chord = tip_chord
        self.root_chord = root_chord
        self.root_twist = root_twist

    def getChordTwist(self):
        r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
        chord_distribution = self.root_chord*(1-r_R)+self.tip_chord # meters
        twist_distribution = self.root_twist*(1-r_R)+self.pitch # degrees
        return chord_distribution, twist_distribution

#airfoil characteristics  
root_airfoil = Airfoil('s818.csv')
mid_airfoil = Airfoil('s816.csv')
tip_airfoil = Airfoil('s817.csv')
polar_alpha_root, polar_cl_root, polar_cd_root = root_airfoil.import_polar_data()
polar_alpha_mid, polar_cl_mid, polar_cd_mid = mid_airfoil.import_polar_data()
polar_alpha_tip, polar_cl_tip, polar_cd_tip = tip_airfoil.import_polar_data()


#blade geometry
delta_r_R = .01
pitch = inps.pitch # degrees
tip_chord = inps.tip_chord
root_chord  = inps.root_chord
root_twist = inps.root_twist
root_boundary_R = 0.4
mid_boundary_R = 0.75
blade = Blade(delta_r_R, pitch, tip_chord, root_chord, root_twist)
chord_distribution, twist_distribution = blade.getChordTwist()


#initial conditions
Uinf = inps.V_RATED # unperturbed wind speed in m/s
TSR = inps.TSR # tip speed ratio
Radius = inps.init_Radius
Omega = Uinf*TSR/Radius
NBlades = 3

TipLocation_R =  1
RootLocation_R =  0.2

BEM = BEMsolver(Uinf, Radius, NBlades, TSR, RootLocation_R, TipLocation_R, chord_distribution, twist_distribution, polar_alpha_root, polar_alpha_mid, polar_alpha_tip,polar_cl_root, polar_cl_mid, polar_cl_tip, polar_cd_root, polar_cd_mid, polar_cd_tip, root_boundary_R, mid_boundary_R)
results = BEM.solveBEM()

r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
dr = (r_R[1:]-r_R[:-1])*Radius
CT = np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2))
CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))

print(f'{CT=},{CP=}')