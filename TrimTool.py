class Atmosphere(self):
    def __init__(self):
        '''Constants'''
        self.R = 287.05
        self.Lapse_rate = -0.0065
        self.T_0 = 288.16
        self.rho_0 = 1.225

    def Atmosphere_state(self, h):
        '''Valid up to 36000 ft'''
        self.T = self.T_0 + self.Lapse_rate*h
        self.rho = 1.225*(self.T/288.16)**(-((self.g/(self.Lapse_rate*self.R))+1))
        self.sigma = self.rho/1.225
        return self.T, self.rho, self.sigma


class Aircraft(self):
    def __init__(self):
        '''Constants'''
        selg.g = 9.81
        self.atmosphere = Atmosphere()
        '''Flight Conditions'''
        self.alt = 6562
        self.gamma = 0 #flight path angle

        '''Aircraft properties'''
        #Main Wing
        self.mass = 6300
        self.cg_pos = 0.29
        self.Wing_area = 25.08
        self.Span = 15.85
        self.MAC = 1.716
        self.Sweep = 0
        self.Wing_theta_0 = 0
        self.Wing_z_position = 0.45
        #Tail
        self.Tail_area = 7.79
        self.Tail_span = 6.6
        self.Tail_leaver_arm = 6.184
        self.Tail_theta_0 = 1.5
        self.Tail_offset = 6.184 # quater chord to quater chord
        self.Tail_z_position = -1.435
        #Body
        self.Diamiter = 1.981
        #Engines
        self.Thrust_line_z = 0.312
        self.Thrust_line_angle = 0
        #Wingbody Aero
        self.Wing_a = 5.19
        self.Cl_max = 1.37
        self.Cm_0 = -0.0711
        self.Cd_0 = 0.03
        self.alpha_0 = -2
        self.Aero_centre = -0.08
        self.Wing_setting_angle = 1
        #Tail Aero
        self.Tail_a = 3.2
        self.Elevator_a = 2.414
        self.W_0 = 2.0
        '''Derived Values'''
        #Wing
        self.AR = lambda b,s:(b**2),s
        self.Semi_span = lambda b:b/2
    '''General Functions'''
    def AR(self, s = self.Span, b=self.Wing_area):
        return (s**2)/b
    '''Tail Functions'''
    def Tail_moment_arm(self):
        '''Cg to tail quater chord'''
        return self.Tail_offset - self.MAC*(self.cg_pos-0.25)
    def Tail_volume(self):
        return (self.Tail_area*self.Tail_moment_arm())/(self.Wing_area*self.MAC)

    def Tail_x_position(self):
        return self.Tail_leaver_arm/self.Span
    def Tail_z_position(self):
        return (self.Wing_z_position-self.Tail_z_position)/self.Span

    def Tail_downwash_angle(self):
        '''Calcuates the induced downwash angle of the tail'''
        '''not finished'''
        return (self.Tail_a/(self.AR()*np.pi**2))*np.sum(np.array([for fi in range(5,85)]))

    '''Drag Factors'''
    def Fuselage_drag_factor(self):
        c = self.Diamiter/self.Span
        return 0.9998+(0.0421*c)-(2.6286*c**2)+(2*c**3)
    def Empirical_constant(self):
        return -3.333e-4*self.Sweep**2 + 6.667e-5*self.Sweep + 0.38
    def Oswald_efficiency_factor(self):
        return 1/(np.pi*self.AR()*self.Cd_0+(1/(0.99self.Fuselage_drag_factor)))
    def Induced_drag_factor(self):
        return 1/(np.pi*self.AR()*self.Oswald_efficiency_factor)

    '''Performace Parameters'''
    def Velocity_min_drag(self, h):
        '''The airspeed that results in the minimum aircraft drag in ms^-1'''
        '''Need to add atmosphere function rho'''
        self.T, self.rho, self.sigma = self.atmosphere.Atmosphere_state(h)
        return (np.sqrt((2*self.mass*self.g)/(self.Wing_area))*(self.Induced_drag_factor()/self.Cd_0)**0.25)
    def Equiv_min_drag(self, h):
        self.T, self.rho, self.sigma = self.atmosphere.Atmosphere_state(h)
        return self.Velocity_min_drag(h)*np.sqrt()

    def V_Stall(self):
        '''Need to add atmosphere function rho'''
        return np.sqrt((2*self.mass*self.g)/())

    def Neutral_point(self):
        pass
    def Cl(self):
    def Cd(self):
    def Cl_wing(self):
    def Cl_T()
    def LD_ratio(self):

    def Trim(self):
        '''Calculates the aircraft trim conditions'''
        '''Finds trim, AoA, Lift, Drag, Thrust, Lift Coefficient, Thrust coefficient, Drag coefficient,
        for V_knts in range(100,250):
            V_ms = V_knts*0.515
            print(V_knts, V_ms)

class Main(self):
    def __init__(self):
        self.mainloop()
    def mainloop(self):


if __name__ == "__main__":
    Main()
