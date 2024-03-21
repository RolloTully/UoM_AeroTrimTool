from sympy.solvers import solve
from sympy import symbols, Eq, cos, sin, nsolve
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
class Atmosphere():
    def __init__(self):
        '''Constants'''
        self.R = 287.05
        self.Lapse_rate = -0.0065
        self.T_0 = 288.16
        self.rho_0 = 1.225
        self.g = 9.81

    def Atmosphere_state(self, h):
        '''Valid up to 36000 ft'''
        self.T = self.T_0 + self.Lapse_rate*h
        self.rho = 1.225*(self.T/288.16)**(-((self.g/(self.Lapse_rate*self.R))+1))
        self.sigma = self.rho/1.225
        return self.T, self.rho, self.sigma


class Aircraft():
    def __init__(self):
        '''Constants'''
        self.g = 9.81
        self.atmosphere = Atmosphere()
        '''Flight Conditions'''
        self.alt = 1829
        self.T, self.rho, self.sigma = self.atmosphere.Atmosphere_state(self.alt)
        print(self.rho)
        input()

        self.gamma = 0 #flight path angle

        '''Aircraft properties'''
        #Main Wing
        self.mass = 8.2
        #self.cg_pos = 0.29
        self.Wing_area = 0.633
        self.Span = 2.7
        self.MAC = 0.24
        self.Sweep = (2.65/180)*np.pi
        self.Wing_theta_0 = -(2/180)*np.pi
        self.Wing_z_position = -0.263
        self.Wing_AR = self.Span/self.MAC
        #Tail
        self.Tail_MAC = 0.133
        self.Tail_span = 0.3125*2
        self.Tail_area = self.Tail_span*self.Tail_MAC
        self.Tail_leaver_arm = 1.055
        self.Tail_theta_0 = (-0/180)*np.pi
        self.Tail_offset = 1.055 # quater chord to quater chord
        self.tail_z_position = -0.431
        #Body
        self.Diamiter = 0.26
        #Engines
        self.Thrust_line_z = -0.25
        self.Thrust_line_angle = (0/180)*np.pi
        #Wingbody Aero
        self.Wing_a = 2*np.pi*(self.Wing_AR/(2+self.Wing_AR))
        self.Cl_max = 1.27
        self.Cm_0 = -0.075
        self.Cd_0 = 0.02
        self.alpha_0 = (-2/180)*np.pi
        self.Aero_centre = -0.08
        self.Wing_setting_angle = (2.291/180)*np.pi
        #Tail Aero
        self.Tail_aspect_ratio = (self.Tail_span**2)/self.Tail_area
        self.Tail_a = ((2*np.pi*self.Tail_aspect_ratio)/(2+self.Tail_aspect_ratio))
        self.Elevator_a = 0.6*self.Tail_a
        self.W_0 = 2.0
        '''Derived Values'''
        #Wing

        self.Semi_span = lambda b : b/2
    '''General Functions'''
    def AR(self):
        return (self.Span**2)/self.Wing_area
    '''Tail Functions'''
    def Tail_moment_arm(self,cg):
        '''Cg to tail quater chord'''
        return self.Tail_offset - self.MAC*(cg)
    def Tail_volume(self,cg):
        return (self.Tail_area*self.Tail_leaver_arm/(self.Wing_area*self.MAC))

    def Tail_x_position(self):
        return self.Tail_leaver_arm/self.Span
    def Tail_z_position(self):
        return (self.Wing_z_position-self.tail_z_position)/self.Span
    '''Drag Factors'''
    def Fuselage_drag_factor(self):
        c = self.Diamiter/self.Span
        return 0.9998+(0.0421*c)-(2.6286*c**2)+(2*c**3)
    def Empirical_constant(self):
        return -3.333e-4*self.Sweep**2 + 6.667e-5*self.Sweep + 0.38
    def Oswald_efficiency_factor(self):
        return 1/(np.pi*self.AR()*self.Empirical_constant()*self.Cd_0+(1/(0.99*self.Fuselage_drag_factor())))
    def Induced_drag_factor(self):
        return 1/(np.pi*self.AR()*self.Oswald_efficiency_factor())

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
        self.Temperature, self.Density, self.Density_ratio = self.atmosphere.Atmosphere_state(self.alt)
        print(np.sqrt((2*self.mass*self.g)/(self.Density*self.Wing_area*self.Cl_max)))
        return np.sqrt((2*self.mass*self.g)/(self.Density*self.Wing_area*self.Cl_max))

    def Tail_Downwash(self):
        self.x = self.Tail_offset/self.Span
        self.z = (self.Wing_z_position-self.Tail_z_position())/self.Span
        self.d1 = self.Wing_a/(self.Wing_AR*np.pi**2)
        self.d2 = lambda f: (0.5*np.cos((f*np.pi)/180))**2
        self.d3 = lambda f: np.sqrt(self.x**2+(0.5*np.cos((f*np.pi)/180))**2+self.z**2)
        self.d4 = lambda f: (self.x+np.sqrt((self.x**2)+(0.5*np.cos((f*np.pi)/180))**2+self.z**2))
        self.d5 = lambda f: (0.5*np.cos((f*np.pi)/180))**2+self.z**2
        self.d6 = self.x/((self.x**2)+(self.z**2))
        return self.d1*np.sum(np.array([(self.d2(f)/self.d3(f))*((self.d4(f)/self.d5(f))+self.d6)*(np.pi/180) for f in range(5,175)]))

    def Trim(self,  cg_pos, FPA, V_inf):
        self.Temperature, self.Density, self.Density_ratio = self.atmosphere.Atmosphere_state(self.alt)
        self.alpha_e, self.C_tau, self.Cd, self.C_lt, self.C_lw, self.Cl = symbols('ae, Ctau, Cd, Clt Clw Cl')

        self.eq1 = Eq(2*((self.mass*self.g)/(self.Density*self.Wing_area*V_inf**2))*sin(self.alpha_e+FPA), self.C_tau*cos(0)-self.Cd*cos(self.alpha_e)+self.Cl*sin(self.alpha_e))
        self.eq2 = Eq(2*((self.mass*self.g)/(self.Density*self.Wing_area*V_inf**2))*cos(self.alpha_e+FPA), self.C_tau*sin(0)+self.Cd*sin(self.alpha_e)+self.Cl*cos(self.alpha_e))
        self.eq5 = Eq(0, (self.Cm_0+(cg_pos-self.Aero_centre)*self.C_lw)-self.Tail_volume(cg_pos)*self.C_lt+self.C_tau*(self.Tail_z_position()/self.MAC))
        self.eq3 = Eq(self.Cd, self.Cd_0+self.Induced_drag_factor()*self.Cl**2)
        self.eq4 = Eq(self.Cl, self.C_lw + self.C_lt*(self.Tail_area/self.Wing_area))

        self.eq6 = Eq(self.C_lw, self.Wing_a*(self.alpha_e+self.Wing_setting_angle-self.Wing_theta_0))

        self.results = nsolve((self.eq1,self.eq2,self.eq3,self.eq4,self.eq5,self.eq6),(self.alpha_e, self.C_tau, self.Cd, self.C_lt, self.C_lw, self.Cl),(0.7,0.5,0.02,0.4,0.1,0.1))
        print("V infinity: ",V_inf, " Alpha e: ", str(180*(self.results[0]/np.pi))[0:8], " C tau: ", str(self.results[1])[0:8], " Cd: ", str(self.results[2])[0:8], "C lt: ", str(self.results[3])[0:8], " C lw: ", str(self.results[4])[0:8], " Cl: ", str(self.results[5])[0:8])
        return self.results
    def V_max_range(self):
        self.CL_Max_SAR = np.sqrt(self.Cd_0/(3*self.Induced_drag_factor()))
        self.V_Max_SAR = np.sqrt((2*self.mass*self.g)/(self.Wing_area*self.CL_Max_SAR*self.rho))
        return self.CL_Max_SAR, self.V_Max_SAR
class Main():
    def __init__(self):
        self.aircraft = Aircraft()
        self.mainloop()
    def mainloop(self):
        '''Performance Parameters'''

        self.CL_Max_SAR, self.V_Max_SAR = self.aircraft.V_max_range()
        print("V Max SAR", (self.V_Max_SAR/0.515)*np.sqrt(self.aircraft.sigma))
        print(self.aircraft.V_Stall())
        input()
        '''Arrays to store data'''
        self.V_array_composite = []
        self.alpha_e_array_composite = []
        self.c_tau_array_composite = []
        self.Cd_array_composite = []
        self.C_lt_array_composite = []
        self.C_lw_array_composite = []
        self.Cl_array_composite  = []
        self.Cg_array_composite = []
        self.cg_array = []
        self.other_array = []
        '''Defining colours so it looks pretty'''
        self.colours=['lightcoral','darkorange','darkgreen','navy','mediumvioletred','darkturquoise','teal','darkolivegreen','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal','teal']

        for i,kn in tqdm(enumerate([x/10 for x in range(0,10,1)])):
            print("Static margin: ",kn)
            self.hn = self.aircraft.Aero_centre+self.aircraft.Tail_volume(0)*(self.aircraft.Tail_a/self.aircraft.Wing_a)*(1-self.aircraft.Tail_Downwash())
            self.h = -(kn-self.hn)
            print("=============================================================================")
            '''Parsing over different cg positions'''
            self.V_array = []
            self.alpha_e_array = []
            self.c_tau_array = []
            self.Cd_array = []
            self.C_lt_array = []
            self.C_lw_array = []
            self.Cl_array  = []
            self.cg_position = self.aircraft.Aero_centre + self.aircraft.MAC*kn
            print("Cg position: ",self.cg_position, end=' ')
            for v in range(int(self.aircraft.V_Stall()*0.9),38,1):
                '''Parses over different value of true airspeed'''
                self.results = self.aircraft.Trim(self.h,0, v)
                self.V_array.append(v)
                self.alpha_e_array.append(self.results[0])
                self.c_tau_array.append(self.results[1])
                self.Cd_array.append(self.results[2])
                self.C_lt_array.append(self.results[3])
                self.C_lw_array.append(self.results[4])
                self.Cl_array.append(self.results[5])
            print(self.aircraft.Tail_volume(self.cg_position))
            print(self.aircraft.Elevator_a)
            print(self.cg_position)
            #input()
            self.other_array.append(-self.cg_position/(self.aircraft.Tail_volume(self.cg_position)*self.aircraft.Elevator_a))
            self.V_array_composite.append(self.V_array)
            self.alpha_e_array_composite.append(self.alpha_e_array)
            self.c_tau_array_composite.append(self.c_tau_array)
            self.Cd_array_composite.append(self.Cd_array)
            self.C_lt_array_composite.append(self.C_lt_array)
            self.C_lw_array_composite.append(self.C_lw_array)
            self.Cl_array_composite .append(self.Cl_array)
            self.cg_array.append(self.cg_position)

        '''we now have a bunch of arrays of array of data, yes its horrible and messy but i also hate matplotlib so yk 🤷'''
        self.V_array_composite = np.array(self.V_array_composite)
        self.alpha_e_array_composite = np.array(self.alpha_e_array_composite)
        self.c_tau_array_composite = np.array(self.c_tau_array_composite)
        self.Cd_array_composite = np.array(self.Cd_array_composite)
        self.C_lt_array_composite = np.array(self.C_lt_array_composite)
        self.C_lw_array_composite = np.array(self.C_lw_array_composite)
        self.Cl_array_composite  = np.array(self.Cl_array_composite)
        '''Setting up like a billion figures'''

        '''Lift to Drag Plot'''
            self.LD_figure = plt.figure()
            self.ax1 = self.LD_figure.add_subplot(111)
            self.ax1.set_title("L/D Ratio")
            self.ax1.set_xlabel("Cruise Velocity knts")
            self.ax1.set_ylabel("L/D Ratio")
            self.ax1.vlines(self.V_Max_SAR/0.515,7,20, label = "V Max SAR", color = 'Blue')#Op Point
            self.ax1.vlines(self.aircraft.V_Stall()/0.515,7,20, label = "Stall Speed", color = "Red")
            self.ax1.vlines((self.aircraft.V_Stall()/0.515)*1.2,7,20, label = "Minimum Operating Speed", color = "Blue", linestyles = 'dashed')
            self.ax1.vlines(69.6544,7,20, label = 'V N0', color = "Red",linestyles='dashed')
            for i in range(0,len(self.V_array_composite)):
                #print(i, self.V_array_composite[i/1.94384 self.C_lw_array[i]/self.Cd_array_composite[i])
                self.line, = self.ax1.plot(self.V_array_composite[i]/0.515, self.C_lw_array_composite[i]/self.Cd_array_composite[i], c = self.colours[i])
                self.line.set_label("Static Margin: "+str(i*10)+'%')
            self.ax1.legend()
            plt.show()

            '''Drag Polar Plot'''
            self.Drag_polar_plot = plt.figure()
            self.ax2 = self.Drag_polar_plot.add_subplot(111)
            self.ax2.set_title("Drag Polar")
            self.ax2.set_xlabel("Drag Coefficient (Cd)")
            self.ax2.set_ylabel("Lift Coefficient (Cl)")
            for i in range(0,len(self.V_array_composite)):
                #print(i, self.V_array_composite[i]/1.94384, self.C_lw_array[i]/self.Cd_array_composite[i])
                self.line, = self.ax2.plot(self.Cd_array_composite[i], self.Cl_array_composite[i], c = self.colours[i])
                self.line.set_label("Static Margin = "+ str(i/10))
            self.ax2.legend()
            plt.show()

            '''Total drag plot'''
            self.Drag_plot = plt.figure()
            self.ax3 = self.Drag_plot.add_subplot(111)
            self.ax3.set_title("Total Drag")
            self.ax3.set_xlabel("Cruise Velocity (Knts)")
            self.ax3.set_ylabel("Total Drag (N)")
            for i in range(0,len(self.V_array_composite)):
                #print(i, self.V_array_composite[i], self.C_lw_array[i]/self.Cd_array_composite[i])
                self.line, = self.ax3.plot(self.V_array_composite[i]*1.94384, 0.5*self.aircraft.Density*self.aircraft.Wing_area*self.Cd_array_composite[i]*self.V_array_composite[i]**2, c = self.colours[i])
                self.line.set_label("Static Margin = "+ str(i/10))
            self.ax3.legend()
            plt.show()



        '''Elevator angle plot'''
        self.Elevator_angle = plt.figure()
        self.ax4 = self.Elevator_angle.add_subplot(111)
        self.ax4.set_title("Trim Elevator angle")
        self.ax4.set_xlabel("Cruise Velocity knts")
        self.ax4.set_ylabel("Elevator Deflection")
        self.ax4.vlines(self.V_Max_SAR/0.515,-15,10, label = "V Max SAR", color = 'Blue')#Op Point
        for i in range(0,len(self.V_array_composite)):
            self.Trim_elevator = (self.C_lt_array_composite[i]/self.aircraft.Elevator_a)-(self.aircraft.Tail_a/self.aircraft.Elevator_a)*((self.alpha_e_array_composite[i]+self.aircraft.Wing_setting_angle)*(1-self.aircraft.Tail_Downwash())+0-self.aircraft.Wing_setting_angle-((2*np.pi)/180) )
            self.line, = self.ax4.plot(self.V_array_composite[i]*1.94384,(self.Trim_elevator/np.pi)*180, c = self.colours[i])
            self.line.set_label("Static Margin = "+ str(i/10))
        self.ax4.legend()
        plt.show()




if __name__ == "__main__":
    Main()
