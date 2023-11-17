# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 17:42:43 2023

@author: gavin
"""
import math
import scipy.interpolate
import copy
import sys

#classes below are for calculating properties given values
class rebar:
    
    def __init__(self,usage):
        #units in mm
        self.usage = usage
        self.size_list = [6,8,10,12,16,20,25,32,40]
        self.size_counter = 0
        self.no = 2 #minimum number of rebars is 2
        #MPa
        self.E = 200*10**3
        #MPa
        self.fyk = 500
        self.safety_factor = 1.15
        #approximated ultimate steel strength
        self.fuk = 600
        
    def increase_no(self):
            self.no = self.no + 1
    
    def decrease_no(self):
            self.no = self.no - 1
        
    @property
    def d(self):
        return self.size_list[self.size_counter]
    
    def increase_d(self):
        self.size_counter = self.size_counter + 1  
    
    def decrease_d(self):
       self.size_counter = self.size_counter - 1  

    @property
    def area(self):
        return 0.25*(self.d**2)*math.pi*self.no
    
    @property
    def fyd(self):
        return self.fyk/self.safety_factor
    
    @property
    def yield_strain(self):
        return self.fyd/self.E
    
    @property
    def fud(self):
        return self.fuk/self.safety_factor
    
    def stress(self,strain):
        if strain >= self.yield_strain:
            return self.fyd
        else:
            return strain*self.E

class concrete:   
    
    def __init__(self,grade):
        #in MPa
        self.fck = grade
        #kN/m^3
        self.density= 25
        self.failure_strain = 0.0035
        self.safety_factor = 1.5
    
    @property
    def fck_cube(self):
        return self.fck/0.8
    
    @property
    def E_cm(self):
        #in GPa
        return 22*((self.fck+8)/10)**0.3
    
    @property
    def fcd(self):
        return self.fck/self.safety_factor

class shear_reinforcement:
    def __init__(self):
        self.needed = True #default to needed
        #MPa
        self.E = 200*10**3
        #MPa
        self.fyk = 500
        self.safety_factor = 1.15
        self.Asl = 0
        
        #possible legs
        self.leg_list = [2,4]
        self.leg_counter = 0
        
        #possible thetas
        self.theta_list = range(22,46)
        self.theta_counter = 0
        
        #possible spacings
        self.spacing_list = [80,90,100,125,150,175,200,225,250,275,300]
        self.spacing_counter = 0
        
        #possible bar diameters
        self.size_list = [0,8,10,12,16]
        self.size_counter = 0
        
    def increase_leg(self):
            self.leg_counter = self.leg_counter + 1
    
    def decrease_leg(self):
            self.leg_counter = self.leg_counter - 1
    
    def increase_theta(self):
            self.theta_counter = self.theta_counter + 1
    
    def decrease_theta(self):
            self.theta_counter = self.theta_counter - 1
            
    def increase_spacing(self):
            self.spacing_counter = self.spacing_counter + 1
    
    def decrease_spacing(self):
            self.spacing_counter = self.spacing_counter - 1
            
    def increase_size(self):
            self.size_counter = self.size_counter + 1
    
    def decrease_no(self):
            self.size_counter = self.no - 1
    
    @property
    def spacing(self):
        return self.spacing_list[self.spacing_counter]
    
    @property
    def size(self):
        return self.size_list[self.size_counter]
    
    @property
    def leg(self):
        return self.leg_list[self.leg_counter]
    
    @property
    def theta(self):
        return self.theta_list[self.theta_counter]
    
    @property
    def Asw_s(self):
        if self.needed == True:
            return (self.leg*math.pi*0.25*self.size**0.5)/self.spacing
        else:
            return 0
    
    @property
    def fywd(self):
        return self.fyk/self.safety_factor
    
    @property
    def fyld(self):
        return self.fyk/self.safety_factor
        
class RC_rectangle:
    def __init__(self,C_grade):
        self.l = 0
        self.b_list = range(10,510,10)
        self.b_counter = 0
        self.d = 0
        self.h_list = range(10,1210,10)
        self.h_counter = 0
        self.dp = 0
        self.wk = 0.3 #if no wk provided, assume 0.3
        self.C_grade = C_grade
        self.t_rebar = rebar("tension")
        self.c_rebar = rebar("compression")
        self.concrete = concrete(C_grade)
        self.shear_reinforcement = shear_reinforcement()
        
        #set these values using class methods below
        self.c_concrete_strain = "Not set"
        
        #calculate these values using class methods below
        self.c_rebar_strain = "Not calculated"
        self.t_rebar_strain = "Not calculated"
        self.c_rebar_stress = "Not calculated"
        self.t_rebar_stress = "Not calculated"
        self.NA = "Not calculated"
        self.MRd = "Not calculated"
        
        #shear related values
        self.VEd = 0 #default to no shear loading required
        
        #reference tables
        self.max_bar_d_table = {0.4:{160:40,200:32,240:20,280:16,320:12,360:10,400:8,450:6},0.3:{160:32,200:25,240:16,280:12,320:10,360:8,400:6,450:5},0.2:{160:25,200:16,240:12,280:8,320:6,360:5,400:4,450:0}}
        
        self.max_bar_spacing_table = {0.4:{160:300,200:300,240:250,280:200,320:150,360:100},0.3:{160:300,200:250,240:200,280:150,320:100,360:50},0.2:{160:200,200:150,240:100,280:50,320:0,360:0}}
        
        
    @property
    def b(self):
        return self.b_list[self.b_counter]
    
    @property
    def h(self):
        return self.h_list[self.h_counter]    
    
    @property
    def total_area(self):
        return self.b*self.h
    
    @property
    def t_rebar_area(self):
        return self.t_rebar.area
    
    @property
    def c_rebar_area(self):
        return self.c_rebar.area
        
    @property
    def concrete_area(self):    
        return self.total_area - self.t_rebar_area - self.c_rebar_area

    def set_values_collapse(self, tolerance):
        #tolerance is in N
        #set the concrete strain, calculate the rest
        #done together because it results come as a set
        self.c_concrete_strain = self.concrete.failure_strain
        
        #always initiate by assuming NA is half of h
        self.NA = self.h/2
        self.c_rebar_strain = "Not calculated"
        self.t_rebar_strain = "Not calculated"
        self.c_rebar_stress = "Not calculated"
        self.t_rebar_stress = "Not calculated"
        #MRd will be in Nmm
        self.MRd = "Not calculated"
        
        
        #guess and check NA
        NA = self.NA
        dp = self.dp
        d = self.d
        top_strain = self.c_concrete_strain
        strain_sc = top_strain*(dp-NA)/NA
        strain_st = top_strain*(d-NA)/NA
        
        #calculate stress resultants
        #abs is used because positive values are used, even in compression
        stress_sc = self.c_rebar.stress(abs(strain_sc))
        stress_st = self.t_rebar.stress(strain_st)
        
        T = self.t_rebar_area*stress_st
        C_s = self.c_rebar_area*stress_sc
        C_c = 0.567*self.concrete.fck*self.b*0.8*NA
        C = C_s + C_c
        
        #Guess and check part
        #if within tolerance, set rebar strains
        NA_LBound = 0
        NA_UBound = self.h
        #NA should not lie within the steel region
        while abs(T-C) > tolerance and (self.c_rebar_area == 0 or NA_UBound > self.dp) and NA_LBound < self.d:
            if abs(T-C) <= tolerance:
            #honestly not necessary to be inside the loop but makes it cleaner
                self.c_rebar_strain = strain_sc
                self.t_rebar_strain = strain_st
            else:
                if C > T:
                    #If C > T, reduce NA
                    NA_UBound = NA
                    NA = (NA + NA_LBound)/2
                elif T > C:
                    #If T > C, increase NA
                    NA_LBound = NA
                    NA = (NA + NA_UBound)/2
            
            
            #recalculate
            strain_sc = top_strain*(dp-NA)/NA
            strain_st = top_strain*(d-NA)/NA
            stress_sc = self.c_rebar.stress(abs(strain_sc))
            stress_st = self.t_rebar.stress(strain_st)
            
            T = self.t_rebar_area*stress_st
            C_s = self.c_rebar_area*stress_sc
            C_c = 0.567*self.concrete.fck*self.b*0.8*NA
            C = C_s + C_c
        
        if NA < self.d and (self.c_rebar_area == 0 or NA > self.dp):
            self.NA = NA
            self.c_rebar_strain = strain_sc
            self.t_rebar_strain = strain_st
            self.c_rebar_stress = stress_sc
            self.t_rebar_stress = stress_st
            #MRd will be in Nmm
            self.MRd = T*(d-0.4*NA) + C_s*(0.4*NA-dp)
        else:
            #setting all to zero = no load carrying capacity
            self.NA = 0
            self.c_rebar_strain = 0
            self.t_rebar_strain = 0
            self.c_rebar_stress = 0
            self.t_rebar_stress = 0
            #MRd will be in Nmm
            self.MRd = 0
    
    ###############################
    #Shear related functions below#
    ###############################
    
    @property
    def k(self):
        return min(1+(200/self.d)**0.5,2)
    
    @property
    def rho1(self):
        return min(self.t_rebar_area/(self.b*self.d),0.02)
    
    @property
    def nu(self):
        return 0.6*(1-self.concrete.fck/250)
    
    @property
    def z(self):
        return 0.9*self.d
    
    @property
    def VRd_c(self):
        value = (0.18/self.concrete.safety_factor)*self.k*self.b*self.d*(100*self.rho1*self.concrete.fck)**(1/3)
        limit = 0.035*(self.k**(3/2))*(self.concrete.fck**0.5)*self.b*self.d
        return max(value,limit)
    
    @property
    def VRd_max(self):
        numerator = self.nu*self.concrete.fcd*self.b*self.z
        denominator = math.tan(self.shear_reinforcement.theta)
        return numerator/denominator
    
    @property
    def Asw_req(self):
        cot = 1/(math.tan(math.radians(self.shear_reinforcement.theta)))
        minimum = (0.08*self.b*self.concrete.fck**0.5)/self.shear_reinforcement.fyk
        req = self.VEd/(self.shear_reinforcement.fywd*self.z*cot)
        return max(minimum,req)
    
    @property
    def Asl_req(self):
        cot = 1/(math.tan(math.radians(self.shear_reinforcement.theta)))
        return (self.VEd*cot)/(2*self.shear_reinforcement.fyld)
    
    @property
    def shear_area(self):
        SR = self.shear_reinforcement
        return self.Asl_req + SR.Asw_s
    
    
    #############################
    #SLS related functions below#
    #############################
    
    #stresses below are all assumed to be 310 (conservatively) Lecture 3, Pg 16
    @property
    def max_bar_d(self):
        wk = self.wk
        stress = 310
        limit_dict = self.max_bar_d_table[wk]
        X = sorted(list(limit_dict.keys()))
        Y = []
        for x in X:
            Y.append(limit_dict[x])
        y_interp = scipy.interpolate.interp1d(X,Y)
        
        #only do the fck multiplier because it is the only reasonably accountable one
        return ((0.3*self.concrete.fck**(2/3))/2.9)*y_interp(stress)
    
    @property 
    def max_bar_spacing(self):
        wk = self.wk
        stress = 310
        limit_dict = self.max_bar_spacing_table[wk]
        X = sorted(list(limit_dict.keys()))
        Y = []
        for x in X:
            Y.append(limit_dict[x])
        y_interp = scipy.interpolate.interp1d(X,Y)
        
        return y_interp(stress)
        return
    

#function to skip all the values smaller than min based on the list properties
def min_counter(lst, val):
    #note that lst is assumed to be sorted
    i = 0
    #edge case where value to be found is greater than final value
    if val > lst[len(lst)-1]:
        return "All values too small"
    
    while val > lst[i]:
        i = i + 1
    return i
        
def min_rebar_single(RC, M, min_spacing,tolerance):
    Mp = 0.167*RC.concrete.fck*RC.b*RC.d**2
    #Kp = Mp/(RC.b*RC.concrete.fck*RC.d**2)
    #K = M/(RC.b*RC.concrete.fck*RC.d**2)
    
    initial_t_no = RC.t_rebar.no
    initial_c_no = RC.c_rebar.no
    initial_t_size_counter = RC.t_rebar.size_counter        
    initial_c_size_counter = RC.c_rebar.size_counter
    
    max_spacing = RC.max_bar_spacing
    max_d = RC.max_bar_d
    
    #[no, d]    
    min_t_rebar = []
    min_c_rebar = []
    min_t_area = 0
    min_c_area = 0    
    
    RC.c_rebar.no = 0
    min_c_rebar = [0,0]
    for i in range(initial_t_size_counter, len(RC.t_rebar.size_list)):
        RC.t_rebar.size_counter = i
        
        #skip the sizes that are larger than max_d
        if RC.t_rebar.d > max_d:
            continue
        
        RC.set_values_collapse(tolerance)
        #assume that distance from side is same as distance from edge
        spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
        
        #in case the initial values are suitable
        if RC.MRd >= M and spacing >= min_spacing and spacing <= max_spacing:
            min_t_rebar = [RC.t_rebar.no, RC.t_rebar.d]
            min_t_area = RC.t_rebar.area
            
        #skip to the rebar number that fulfils max spacing requirements
        while spacing > max_spacing:
            RC.t_rebar.no = RC.t_rebar.no + 1
            spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
        
        while RC.MRd < M and spacing >= min_spacing:
            RC.t_rebar.no = RC.t_rebar.no + 1
            spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
            RC.set_values_collapse(tolerance)
            

        #after looping through all the possible spacings, either
        # 1) Smallest n has been found
        # 2) Spacing limits exceeded
        if RC.MRd >= M and (min_t_area == 0 or RC.t_rebar.area < min_t_area):
            #Use < for area so that the smallest diameters + most bars will always be chosen if tied. This allows better bonding.
            min_t_rebar = [RC.t_rebar.no, RC.t_rebar.d]
            min_t_area = RC.t_rebar.area
        
            
        #Reset rebar count
        RC.t_rebar.no = initial_t_no
    #Reset rebar size
    RC.t_rebar.size_counter = initial_t_size_counter

    return [min_t_rebar,min_c_rebar]


def min_rebar_double(RC,M,min_spacing,tolerance):
    Mp = 0.167*RC.concrete.fck*RC.b*RC.d**2
    #Kp = Mp/(RC.b*RC.concrete.fck*RC.d**2)
    #K = M/(RC.b*RC.concrete.fck*RC.d**2)
    
    initial_t_no = RC.t_rebar.no
    initial_c_no = RC.c_rebar.no
    initial_t_size_counter = RC.t_rebar.size_counter        
    initial_c_size_counter = RC.c_rebar.size_counter
    
    max_spacing = RC.max_bar_spacing
    max_d = RC.max_bar_d
    
    #[no, d]    
    min_t_rebar = []
    min_c_rebar = []
    min_t_area = 0
    min_c_area = 0
    
    
    #If M > Mp compression rebar is needed
    #ensure compression rebar exists
    if RC.c_rebar.no < 2:
        RC.c_rebar.no = 2
        initial_c_no = RC.c_rebar.no
    if RC.c_rebar.area == 0:
        RC.c_rebar.size_counter = RC.c_rebar.size_counter + 1
        initial_c_size_counter = RC.c_rebar.size_counter

        
    for i in range(initial_t_size_counter, len(RC.t_rebar.size_list)):
        RC.t_rebar.size_counter = i
        
        
        #skip the sizes that are larger than max_d
        if RC.t_rebar.d > max_d:
            continue
        
        t_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
        c_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.c_rebar.no-1)
        #set the initial values for the loops
        RC.set_values_collapse(tolerance)
        
        #in case the initial values are suitable
        if RC.MRd >= M and t_spacing >= min_spacing and t_spacing <= max_spacing and c_spacing >= min_spacing:
            min_t_rebar = [RC.t_rebar.no, RC.t_rebar.d]
            min_t_area = RC.t_rebar.area
            min_c_rebar = [RC.c_rebar.no, RC.c_rebar.d]
            min_c_area = RC.c_rebar.area 
        
        #skip to the rebar number that fulfils max spacing requirements
        while t_spacing > max_spacing:
            RC.t_rebar.no = RC.t_rebar.no + 1
            t_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
        
        while RC.MRd < M and t_spacing >= min_spacing:
            RC.t_rebar.no = RC.t_rebar.no + 1
            t_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.t_rebar.no-1)
            
            for j in range(initial_c_size_counter, len(RC.c_rebar.size_list)):
                RC.c_rebar.size_counter = j
                RC.set_values_collapse(tolerance)
                
                min_total_area = min_t_area + min_c_area
                total_area = RC.t_rebar.area + RC.c_rebar.area
                
                #skip the checks if it will not reduce the total area
                if total_area > min_total_area and min_total_area != 0:
                    continue
                
                c_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.c_rebar.no-1)
                
                while RC.MRd < M and c_spacing >= min_spacing:
                    RC.c_rebar.no = RC.c_rebar.no + 1
                    c_spacing = (RC.b - 2*(RC.h-RC.d))/(RC.c_rebar.no-1)
                    RC.set_values_collapse(tolerance)
                          
                #after looping through all the possible spacings, either
                # 1) Smallest n has been found
                # 2) Spacing limits exceeded
                if RC.MRd >= M:
                    min_t_rebar = [RC.t_rebar.no, RC.t_rebar.d]
                    min_t_area = RC.t_rebar.area
                    min_c_rebar = [RC.c_rebar.no, RC.c_rebar.d]
                    min_c_area = RC.c_rebar.area 
                #reset c_no
                RC.c_rebar.no = initial_c_no
            #reset c_size
            RC.c_rebar.size_counter = initial_c_size_counter
        #reset t_no
        RC.t_rebar.no = initial_t_no
    #reset t_size
    RC.t_rebar.size_counter = initial_t_size_counter
    
    return [min_t_rebar,min_c_rebar]

def area_rebar_list(lst):
    #returns area from [[no,d],[no,d]]
    #Placeholder R, tension or compression does not matter
    R = rebar("tension")
    R.no = lst[0][0]
    R.size_counter = min_counter(R.size_list,lst[0][1])
    area1 = R.area
    
    R.no = lst[1][0]
    R.size_counter = min_counter(R.size_list,lst[1][1])
    area2 = R.area
    
    return area1 + area2
    
#find the minimum rebar for flexure
def min_rebar_design(RC, M, min_spacing,tolerance):
    Mp = 0.167*RC.concrete.fck*RC.b*RC.d**2
    #Kp = Mp/(RC.b*RC.concrete.fck*RC.d**2)
    #K = M/(RC.b*RC.concrete.fck*RC.d**2)
    
    if M > Mp:
        result = min_rebar_double(RC, M, min_spacing,tolerance)
        return result           
    else:
        result1 = min_rebar_single(RC, M, min_spacing, tolerance)
        result2 = min_rebar_double(RC, M, min_spacing,tolerance)

        if len(result1[0]) == 0 and len(result2[0]) == 0:
            #both no possibilities, return empty list
            return result2
        
        if len(result1[0]) == 0:
            return result1
        
        if len(result2[0]) == 0:
            return result2            
            
        if area_rebar_list(result2) < area_rebar_list(result1):
        #if tied, do single reinforcement
            return result2
        else:
            return result1
    
     
    
#check if deflection is sufficiently limited
def deflection_check(RC, K):
    #assume that provided reinforcement = required reinforcement
    #assuming that sigma_s = 310MPa (conservative)
    fck = RC.concrete.fck
    rho0 = 0.001*fck**0.5
    rho = RC.t_rebar_area/(RC.b*RC.d)
    rhop = RC.c_rebar_area/(RC.b*RC.d)
    
    if rho <= rho0:
        ld_limit = K*(11+1.5*(fck**0.5)*(rho0/rho)+3.2*(fck**0.5)*(rho0/rho)**1.5)
    elif rho > rho0:
        ld_limit = K*(11+1.5*(fck**0.5)*(rho0/(rho-rhop))+(1/12)*(fck**0.5)*(rhop/rho0)**0.5)
        
    if RC.l > 7:
        ld_limit = ld_limit * 7/RC.l
        
    if RC.l/RC.d > ld_limit:
        return False
    else:
        return True     

def yield_check(RC):
    if RC.t_rebar_stress == RC.t_rebar.fyd:
        return True
    else:
        return False
    
#############
#shear stuff#
#############
def VRd_max_sufficient(RC):
    #first, set theta to maximum (assuming 22deg to 45 deg, it strictly increases VRd_max)
    initial_theta_counter = RC.shear_reinforcement.theta_counter
    RC.shear_reinforcement.theta_counter = len(RC.shear_reinforcement.theta_list) - 1
    
    if RC.VEd > RC.VRd_max:
        #reset the object
        RC.shear_reinforcement.theta_counter = initial_theta_counter
        return False
    else:
        #reset the object
        RC.shear_reinforcement.theta_counter = initial_theta_counter
        return True
    
def optimise_shear_area(RC):
    #return [leg,size,theta,spacing]
    SR = RC.shear_reinforcement
    initial_theta_counter = SR.theta_counter
    initial_leg_counter = SR.leg_counter
    initial_size_counter = SR.size_counter
    initial_spacing_counter = SR.spacing_counter
    VEd = RC.VEd
    
    optimal = []
    optimal_area = 0
    #iterate through theta in reverse since you should maximise it
    for i in range(len(SR.theta_list) - 1,initial_theta_counter - 1, -1):
        SR.theta_counter = i
        for j in range(initial_leg_counter,len(SR.leg_list)):
            SR.leg_counter = j
            for k in range(initial_spacing_counter, len(SR.spacing_list)):
                SR.spacing_counter = k
                for l in range(initial_size_counter, len(SR.size_list)):
                    SR.size_counter = l
                    VRd_max = RC.VRd_max
                    area = RC.shear_area
                    if VRd_max >= VEd:
                        if area < optimal_area or optimal_area == 0:
                            #use < because prioritise lower legs and spacing
                             optimal = [SR.leg,SR.size,SR.theta,SR.spacing]
    #reset values
    SR.theta_counter = initial_theta_counter 
    SR.leg_counter = initial_leg_counter 
    SR.size_counter = initial_size_counter 
    SR.spacing_counter = initial_spacing_counter 
    
    return optimal
    

################
#da actual lööp#
################

#start of the code which optimises the parameters for class RC_rectangle
#first,set the values (e.g. fixed values, min/max values)

#b_min is determined by fire resistance
b_min = 300

#a_min is the minimum distance of the rebars from the edge. Determined by fire resistance and corrosion
a_min = min(35,30)

#h_min is just to limit the search space a little bit
h_min = 100

#next, optimise for moments
RC = RC_rectangle(40)
b_min_counter = min_counter(RC.b_list,b_min)
RC.b_counter = b_min_counter
h_min_counter = min_counter(RC.h_list, h_min)
RC.h_counter = h_min_counter
RC.wk = 0.3

#min_spacing (in mm)
min_spacing = 15

#deflection check K
deflection_K = 1.5

#tolerance in N
NA_tolerance = 1000

#Make d as large as possible, given the current dimensions
RC.d = RC.h - a_min
RC.dp = a_min


#optimisation aim (what's the goal?)
def reward_value(RC):
    return RC.b * RC.h
    #return 0.9*RC.concrete_area + 1.8*(RC.c_rebar_area + RC.t_rebar_area + RC.shear_area)

#looooop (exhaustive search)
#target_M is in Nmm
target_M = 95*10**6
Flexure_ULS_list = [] #list of RC objects that fulfil flexure_ULS + conservative SLS assumptions
for i in range(b_min_counter, len(RC.b_list)):
    RC.b_counter = i
    for j in range(h_min_counter, len(RC.h_list)):
        RC.h_counter = j
        
        #re-adjust depths
        RC.d = RC.h - a_min
        RC.dp = a_min
        
        #find minimum rebar area required
        rebar_design = min_rebar_design(RC, target_M, min_spacing, NA_tolerance)
        t_rebar_design = rebar_design[0]
        c_rebar_design = rebar_design[1]
        
        #if either one is a blank list, it means that no appropriate rebar design is possible
        if len(t_rebar_design) != 0 and len(c_rebar_design) != 0:
            RC.t_rebar.no = t_rebar_design[0]
            
            if t_rebar_design[1] != 0:
                RC.t_rebar.size_counter = min_counter(RC.t_rebar.size_list, t_rebar_design[1])
            else:
                RC.t_rebar.size_counter = 0 #set to smallest size if rebar not needed anyway
            RC.c_rebar.no = c_rebar_design[0]
            
            if c_rebar_design[1] != 0:
                RC.c_rebar.size_counter = min_counter(RC.c_rebar.size_list, c_rebar_design[1])
            else:
                RC.c_rebar.size_counter = 0 #set to smallest size if rebar not needed anyway
   
            #just to update the values properly
            RC.set_values_collapse(NA_tolerance)
            
            if deflection_check(RC, deflection_K) == True and yield_check(RC) == True:
                #Ok this should technically be inside the rebar design function, so a few edge cases might be missed. But i am running out of time.
                RC_copy = copy.deepcopy(RC)
                Flexure_ULS_list.append(RC_copy)

#print(len(Flexure_ULS_list))
#now, use the Flexure_ULS_list to design for shear
#VEd should be in N
VEd = 72*10**3
ULS_list = []
for RC in Flexure_ULS_list:
    RC.VEd = VEd
    #No need to bother with the RC that can never be designed
    if VRd_max_sufficient(RC):
        optimal_shear = optimise_shear_area(RC)
        if len(optimal_shear) != 0:
            RC_new = copy.deepcopy(RC)
            SR = RC_new.shear_reinforcement
            SR.leg_counter = min_counter(SR.leg_list, optimal_shear[0])
            SR.size_counter = min_counter(SR.size_list, optimal_shear[1])
            SR.theta_counter = min_counter(SR.theta_list, optimal_shear[2])
            SR.spacing_counter = min_counter(SR.spacing_list, optimal_shear[3])
            
            ULS_list.append(RC_new)

#check the reward of all the ULS valid RC
best_reward = 0

#technically the reward value is not a reward because we are minimising it :D
for RC in ULS_list:
    reward = reward_value(RC)
    if best_reward == 0 or reward < best_reward:
        best_reward = reward
        best_RC = RC
        
#print("The best reward is: " + str(best_reward))
print(best_RC.b, best_RC.h)
print(best_RC.t_rebar_strain,best_RC.c_rebar_strain)
print(best_RC.t_rebar_stress,best_RC.c_rebar_stress)
print(best_RC.t_rebar_area, best_RC.c_rebar_area)

print(best_RC.NA, best_RC.MRd)
print(best_RC.t_rebar.no, best_RC.t_rebar.d)
print(best_RC.c_rebar.no, best_RC.c_rebar.d)
print(best_RC.shear_reinforcement.size, best_RC.shear_reinforcement.spacing,best_RC.Asl_req)
   

    # C_c = 0.567*best_RC.concrete.fck*best_RC.b*0.8*best_RC.NA
    # C_s = best_RC.c_rebar_area*best_RC.c_rebar_stress
    # T = best_RC.t_rebar_area*best_RC.t_rebar_stress
    # print(C_c, C_s, C_c + C_s, T)



#test
# RC_test = RC_rectangle(25)
# RC_test.b_counter = min_counter(RC_test.b_list,300)
# RC_test.h_counter = min_counter(RC_test.h_list, 500)
# RC_test.d = 450
# RC_test.dp = 50
# RC_test.t_rebar.no = 3
# RC_test.c_rebar.no = 2
# RC_test.t_rebar.size_counter = min_counter(RC_test.t_rebar.size_list, 32)
# RC_test.c_rebar.size_counter = min_counter(RC_test.c_rebar.size_list, 12)
# RC_test.NA = 300

# RC_test.set_values_collapse(1000)
# print(RC_test.NA)
# print(RC_test.MRd)


#print(h_min_counter)
print("end")


# def find_min_rebars(RC, M):
#     Mp = 0.167*RC.concrete.fck*RC.b*RC.d**2
#     Kp = Mp/(RC.b*RC.concrete.fck*RC.d**2)
#     K = M/(RC.b*RC.concrete.fck*RC.d**2)
    
#     if M > Mp:
#         A_sp = ((K-Kp)*RC.b*RC.concrete.fck*RC.d**2)/(RC.c_rebar_stress*(RC.d-RC.dp))
#         A_s = (Kp*RC.b*RC.concrete.fck*RC.d**2)/(RC.t_rebar_stress*0.82*RC.d)
#     else:
#         A_sp = 0
#         A_s = M/(RC.t_rebar_stress*RC.z)
#     return [A_s, A_sp]

# #find minimum design of rebar that fits area
# def find_min_rebar_design(RC,rebar,area,min_spacing):
#     #returns output [no., d]
#     #min spacing is limited by vibrator size, can be arbitrary
#     if area == 0:
#         return [0,0]
    
#     initial_size_counter = rebar.size_counter
#     #initial_size = rebar.size
#     initial_no = rebar.no
    
#     #size_counter = initial_size_counter
#     no = initial_no
    
#     max_spacing = RC.max_bar_spacing
#     max_d = RC.max_bar_d
    
#     min_area = 0
#     min_pair = []
    
#     if rebar.area == 0:
#         i_start = initial_size_counter + 1
#     else:
#         i_start = initial_size_counter
    
#     for i in range(i_start, len(rebar.size_list)):
#         rebar.no = initial_no
#         while rebar.area < area:
#             rebar.no = rebar.no + 1
#         #assume that distance from side is same as distance from edge
#         spacing = (RC.b - 2*(RC.h-RC.d))/(no-1)
        
#         if rebar.usage == "compression":
#             #compression only consider min
#             if spacing >= min_spacing:
#                 if min_area == 0 or rebar.area < min_area:
#                     min_area = rebar.area
#                     min_pair = [rebar.no,rebar.d]
#                 elif rebar.area == min_area:
#                     #prioritise the one with more rebars for better bonding
#                     if min_pair[0] < rebar.no:
#                         min_pair = [rebar.no,rebar.d]
#         elif rebar.usage == "tension":
#             #tension consider min and max
#             if spacing >= min_spacing and spacing <= max_spacing and rebar.d <= max_d:
#                 if min_area == 0 or rebar.area < min_area:
#                     min_area = rebar.area
#                     min_pair = [rebar.no,rebar.d]
#                 elif rebar.area == min_area:
#                     #prioritise the one with more rebars for better bonding
#                     if min_pair[0] < rebar.no:
#                         min_pair = [rebar.no,rebar.d]
#     #reset the rebar. Required to do this because I did not think through the whole OOP thing :D
#     rebar.size_counter = initial_size_counter
#     rebar.no = initial_no
    
#     return min_pair
        
