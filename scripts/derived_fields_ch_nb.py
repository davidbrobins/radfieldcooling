# Module to write the functions necessary to create n_b, Gamma, Lambda derived fields
# Import as 'from derived_fields_ch_nb import *
# Still need to actually create the derived fields with the commented out lines below
import yt
import numpy as np
from scipy.io import FortranFile
from yt.units import gram, second, erg, centimeter
import CF
table=FortranFile('cf_table.I2.dat', 'r', header_dtype=np.int32) #open cf_table.I2.dat and note 4-byte header and footer
constants = table.read_record([('lt','i4'),('ld','i4'),('np','i4',3),('lp4','i4'), ('qmin','f4',3), ('q1', 'f4'), ('qmax','f4',3),('q2', 'f4'),('lx', 'i4'),('xmin', 'f4'),('xmax', 'f4')]) #read in the values of several constants (int/float noted from frt_cf3.F)
altval=table.read_record(np.float32) #Read in array of real temperature values
indx=table.read_record(np.int32).reshape((24,21,16), order='F') #Read in and reshape 3-D array of indices
data=np.zeros((6,81,13,3789)) #Create array to hold data values
for i in range(3789):
    data[:,:,:,i]=table.read_record(np.float32).reshape((6,81,13), order='F') #Read in 3-D array of data values at given 4th index
table.close()
CF.frtinitcf(0,constants['lt'][0], constants['ld'][0],constants['np'][0], constants['lp4'][0],constants['qmin'][0],constants['q1'][0],constants['qmax'][0],constants['q2'][0],constants['lx'][0],constants['xmin'][0],constants['xmax'][0], altval, indx, data) #initialize CF with values from cf_table.I2.dat

proton_mass = 1.67262192369e-24*gram #proton mass in grams
def rho_to_n_b(field, data):
    return data["density"]/proton_mass #Convert from mass density to baryon number density by dividing by proton mass (approximate in general, exact for all H
yt.add_field(('gas', 'baryon_number_density'), function=rho_to_n_b, units='1/cm**3')
def cooling_func_from_array(input_array): #get cooling function from array of input values (for use in derived field calculation)
    T=input_array[0]
    n_b=input_array[1]
    Z=input_array[2]
    P_LW=input_array[3]
    P_HI=input_array[4]
    P_HeI=input_array[5]
    P_CVI=input_array[6] #extract the 7 input floats from the array
    (cfun,hfun,ierr)=CF.frtgetcf(T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI) #call f2py function to get cooling and heating functions for these inputs
    return cfun #output cooling function
def cooling_rate(field, data): #derived field function to get cooling function
    data_array = [] #create array to house the needed data
    T=data['gas', 'temperature']
    data_array.append(T.to_value()) #append input fields as numpy arrays
    n_b=data['gas',"baryon_number_density"]
    data_array.append(n_b.to_value())
    Z=data['gas',"metallicity"]
    data_array.append(Z.to_value())
    P_LW=data['artio', 'RT_DISK_VAR_0']
    data_array.append(P_LW.to_value()) #The 4 photoionization rates should be in 1/s, but yt outputs them as unitless
    P_HI=data['artio', 'RT_DISK_VAR_1']
    data_array.append(P_HI.to_value())
    P_HeI=data['artio', 'RT_DISK_VAR_2']
    data_array.append(P_HeI.to_value())
    P_CVI=data['artio', 'RT_DISK_VAR_3']
    data_array.append(P_CVI.to_value())
    return np.apply_along_axis(cooling_func_from_array, 0, data_array)*erg*centimeter**3/second #Calculate cooling_func_from_array using ith element from each input array
yt.add_field(('gas', 'cooling_rate'), function=cooling_rate, units='erg*cm**3/s')
def heating_func_from_array(input_array): #get heating function from array of input values (for use in derived field calculation)                                     
    T=input_array[0]
    n_b=input_array[1]
    Z=input_array[2]
    P_LW=input_array[3]
    P_HI=input_array[4]
    P_HeI=input_array[5]
    P_CVI=input_array[6] #extract the 7 input floats from the array                                                                                                   
    (cfun,hfun,ierr)=CF.frtgetcf(T, n_b, Z, P_LW, P_HI, P_HeI, P_CVI) #call f2py function to get cooling and heating functions for these inputs                       
    return hfun #output heating function                                                                                                      
def heating_rate(field, data): #derived field function to get heating function                                                                                   
    data_array = [] #create array to house the needed data                                                                                                            
    T=data['gas', 'temperature']
    data_array.append(T.to_value()) #append input fields as numpy arrays                                                                                              
    n_b=data['gas',"baryon_number_density"]
    data_array.append(n_b.to_value())
    Z=data['gas',"metallicity"]
    data_array.append(Z.to_value())
    P_LW=data['artio', 'RT_DISK_VAR_0']
    data_array.append(P_LW.to_value()) #The 4 photoionization rates should be in 1/s, but yt outputs them as unitless                                          
    P_HI=data['artio', 'RT_DISK_VAR_1']
    data_array.append(P_HI.to_value())
    P_HeI=data['artio', 'RT_DISK_VAR_2']
    data_array.append(P_HeI.to_value())
    P_CVI=data['artio', 'RT_DISK_VAR_3']
    data_array.append(P_CVI.to_value())
    return np.apply_along_axis(heating_func_from_array, 0, data_array)*erg*centimeter**3/second #Calculate heating_func_from_array using ith element from each input array  
yt.add_field(('gas', 'heating_rate'), function=heating_rate, units='erg*cm**3/s')
