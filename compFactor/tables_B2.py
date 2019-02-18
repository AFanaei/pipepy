from numpy import array , zeros
# These are Characterization parameters for Calculation of state compression factor/Composition analysis
# Identification number for component list as a Dic type
IdNum = {0:'Methane' , 1:'Nitrogen' , 2:'Carbon dioxide' , 3:'Ethane' , 4:'Propane' , 5:'Water' , 6:'Hydrogen sulfide' , 7:'Hydrogen',
         8:'Carbon monoxide' , 9:'Oxygen' , 10:'iso_Butane', 11:'n-Butane', 12:'iso_Pentane', 13:'n_Pentane' , 14:'n_Hexane', 15:'n_Heptane',
         16:'n_Octane' , 17:'n_Nonane' , 18:'n_Decane', 19:'Helium', 20:'Argon' }
# Molar mass [kg.kmol^-1]
M = array([16.0430, 28.0135, 44.0100, 30.0700, 44.0970, 18.0153, 34.0820, 2.0159, 28.0100, 31.9988, 58.1230, 58.1230, 72.1500, 72.1500,
     86.1770, 100.2040, 114.2310, 128.2580, 142.2850, 4.0026, 39.9480])
# Energy parameter [K]
E = array([151.318300, 99.737780, 241960600, 244.166700, 298.118300, 514.015600, 296.355000, 26.957940, 105.534800, 122.766700, 324.068900,
     337.638900, 365.599900, 370.628300, 402.636293, 427.722630, 450.325022, 470.840891, 489.558373, 2.610111, 119.629900 ])
# Size parameter [(m^3/kmol)^(1/3)]
K = array([.4619255, .4479153, .4557489, .5279209, .5837490, .3825868, .4618263, .3514916, .4533894, .4186954, .6406937, .6341423, .6738577,
     .6798307, .7175118, .7525189, .7849550, .8152731, .8437826, .3589888, .4216551 ])
# Orientation parameter
G = array([0.0, .027815, .189065, .079300, .141239, .332500, .088500, .034369, .038953, .021000, .256692, .281835, .332267, .366911, .289731,
     .337542, .383381, .427354, .469659, 0.0, 0.0 ])
# Quadru pole parameter
Q = zeros(21) ; Q[2] = .690000 ; Q[5] = 1.067750 ; Q[6] = .633276  
# High_temp.parameter
F = zeros(21) ; F[7] = 1.0 
# Dipole parameter
S = zeros(21) ; S[5] = 1.582200 ; S[6] = .390000 
# Association parameter
W = zeros(21) ; W[5] = 1.0