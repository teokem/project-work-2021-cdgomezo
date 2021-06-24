from numpy import *

def modelCO2(initCond, fluxes):

    # Constants
    M = 5.137e21     # mass of the atmosphere [g]
    Ma = 28.964      # molar mass of air [g/mol]
    MC = 12          # molar mass of carbon [g/mol]
    r = MC*(M/2/Ma)  # conversion from molar fraction to gC
    dt = 1.
    tau = 0.5

    # Create a list containing just c0
    conc = zeros((len(fluxes['time'])+1, 2))
    conc[0,:] = array((initCond['co2']['nh'], initCond['co2']['sh']))

    for i in range(len(fluxes['time'])):
        
        conc[i+1, 0] = conc[i, 0] + (tau * (-conc[i, 0] + conc[i, 1]) + 
                                     (fluxes['co2']['nh']['bio'][i] + fluxes['co2']['nh']['oce'][i] + fluxes['co2']['nh']['ff'][i]) / r) * dt
        conc[i+1, 1] = conc[i, 1] + (tau * (conc[i, 0] - conc[i, 1]) + 
                                     (fluxes['co2']['sh']['bio'][i] + fluxes['co2']['sh']['oce'][i] + fluxes['co2']['sh']['ff'][i]) / r) * dt
    
    obs = {
        'units' : initCond['units'],
        'co2': {
            'nh': {
                'time': fluxes['time'],
                'value': conc[1:,0].reshape(-1,1)
            },
            'sh': {
                'time': fluxes['time'],
                'value': conc[1:,1].reshape(-1,1)
            }
        }
    }

    return obs