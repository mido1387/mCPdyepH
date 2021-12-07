"""
Spectrophotometer meta-Cresol Purple dye pH calculation
Calculations from Liu et al, Environmental Science & Technology 2011 45 (11), 4862-4868
Works for Salinity: 20 < S < 40 psu
Works for Temperature: 278 < T < 308  K
Author: Michael R Dooley
"""

def dye_pH(A578, A434, A730, Vdye, T, S, dpa, dpb):

    """
    pH of sample
    % function es = dye_pH(A578, A434, A730, Vdye, T, S),
    %- INPUT:
    %-  1) Absorbance at 578 nm
    %-  2) Absorbance at 434 nm
    %-  3) Absorbance at 730 nm
    %-  4) Volume of dye added in micro liters
    %-  5) Temperature in Kelvin
    %-  6) Salinity in psu
    %-  7) Dye perturbation line constant a
    %-  8) Dye perturbation line constant b
    %- OUTPUT:
    %-  1) pH of sample.  
    %- 
    """
    
    import numpy as np
    
    A578corr=A578-A730
    A434corr=A434-A730
    
    if dpa != 0:
        Absratio=(A578corr/A434corr) - Vdye * (dpa+dpb*(A578corr/A434corr))
        e1=-0.007762+4.5174*(10**-5)*T
        e2e3=-0.020813+2.60262*(10**-4)*T+1.0436*(10**-4)*(S-35)
        a=-246.64209+0.315971*S+2.8855*(10**-4)*S**2
        b=7229.23864-7.098137*S-0.057034*S**2
        c=44.493382-0.052711*S
        d=0.0781344
        pH=a+(b/T)+c*np.log(T)-d*T+np.log10((Absratio-e1)/(1-Absratio*e2e3))
    else:
        Absratio=(A578corr/A434corr)
        e1=-0.007762+4.5174*(10**-5)*T
        e2e3=-0.020813+2.60262*(10**-4)*T+1.0436*(10**-4)*(S-35)
        a=-246.64209+0.315971*S+2.8855*(10**-4)*S**2
        b=7229.23864-7.098137*S-0.057034*S**2
        c=44.493382-0.052711*S
        d=0.0781344
        pH=a+(b/T)+c*np.log(T)-d*T+np.log10((Absratio-e1)/(1-Absratio*e2e3))
    return pH