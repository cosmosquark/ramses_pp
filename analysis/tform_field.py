import yt
import numpy as np
from ramses_pp.analysis import tform


def get_friedman(ds):
    alpha=1e-6
    axpmin=1e-3
    ntable=1000

    if not "friedman" in ds.parameters:
        axp_out, hexp_out, tau_out, t_out, age_tot = tform.friedman(ds.omega_matter,
                                                                       ds.omega_lambda,
                                                                       ds.omega_curvature,
                                                                       alpha,
                                                                       axpmin,
                                                                       ntable)
        hubble_parameter = ds.arr(ds.hubble_constant, "100*km/s/Mpc")
        aexp = ds.parameters["aexp"]
        #Find neighbouring expansion factor
        i = 1
        while ( (axp_out[i] > aexp ) and (i < ntable) ):
            i+=1

        #Interpolate time
        time_simu = t_out[i] * (aexp - axp_out[i-1])/(axp_out[i]-axp_out[i-1]) + \
            t_out[i-1]*(aexp - axp_out[i])/(axp_out[i-1] - axp_out[i])

        time_simu = ds.arr(time_simu,"dimensionless")
        age_tot = ds.arr(age_tot,"dimensionless")
        age_simu = (time_simu+age_tot)/ hubble_parameter.in_units("m /(m*Gyr)") / ds.arr(1,"dimensionless")

        friedman = {
            'axp_out':ds.arr(axp_out,"dimensionless"),
            'hexp_out':ds.arr(hexp_out,"1/Gyr"),
            'tau_out':ds.arr(tau_out,"Gyr"),
            't_out':ds.arr(t_out,"Gyr"),
            'age_tot':age_tot,
            'age_simu':age_simu,
            'time_simu':time_simu,
            }        
        
        ds.parameters["friedman"] = friedman
        return

    else: # we already got the friedmann stuff
        return


def setup_tform_field(ad):
    def _tform(field, data):
        hubble_parameter = data.ds.arr(data.ds.hubble_constant, "100*km/s/Mpc")
        tform_val = tform.tform(data["particle_age"].in_units("Gyr").value, 
                              data.ds.parameters["friedman"]["tau_out"].in_units("Gyr").value, 
                              data.ds.parameters["friedman"]["t_out"].in_units("Gyr").value, 
                              data.ds.parameters["friedman"]["time_simu"].in_units("dimensionless").value, 
                              hubble_parameter.in_units("m /(m*Gyr)").value)
        tform_yt = data.ds.arr(tform_val,"Gyr")
        return tform_yt

    ad.ds.add_field(name = ("all", "particle_tform"),function=_tform, particle_type=True, units="Gyr")



    
