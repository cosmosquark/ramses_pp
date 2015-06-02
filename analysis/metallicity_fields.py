import yt
from yt import add_field
import numpy as np

metal_ref = "asplund"

def _compute_particle_met(field, data):
	"""
	Computes the metallicities from either a YT dataset
	Or a dictionary of YT Arrays
	"""

	## Anders & Grevesse 1989 as default values

	if metal_ref == None or metal_ref == "anders":
		sol_abund = {"particle_H":0.706, "particle_He":0.275, "particle_C":3.03e-3, "particle_N":1.11e-3, "particle_O":9.59e-3, "particle_Ne":0.00000001, "particle_Mg":5.15e-4, "paticle_Si":6.53e-4, "particle_Fe":1.17e-3,"particle_metallicity":0.019}

	else:
		# asplund data
		sol_abund = {"particle_H":0.715, "particle_He":0.270, "particle_C": 0.0024, "particle_N":0.000728 , "particle_O":0.006, "particle_Ne":0.0013, "particle_Mg":0.000742, "particle_Si":0.0007, "particle_Fe":0.00135,  "particle_metallicity":0.0143}

	print field.name[1], "test"
	if field.name[1] == "particle_FeH":
		metal = np.log10(data['all',"particle_Fe_fraction"].value / data['all',"particle_H_fraction"].value) - np.log10(sol_abund["particle_Fe"] / sol_abund["particle_H_fraction"])
	else:
		string = str(field.name[1])[:-2] #  strips the Fe from the end of the field
		metal = np.log10(data['all',string].value / data['all','particle_Fe_fraction'].value) - np.log10(sol_abund[string] / sol_abund["particle_Fe_fraction"])
	
	print metal
	return metal


def add_particle_fields(ds,He=False):

	ds.add_field(name = ('all','particle_FeH'),function=_compute_particle_met, units="dimensionless", particle_type=True)
	ds.add_field(name = ('all','particle_MgFe'),function=_compute_particle_met, units="dimensionless", particle_type=True)
	ds.add_field(name = ('all','particle_OFe'),function=_compute_particle_met, units="dimensionless", particle_type=True)
	ds.add_field(name = ('all','particle_CFe'),function=_compute_particle_met, units="dimensionless", particle_type=True)
	ds.add_field(name = ('all','particle_NFe'),function=_compute_particle_met, units="dimensionless", particle_type=True)
	ds.add_field(name = ('all','particle_SiFe'),function=_compute_particle_met, units="dimensionless", particle_type=True)

	return ds
