#!/bin/env python3

def calc_rg(trj):
                    mass = []
                    for at in trj.topology.atoms:
                        mass.append(at.element.mass)
                    mass_CA = len(mass)*[0.0]
                    # put the CA entries equal to 1.0
                    for i in trj.topology.select("name CA"):
                        mass_CA[i] = 1.0
                    # calculate CA radius of gyration
                    rg_CA = md.compute_rg(trj, masses=np.array(mass_CA))
                    return rg_CA


