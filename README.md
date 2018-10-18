# openmm_metadynamics
This is a toy metadynamics example using OpenMM. It takes advantage of the `CustomCentroidBondForce` which can use a tabulated function as input. This means the position of multiple molecules can add Hills to the bias at the same time. The power of this is that you can have multiple ligands with degenerate CVs (i.e. 10 ligands, all biasing x y z coordinates) and get extra sampling with a single system. Conceptually this is a blend of metadynamics and 'flooding' 

This would be an openmm version of this, which is only in plumed and not in an official release: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.8b00448
