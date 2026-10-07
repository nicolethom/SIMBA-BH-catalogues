This repository includes python scripts to create additional catalogues for analysis of black hole demographics and histories in the Simba Simulations.
This includes the following:
1) bhallanalyse.py : crossmatches the bhALL.hdf5 file to any or all caesar galaxy catalogues.
                     Black holes are identified visa the ID of a galaxy's central SMBH and the event with scale factor closest in event to the snapshots redshift is identified as the selected black hole's properties
                     Output is a/a set of catalogues with black hole details as per bhALL and corresponding to the caesar snapshot number and in order of the caesar galaxy catalogue (i.e listed by group ID)
2) bhar_hist.py: Computes a time averaged black hole accretion rates and histories (using high time resolution bhALL file) for all central SMBHs in the caesar catalogues.
                 Time averaged quantities are written to the corresponding catalogues created by bhallanalyse.
3) radio_luminosity.py: Computes the 1.4GHz radio luminosities for all galaxies in caesar for star formation and AGN (using the KJF08 relation).
                        P1.4 is then written to the above-mentioned black hole catalogues.                       
