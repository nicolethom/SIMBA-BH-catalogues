import numpy as np
import h5py as h5
import readgadget as rg
import os
import time
start = time.time()
box = 'm100n1024'
snaps=np.arange(127,128,1)
snapdir = '/idia/data/laduma/SIMBA/%s/s50/'%box #snapshot directory
bhalldir = snapdir+'blackhole_details/bhALL.hdf5' #bhall file
#source directories
for snap in snaps:
    cfile = snapdir+'Groups/%s_%03i.hdf5'%(box,snap) #caesar file
    rgsnap = snapdir+'snap_%s_%03i.hdf5'%(box,snap) #snapshot

    #directory to save bh file
    outdir = '/users/nthomas/data/Groups/'
    #-----------------caesar info-----------------------#

    f = h5.File(cfile,'r')

    ms = np.asarray(f['galaxy_data/dicts/masses.stellar'])
    mbh = np.asarray(f['galaxy_data/dicts/masses.bh'])
    mbhdot = np.asarray(f['galaxy_data/bhmdot'])
    bhlist_start = np.asarray(f['galaxy_data/bhlist_start'])
    bhlist_end = np.asarray(f['galaxy_data/bhlist_end'])
    bhlist = np.asarray(f['galaxy_data/lists/bhlist'])
    z = np.asarray(h5.AttributeManager(f['simulation_attributes/parameters'])['Redshift'])
    z = np.round(z,2)
    print('redshift=%1.2f'%(z))

    # this bit matches all the black hole indices to each galaxy in caesar
    # ie selects black holes from the bhlist 
    print('matching BHs IDs for %g caesar galaxies'%(len(ms)))
    bhInd_per_gal = []

    for bh in range(len(bhlist_start)):
        bhs = bhlist[bhlist_start[bh]:bhlist_end[bh]]
        bhInd_per_gal.append(bhs)

    bhInd_per_gal = np.asarray(bhInd_per_gal)

    f.close()

    #---------------snapshot inf0----------------#

    bhs = rg.readsnap(rgsnap,'ParticleIDs','bndry',units=1,suppress=1)
    bhmass = rg.readsnap(rgsnap,'BH_Mass','bndry',suppress=1)*1e10/0.68

    # mapping central black holes to caesar galaxies
    # this bit matches the black holes previously identified by index
    # to the black holes in the simba snapshot and allocates the 
    # appropriate black hole IDs to their host galaxy in caesar
    print('matching BHs from snapshot to caesar galaxies')
    bhID_per_gal=[]
    bhmasses_per_gal=[]
    for bh in range(len(bhInd_per_gal)):
        bhID_per_gal.append(bhs[bhInd_per_gal[bh]])
        bhmasses_per_gal.append(bhmass[bhInd_per_gal[bh]])

    bhmasses_per_gal = np.asarray(bhmasses_per_gal)
    bhID_per_gal = np.asarray(bhID_per_gal)

    print('identifying central SMBH')
    # it then identifies the largest black hole of the associated black holes
    # as the central supermassive black hole
    # bhID_per_gal takes value the ID of the most massive BH
    for bhmass in range(len(bhmasses_per_gal)):
        if len(bhmasses_per_gal[bhmass])>0:
            maxbhmass = np.max(bhmasses_per_gal[bhmass])
            maxbhpos = np.where(bhmasses_per_gal[bhmass]==maxbhmass)[0][0]
            bhID_per_gal[bhmass]=int(bhID_per_gal[bhmass][maxbhpos])
            bhmasses_per_gal[bhmass]=bhmasses_per_gal[bhmass][maxbhpos]
        else:
            bhmasses_per_gal[bhmass]=0.
            bhID_per_gal[bhmass]=np.nan
    # if a galaxy has no black holes
    # then of course it has no black hole mass nor ID!

    #---------------bhALL info-------------#

    # We want to take all the info from the black hole details file
    # and match it so the respective SMBHs
    f = h5.File(bhalldir)
    f_keys = np.asarray(['BH_Mass', 'ID', 'Jgas', 'Jstar', 'Mass', 'Mass_AlphaDisk', 'Mdot', 'Mdot_alphadisk', 'Mdot_bondi', 'Mgas', 'MgasBulge', 'Mhot', 'Mstar', 'MstarBulge', 'R0', 'Sfr', 'a', 'dt', 'e', 'p', 'rho', 'v'])

    for param in f_keys:
        var1 = param
        exec(var1+"=np.asarray(f[param])")

    f.close()

    # To do so we need to find the closest black hole event to the snapshot redshift
    print('finding BH events')
    z_lim1,z_lim2 = z+0.01,z-0.01
    a_ = (1./(1.+z))
    a_lim1 = (1./(1.+z_lim1))
    a_lim2 = (1./(1.+z_lim2))

    local = (a>a_lim1)&(a<a_lim2)
    print(local)
    # First we select all the black hole events within a redshift bin

    local_BH_Mass = BH_Mass[local]
    local_ID = ID[local]
    local_Jgas = Jgas[local]
    local_Jstar = Jstar[local]
    local_Mass = Mass[local]
    local_Mass_AlphaDisk = Mass_AlphaDisk[local]
    local_Mdot = Mdot[local]
    local_Mdot_alphadisk = Mdot_alphadisk[local]
    local_Mdot_bondi = Mdot_bondi[local]
    local_Mgas = Mgas[local]
    local_MgasBulge = MgasBulge[local]
    local_Mhot = Mhot[local]
    local_Mstar = Mstar[local]
    local_MstarBulge = MstarBulge[local]
    local_R0 = R0[local]
    local_Sfr = Sfr[local]
    local_a = a[local]
    local_dt = dt[local]
    local_e = e[local]
    local_p = p[local]
    local_rho = rho[local]
    local_v = v[local]

    #GET READY TO FIND THE SMBH!

    BH_Mass = np.zeros(len(ms))
    ID= np.zeros(len(ms))
    Jgas= np.zeros((len(ms),3))
    Jstar= np.zeros((len(ms),3))
    Mass= np.zeros(len(ms))
    Mass_AlphaDisk= np.zeros(len(ms))
    Mdot= np.zeros(len(ms))
    Mdot_alphadisk= np.zeros(len(ms))
    Mdot_bondi= np.zeros(len(ms))
    Mgas= np.zeros(len(ms))
    MgasBulge= np.zeros(len(ms))
    Mhot= np.zeros(len(ms))
    Mstar= np.zeros(len(ms))
    MstarBulge= np.zeros(len(ms))
    R0= np.zeros(len(ms))
    Sfr= np.zeros(len(ms))
    a= np.zeros(len(ms))
    dt= np.zeros(len(ms))
    e= np.zeros(len(ms))
    p= np.zeros((len(ms),3))
    rho= np.zeros(len(ms))
    v= np.zeros((len(ms),3))
    count=0

    # for each galaxy given it hosts a SMBH, we locate the SMBH in the redshift slice
    # if there more than one entries, the properties of the SMBH take value of the
    # entry closest to the snapshot redshift
    print('identifying closest BH event')
    for bh in range(len(bhID_per_gal)):
        if np.isfinite(bhID_per_gal[bh])==True:
            pos = np.where(local_ID==bhID_per_gal[bh])[0]
            print(bhID_per_gal[bh])
            print(pos)
            closesta = pos[0]
            if len(pos)>1:
                for sf in range(len(pos)):
                    if abs(a_-local_a[pos[sf]])<abs(a_-local_a[closesta]):
                        closesta = pos[sf]
                    else:
                        pass
            else:
                closesta = pos[0] 
            pos = np.copy(closesta)

            BH_Mass[bh] = local_BH_Mass[pos]
            ID[bh] = int(local_ID[pos])
            Jgas[bh] = local_Jgas[pos]
            Jstar[bh] = local_Jstar[pos]
            Mass[bh] = local_Mass[pos]
            Mass_AlphaDisk[bh] = local_Mass_AlphaDisk[pos]
            Mdot[bh] = local_Mdot[pos]
            Mdot_alphadisk[bh] = local_Mdot_alphadisk[pos]
            Mdot_bondi[bh] = local_Mdot_bondi[pos]
            Mgas[bh] = local_Mgas[pos]
            MgasBulge[bh] = local_MgasBulge[pos]
            Mhot[bh] = local_Mhot[pos]
            Mstar[bh] = local_Mstar[pos]
            MstarBulge[bh]= local_MstarBulge[pos]
            R0[bh]= local_R0[pos]
            Sfr[bh]= local_Sfr[pos]
            a[bh]= local_a[pos]
            dt[bh]= local_dt[pos]
            e[bh]= local_e[pos]
            p[bh]= local_p[pos]
            rho[bh]= local_rho[pos]
            v[bh]= local_v[pos]

        else:
            ID[bh] = np.nan
    # If a galaxy hosts no black hole, it takes black ID nan
    # and property values = 0.

    outfile = outdir+'BH_Bondi_%s_%03i.hdf5'%(box,snap)
    try:
        os.system('rm %s'%outfile)
        print('file already exists, rewriting bh analysis')
    except:
        print('writing analysis')

    g = h5.File(outfile,'a')
    names = np.copy(f_keys)
    datas = [BH_Mass,ID,Jgas,Jstar,Mass,Mass_AlphaDisk,Mdot,Mdot_alphadisk,Mdot_bondi,Mgas,MgasBulge,Mhot,Mstar,MstarBulge,R0,Sfr,a,dt,e,p,rho,v]

    for i in range(len(datas)):
        #print(type(datas[i]))
        g.create_dataset(names[i],data=datas[i])
    g.close()
    print('Crossmatching SMBH event details to caesar galaxies snapshot %i complete!'%(snap))
end = time.time()
end = (end-start)/60.

print('Complete in %f minutes'%end)







