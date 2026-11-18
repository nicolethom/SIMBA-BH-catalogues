"""
Combines bhallanalyse-2.py and bhar_hist.py.

For each snapshot:
  1. match the central SMBH (most massive BH) to each caesar galaxy
  2. find the bhALL event closest to the snapshot redshift for each SMBH
  3. average the accretion rates over 50 and 100 Myr, starting from that event
  4. write everything to one per-galaxy file

bhALL is never loaded whole. It is read in blocks of CHUNK_ROWS rows, and only
the rows for central SMBHs inside the redshift window are kept.
"""
import os
import time

import numpy as np
import h5py as h5
import readgadget as rg
from astropy import units as u
from scipy.integrate import quad
from scipy.optimize import brentq

start = time.time()

#-----------------settings-----------------------#

box = 'm100n1024'
sim = 'simbac'
snaps = np.arange(151,152,1)
outdir = '/idia/projects/simba/nthomas/simbac' #directory to save bh file

sims = {
    'simbac': dict(
        simdir='/idia/data/laduma/SIMBAC/%s/'%box,
        caesar='Groups/simba-c_%03i.hdf5',
        snapshot='snapshot_%03i.hdf5',
        z_group='simulation_attributes/parameters',
        z_key='Redshift',
    ),
    's50': dict(
        simdir='/idia/data/laduma/SIMBA/%s/s50/'%box,
        caesar='Groups/'+box+'_%03i.hdf5',
        snapshot='snap_'+box+'_%03i.hdf5', #not in either original script, check this
        z_group='simulation_attributes',
        z_key='redshift',
    ),
}

# rows of bhALL held in memory at once. Peak use per dataset is about
# CHUNK_ROWS x 8 bytes (x 3 for vectors), so lower this if RAM is still tight
CHUNK_ROWS = 2**22

z_window = 0.01 #redshift half-width around the snapshot
avg_times = {'50MYR':0.05,'100MYR':0.1} #averaging bins in Gyr

H0 = 68*(u.km*(u.s**-1)*(u.Mpc**-1))
th = (1./H0).to('Gyr').value

#-----------------functions-----------------------#

def t_L_int(z):
    den1 = (1+z)
    den2 = np.sqrt((0.3*(1+z)**3)+0.7)
    return 1./(den1*den2)

def lookback(z):
    return th*quad(t_L_int,0,z)[0]

def averaging_limit(z):
    # highest redshift the longest average can reach. z_window alone is too
    # short for this beyond z~0.3, where 0.01 in z is less than 100 Myr
    target = lookback(z+z_window)+max(avg_times.values())
    return brentq(lambda zz: lookback(zz)-target,z,1000.)

def lookback_table(z,z_far):
    # the range is narrow, so quad is run once on a fine grid and events are
    # interpolated from it rather than integrating for every event
    zgrid = np.linspace(z-2*z_window,z_far,2001)
    lbt = np.asarray([lookback(i) for i in zgrid])
    return zgrid,lbt

def runningmedian(x,y,xlolim,ylolim,dt):
    # mean of y in the first bin of width dt, starting at the smallest x
    xp = x[(x>=xlolim)&(y>=ylolim)]
    yp = y[(x>=xlolim)&(y>=ylolim)]
    bin_lo = np.min(x)
    bin_hi = bin_lo+dt
    ysub = yp[(xp>=bin_lo)&(xp<=bin_hi)]
    ymean = np.mean(ysub)
    return bin_lo+0.5*dt,ymean

def read_caesar(cfile,z_group,z_key):
    with h5.File(cfile,'r') as f:
        print('caesar file read')
        gals = dict(
            mbh=np.asarray(f['galaxy_data/dicts/masses.bh']),
            bhlist_start=np.asarray(f['galaxy_data/bhlist_start']),
            bhlist_end=np.asarray(f['galaxy_data/bhlist_end']),
            bhlist=np.asarray(f['galaxy_data/lists/bhlist']),
            z=np.round(float(f[z_group].attrs[z_key]),2),
        )
    return gals

def central_bh_ids(rgsnap,gals):
    # matches the black holes in each caesar galaxy to the snapshot and takes
    # the most massive as the central SMBH. A galaxy with no BHs gets ID -1
    bhs = rg.readsnap(rgsnap,'ParticleIDs','bndry',units=1,suppress=1)
    bhmass = rg.readsnap(rgsnap,'BH_Mass','bndry',suppress=1)*1e10/0.68
    print('snapshot data read')

    ngal = len(gals['bhlist_start'])
    bhID_per_gal = np.full(ngal,-1,dtype=np.int64)
    for gal in range(ngal):
        inds = gals['bhlist'][gals['bhlist_start'][gal]:gals['bhlist_end'][gal]]
        if len(inds)>0:
            bhID_per_gal[gal] = bhs[inds[np.argmax(bhmass[inds])]]
    return bhID_per_gal

def block_a_ranges(a_dset):
    # one pass over the scale factors, recording the range in each block, so
    # later steps can skip the blocks that cannot hold a snapshot's window
    block_starts = np.arange(0,a_dset.shape[0],CHUNK_ROWS)
    a_min = np.zeros(len(block_starts))
    a_max = np.zeros(len(block_starts))
    for i,lo in enumerate(block_starts):
        a = a_dset[lo:lo+CHUNK_ROWS]
        a_min[i],a_max[i] = a.min(),a.max()
    return block_starts,a_min,a_max

def collect_events(f,blocks,a_lim1,a_lim2,wanted):
    # returns the BH events in the window for the IDs in wanted (sorted),
    # grouped by BH. Row order is kept within each BH
    block_starts,a_min,a_max = blocks
    bh_ind = [np.empty(0,dtype=np.int64)]
    rows = [np.empty(0,dtype=np.int64)]
    a_ev = [np.empty(0)]
    mdot = [np.empty(0)]
    bondi = [np.empty(0)]

    if len(wanted)>0:
        for lo in block_starts[(a_max>a_lim1)&(a_min<a_lim2)]:
            block = slice(lo,lo+CHUNK_ROWS)
            a = f['a'][block]
            local = np.flatnonzero((a>a_lim1)&(a<a_lim2))
            if len(local)==0:
                continue

            ids = f['ID'][block][local].astype(np.int64)
            pos = np.minimum(np.searchsorted(wanted,ids),len(wanted)-1)
            match = wanted[pos]==ids
            local = local[match]
            if len(local)==0:
                continue

            bh_ind.append(pos[match])
            rows.append(local+lo)
            a_ev.append(a[local].astype(np.float64))
            mdot.append(f['Mdot'][block][local].astype(np.float64))
            bondi.append(f['Mdot_bondi'][block][local].astype(np.float64))

    bh_ind = np.concatenate(bh_ind)
    order = np.argsort(bh_ind,kind='stable')
    offsets = np.searchsorted(bh_ind[order],np.arange(len(wanted)+1))
    return (offsets,np.concatenate(rows)[order],np.concatenate(a_ev)[order],
            np.concatenate(mdot)[order],np.concatenate(bondi)[order])

def gather_rows(dset,rows):
    # reads the given sorted, unique rows from a dataset one block at a time,
    # as h5py is very slow at reading long lists of scattered rows
    out = np.zeros((len(rows),)+dset.shape[1:],dtype=dset.dtype)
    for block in np.unique(rows//CHUNK_ROWS):
        lo = block*CHUNK_ROWS
        first,last = np.searchsorted(rows,[lo,lo+CHUNK_ROWS])
        out[first:last] = dset[lo:lo+CHUNK_ROWS][rows[first:last]-lo]
    return out

def analyse_events(events,z,a_,a_near,z_far):
    # for each SMBH finds the event closest to the snapshot (from events with
    # a>a_near), then averages the accretion rates starting from that event.
    # An SMBH with nothing to average stays nan
    offsets,rows,a_ev,mdot,bondi = events
    nbh = len(offsets)-1
    closest_row = np.full(nbh,-1,dtype=np.int64)
    aves = {}
    for name in avg_times:
        aves['BHAR_'+name] = np.full(nbh,np.nan)
        aves['BHAR_BONDI_'+name] = np.full(nbh,np.nan)
    zgrid,lbt_grid = lookback_table(z,z_far)

    for bh in range(nbh):
        ev = slice(offsets[bh],offsets[bh+1])
        a_bh = a_ev[ev]
        near = np.flatnonzero(a_bh>a_near)
        if len(near)==0:
            continue

        closesta = near[np.argmin(np.abs(a_-a_bh[near]))]
        closest_row[bh] = rows[ev][closesta]

        zs = (1./a_bh)-1.
        redshift_late = min(zs[closesta],z)
        if not np.any((zs>=redshift_late)&(zs<redshift_late+z_window)):
            continue

        condi = zs>=redshift_late
        zs = zs[condi]
        BHARs = mdot[ev][condi]
        BHAR_B = bondi[ev][condi]

        lbt = np.interp(zs,zgrid,lbt_grid)
        for name,avg_time in avg_times.items():
            aves['BHAR_BONDI_'+name][bh] = runningmedian(lbt,BHAR_B,0,0,avg_time)[1]
            aves['BHAR_'+name][bh] = runningmedian(lbt,BHARs,0,0,avg_time)[1]

    return closest_row,aves

#-----------------main-----------------------#

cfg = sims[sim]
bhalldir = cfg['simdir']+'blackhole_details/bhALL.hdf5'
print('Directories read')

with h5.File(bhalldir,'r',rdcc_nbytes=0) as f:
    print('bhall file read')
    f_keys = list(f.keys())
    #['BH_Mass', 'ID', 'Jgas', 'Jstar', 'Mass', 'Mass_AlphaDisk', 'Mdot', 'Mdot_alphadisk', 'Mdot_bondi', 'Mgas', 'MgasBulge', 'Mhot', 'Mstar', 'MstarBulge', 'R0', 'Sfr', 'a', 'dt', 'e', 'p', 'rho', 'v'])
    blocks = block_a_ranges(f['a'])

    for snap in snaps:
        print(snap)
        cfile = cfg['simdir']+cfg['caesar']%(snap) #caesar file
        rgsnap = cfg['simdir']+cfg['snapshot']%(snap) #snapshot

        #---------------galaxies and their SMBHs-------------#

        gals = read_caesar(cfile,cfg['z_group'],cfg['z_key'])
        z = gals['z']
        mbh = gals['mbh']
        ngal = len(mbh)
        print('redshift=%1.2f'%(z))
        print('matching BHs IDs for %g caesar galaxies'%(ngal))
        bhID_per_gal = central_bh_ids(rgsnap,gals)

        #---------------BH events near the snapshot-------------#

        # select the black hole events within a redshift bin, extended to
        # higher redshift far enough to cover the longest average
        print('finding BH events')
        a_ = 1./(1.+z)
        a_near = 1./(1.+z+z_window)
        a_lim2 = 1./(1.+z-z_window)
        z_far = averaging_limit(z)
        a_lim1 = 1./(1.+z_far)

        has_bh = bhID_per_gal>=0
        wanted = np.unique(bhID_per_gal[has_bh])
        events = collect_events(f,blocks,a_lim1,a_lim2,wanted)
        print('%i events for %i SMBHs'%(len(events[1]),len(wanted)))

        print('identifying closest BH event and averaging BHAR')
        closest_row,aves = analyse_events(events,z,a_,a_near,z_far)
        del events

        # each galaxy points to its SMBH's place in wanted
        gal_ind = np.flatnonzero(has_bh)
        bh_of_gal = np.searchsorted(wanted,bhID_per_gal[gal_ind])

        missing = closest_row[bh_of_gal]<0
        if missing.any():
            print('WARNING: %i SMBHs have no events in the redshift window'%(missing.sum()))
        found_gal = gal_ind[~missing]
        found_rows = closest_row[bh_of_gal[~missing]]
        uniq_rows,row_of_gal = np.unique(found_rows,return_inverse=True)

        #---------------write the analysis-------------#

        # If a galaxy hosts no black hole, it takes black ID nan. Its property
        # values are nan, as are those of an SMBH with no events in the window
        outfile = os.path.join(outdir,'bhall_details_%03i.hdf5'%(snap))
        with h5.File(outfile,'w') as g:
            for key in f_keys:
                col = gather_rows(f[key],uniq_rows)
                data = np.full((ngal,)+f[key].shape[1:],np.nan)
                data[found_gal] = col[row_of_gal]
                if key=='ID':
                    data[gal_ind] = bhID_per_gal[gal_ind]
                g.create_dataset(key,data=data)

            accreting = mbh[gal_ind]>0.
            for name,bh_ave in aves.items():
                data = np.full(ngal,np.nan)
                data[gal_ind[accreting]] = bh_ave[bh_of_gal[accreting]]
                g.create_dataset(name,data=data)

        print('snap %03i complete'%(snap))

end = time.time()
end = (end-start)/60.

print('Complete in %f minutes'%end)
