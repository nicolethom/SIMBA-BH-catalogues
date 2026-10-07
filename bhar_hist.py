import numpy as np
import h5py as h5
import matplotlib
import pylab as plt
from astropy import units as u
from astropy import constants as cons
from scipy.integrate import quad
from scipy.misc import derivative
import time

start = time.time()

box = 'm100n1024'
suite = 's50'
snapnums = np.arange(140,150,1)
#snapnums=[151]
simdir = '/idia/data/laduma/SIMBA/%s/%s/'%(box,suite)
datadir = '/users/nthomas/data/Groups/'
imgdir = '/users/nthomas/images/Groups/%s/%s/'%(box,suite)
bhalldir = simdir+'blackhole_details/bhALL.hdf5'

H0 = 68*(u.km*(u.s**-1)*(u.Mpc**-1))
th = 1./H0

def t_L_int(z):
    den1 = (1+z)
    den2 = np.sqrt((0.3*(1+z)**3)+0.7)
    return 1./(den1*den2) 

def t_L(x):
    tL = th*x
    tL = tL.to('Gyr')
    return tL
'''
def top_ticks(zticks):
    LBtime_ticks = []
    for tick in range(len(zticks)):
        integral = quad(t_L_int,0,tick)
        lbt = t_L(integral)
        LBtime_ticks.append(lbt)
    LBtime_ticks = np.asarray(LBtime_ticks)
    return LBtime_ticks
'''  
dt = 0.05
def runningmedian(x,y,xlolim,ylolim,cs,dt):
        xp = x[(x>=xlolim)&(y>=ylolim)]
        yp = y[(x>=xlolim)&(y>=ylolim)]
        #nbins = np.ceil(max(x)/0.1)
        
        hist,bin_edges=np.histogram(xp,bins=np.arange(min(x),min(x)+2*dt,dt))
        #print(bin_edges)
        bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
        xsub = xp[xp>=bin_edges[0]]
        ysub = yp[xp>=bin_edges[0]]
        ysub = ysub[xsub<=bin_edges[1]]
        ymean = np.log10(np.mean(10**ysub))
        return bin_cent[0],ymean
        '''
        ymed = []
        ymean = []
        ysigma = []
        for i in range(0,len(bin_edges[:-1])):
                xsub = xp[xp>bin_edges[i]]
                ysub = yp[xp>bin_edges[i]]
                ysub = ysub[xsub<bin_edges[i+1]]
                ymed.append(np.median(10**ysub))
                ymean.append(np.mean(10**ysub))
                ysigma.append(np.std(10**ysub))    
        bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ymean = np.asarray(ymean)
        ymed = np.asarray(ymed)
        ysiglo = np.maximum(ymed-ysigma,ymed*0.1)
        ysiglo = np.log10(ymed)-np.log10(ysiglo)
        ysighi = np.log10(ymed+ysigma)-np.log10(ymed)
        #ymean = np.log10(ymean)
        '''    
        
for snapnum in snapnums:     
    print(snapnum)
    f = h5.File(simdir+'Groups/%s_%03i.hdf5'%(box,snapnum))
    ms = np.asarray(f['galaxy_data/dicts/masses.stellar'])
    mbh = np.asarray(f['galaxy_data/dicts/masses.bh'])
    sfr = np.asarray(f['galaxy_data/sfr'])
    galID = np.asarray(f['galaxy_data/GroupID'])
    redshift = np.asarray(h5.AttributeManager(f['simulation_attributes'])['redshift'])
    print('z=%2.2f'%redshift)
    redshift = np.round(redshift,2)
    f.close()
    ssfr = 1e9*(sfr/ms)

    bfile = h5.File(simdir+'Groups/BH_Bondi_%s_%03i.hdf5'%(box,snapnum),'r')
    
    bhID = np.asarray(bfile['ID'])
    bhdot = np.asarray(bfile['Mdot'])
    bh_bondi = np.asarray(bfile['Mdot_bondi'])
    scale_fac = np.asarray(bfile['a'])
    bh_gt = bhdot-bh_bondi
    bfile.close()

    z_in_file = (1./scale_fac)-1.
    
    fbondi = bh_bondi/bhdot

    G = cons.G
    mp = cons.m_p
    sigT = cons.sigma_T
    c = cons.c
    eta = 0.1

    factor = 4.*np.pi*G*mp/(eta*sigT*c)
    factor = factor.to('yr^-1')
    dotedd = factor*mbh*u.Msun
    fedd = bhdot/dotedd
    fedd = fedd.value


    BHID_rg = np.copy(bhID)#[rgs]

    mstar = ms#[rgs]
    mbhs = mbh#[rgs]
    galIDs = galID#[rgs] 

    bhall = h5.File(bhalldir,'r')

    BHIDs = np.asarray(bhall['ID'])
    BHAR = np.asarray(bhall['Mdot'])
    BHAR_BONDI = np.asarray(bhall['Mdot_bondi'])
    BHAR_GT = BHAR-BHAR_BONDI

    a = np.asarray(bhall['a'])
    z = (1./a)-1.
    

    bhall.close()

    IDs = np.copy(BHID_rg)
    lIDs = len(IDs)
    print('%i galaxies'%lIDs)
    mss = np.copy(mstar)
    bhm = np.copy(mbhs)

    ave_lbt_B = np.zeros(lIDs)
    ave_bhar_B = np.zeros(lIDs)
    ave_lbt = np.zeros(lIDs)
    ave_bhar= np.zeros(lIDs)
    
    tot_bhs = np.copy(lIDs)

    '''
    fig,ax = plt.subplots(figsize=(12,6))
    ax.set_xlabel('Lookback Time [Gyr]',fontsize=14)
    ax.set_ylabel('BHAR [Mdot/yr]',fontsize=14)
    ax.set_ylim(-8,5)
    ax2 = ax.twiny()
    ax2.set_xlabel('redshift z',fontsize=14)
    '''

    
    for rg in range(lIDs):
        if bhm[rg]>0.:
            if z_in_file[rg]<redshift:
                redshift_late=z_in_file[rg]
            elif z_in_file[rg]>=redshift:
                redshift_late=redshift
            bh_id = np.where(BHIDs==IDs[rg])[0]
            zs = z[bh_id]
            condi = (zs>=redshift_late)&(zs<redshift_late+0.01)
            BHARs = BHAR[bh_id]
            BHAR_B = BHAR_BONDI[bh_id]
            zs = zs[condi]
            BHARs = BHARs[condi]
            BHAR_B = BHAR_B[condi]
            order = np.argsort(zs)
            '''
            zticks = np.array([0,0.25,0.5])
            zticks = np.append(zticks,np.arange(1,np.ceil(max(zs))+0.2,1))
            notick = np.array([5,6,8])
            for tick in notick:
                b = np.where(zticks==tick)[0]
                zticks = np.delete(zticks,b)

            LBtime_ticks = []
            for tick in range(len(zticks)):
                integral = quad(t_L_int,0,zticks[tick])
                lbt = (integral*th).to('Gyr')
                LBtime_ticks.append(lbt.value)
            lbt_zs = np.asarray(LBtime_ticks)[:,0]
            '''
            zs = zs[order]
            BHARs = BHARs[order]
            BHAR_B = BHAR_B[order]

            lbt_int = np.asarray([quad(t_L_int,0,i)[0] for i in zs])
            lbt = (th*lbt_int).to('Gyr')
            lbt = lbt.value
            
            nonz_bhar = np.isfinite(np.log10(BHARs))
            nonz_bondi = np.isfinite(np.log10(BHAR_B))


            if len(BHARs)>1:
                ave_lbt_B[rg], ave_bhar_B[rg] = runningmedian(lbt,np.log10(BHAR_B+10**-8),0,-8,'red',dt)
                ###ave_lbt_G[rg], ave_bhar_G[rg]= runningmedian(lbt,np.log10(BHAR_G+10**-8),0,-8,'blue')
                ave_lbt[rg], ave_bhar[rg]= runningmedian(lbt,np.log10(BHARs+10**-8),0,-8,'grey',dt)
            elif len(BHARs)==1:
                ave_bhar_B[rg]= BHAR_B[0]
                ave_bhar[rg] = BHARs[0]
            else:
                pass

        else:
            pass
        '''
        ax.set_xlim(-0.05,max(lbt)+0.2)
        ax.set_xticks(np.arange(0,np.ceil(max(lbt))+0.2,2))

        ax2.set_xlim(-0.05,max(lbt)+0.2)
        ax2.set_xticks(lbt_zs)
        ax2.set_xticklabels(zticks)

        #ax.plot(lbt[0:3],np.log10(BHAR_B)[0:3],c='yellow',marker='x',markersize=5,linestyle='None')
        #ax.plot(lbt[0:3],np.log10(BHAR_G)[0:3],c='green',marker='x',markersize=5,linestyle='None')

        ax.annotate('Mbh=%3.2e Msun,Mstar=%3.2e Msun'%(bhm[rg],mss[rg]),xy=(0.05,0.85), xycoords='axes fraction',size=12,bbox=dict(boxstyle="round", fc="w"))
        '''

    #ax.legend(loc=1,prop={'size':12})
    #fig.savefig(imgdir+'bhar_z_RG.png',bbox_inches='tight')

    
    outfile = simdir+'Groups/BH_Bondi_%s_%03i.hdf5'%(box,snapnum)
    g = h5.File(outfile,'a')

    #names = ['AVE_BHAR_10MYR','AVE_BHAR_B_10MYR','AVE_BHAR_20MYR','AVE_BHAR_B_20MYR','AVE_BHAR_30MYR','AVE_BHAR_B_30MYR','AVE_BHAR_50MYR','AVE_BHAR_B_50MYR','AVE_BHAR_70MYR','AVE_BHAR_B_70MYR','AVE_BHAR_100MYR','AVE_BHAR_B_100MYR']#,'AVE_BHAR_GT_100MYR']
    #datas = [ave_bhar[0],ave_bhar_B[0],ave_bhar[1],ave_bhar_B[1],ave_bhar[2],ave_bhar_B[2],ave_bhar[3],ave_bhar_B[3],ave_bhar[4],ave_bhar_B[4],ave_bhar[5],ave_bhar_B[5]]#,ave_bhar_G]
    names = ['BHAR_50MYR','BHAR_BONDI_50MYR']
    datas = [ave_bhar,ave_bhar_B]
    for i in range(len(datas)):
        data = g[names[i]]
        data[...] = datas[i]
        #g.create_dataset(names[i],data=datas[i])
    g.close()  
    
    print('snap %03i complete'%(snapnum))
    
end = time.time()
end = (end-start)/60.

print('Complete in %f minutes'%end)