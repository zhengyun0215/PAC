from mpi4py import MPI
from astropy.table import Table
from sklearn.neighbors import BallTree
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import numpy as np
cosmo = FlatLambdaCDM(H0 = 67.74,Om0 = 0.3089,Ob0 = 0.0486)


mh = [11.0,10.5,10.0,9.5]
ml = [10.5,10.0,9.5,9.0]
#ml0 = np.linspace(9.5,10.5,6)
ms = [11.3,11.5,11.7]
z = [0.025,0.075,0.125,0.175]
zs = ['0p025','0p075','0p125','0p175']
#field_len = 10
photo_len = len(ml)
spec_len = len(ms)
z_len = len(z)
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


spec = rank // (z_len*photo_len)
photo = (rank % (z_len*photo_len)) // z_len
k = (rank % (z_len*photo_len)) % z_len



t_spec = Table.read('./spec.fits')
t_spec = t_spec[(t_spec['mass']>=ms[spec])&(t_spec['mass']<ms[spec]+0.2)]
t_spec = t_spec[(t_spec['Z']>=z[k]-0.025)&(t_spec['Z']<=z[k]+0.025)]


if len(t_spec)>0:
    print("Let's go on!!!!")

    t_photo = Table.read('./photo.fits')
    
    ############color cut###########
    t_photo = t_photo[(t_photo['u_r'+zs[k]]>(0.11*t_photo['mass'+zs[k]] - 0.45 * 0.1 + 0.94))]
    t_photo = t_photo[(t_photo['mass'+zs[k]]>=ml[photo])&(t_photo['mass'+zs[k]]<mh[photo])]

    t_photo = t_photo[('ra','dec')]
    t_ran1 = Table.read('./random1.fits')
    t_ran2 = Table.read('./random2.fits')


    print("The spec is:")
    print(spec)
    print("The photo is:")
    print(photo)
    print("The k is:")
    print(k)
    print("The massbin for central is")
    print(ms[spec])
    print("The redshift bin is:")
    print(z[k]-0.025,z[k]+0.025)
    print("The massbin for sat is:")
    print(ml[photo],mh[photo])
    print("The spec length:")
    print(len(t_spec))
    print("The photo catalogue length:")
    print(len(t_photo))

    
    tree_ran2 = BallTree(np.array([t_ran2['ra']*np.pi/180, t_ran2['dec']*np.pi/180]).T, leaf_size=16, metric='haversine')
    tree_photo = BallTree(np.array([t_photo['ra']*np.pi/180, t_photo['dec']*np.pi/180]).T, leaf_size=16, metric='haversine')

    r1 = cosmo.comoving_distance(z[k]).value*cosmo.h
    r = np.logspace(-2, 2, 21)
    r = r[4:]
    r_need = 10**((np.log10(r[1:])+np.log10(r[:-1]))/2.0)
    CD = cosmo.comoving_distance(t_spec['Z']).value*cosmo.h
    center_coor = np.array([t_spec['ra']*np.pi/180, t_spec['dec']*np.pi/180]).T
    center_ran = np.array([t_ran1['ra']*np.pi/180, t_ran1['dec']*np.pi/180]).T
    N_ratio1 = len(t_ran1)/len(t_spec)
    N_ratio2 = len(t_ran2)/len(t_photo)
    CD_ran = np.random.choice(CD,len(t_ran1),replace=True)
    D1R2 = np.array([tree_ran2.query_radius(center_coor, r[a]/CD, count_only=True) for a in range(len(r))])/N_ratio2
    D1D2 = np.array([tree_photo.query_radius(center_coor, r[a]/CD, count_only=True) for a in range(len(r))])
    D2R1 = np.array([tree_photo.query_radius(center_ran, r[a]/CD_ran, count_only=True) for a in range(len(r))])/N_ratio1
    R1R2 = np.array([tree_ran2.query_radius(center_ran, r[a]/CD_ran, count_only=True) for a in range(len(r))])/N_ratio1/N_ratio2
    D1R2 = D1R2.T
    D1D2 = D1D2.T
    D2R1 = D2R1.T
    R1R2 = R1R2.T
    index_spec = np.arange(len(t_spec))
    index_ran = np.arange(len(t_ran1))
    jack=50
    index_spec = np.array_split(index_spec,jack)
    index_ran = np.array_split(index_ran,jack)
    cl = np.zeros((jack,len(r_need)))
    area = len(t_ran2)/2500.*(np.pi/180)**2
    for i in range(jack):
        D1R2_jack = np.delete(D1R2,index_spec[i],0)
        D1R2_mean = np.sum(D1R2_jack, 0)
        D1R2_mean = D1R2_mean[1:]-D1R2_mean[:-1]
        D1D2_jack = np.delete(D1D2,index_spec[i],0)
        D1D2_mean = np.sum(D1D2_jack, 0)
        D1D2_mean = D1D2_mean[1:]-D1D2_mean[:-1]
        D2R1_jack = np.delete(D2R1,index_ran[i],0)
        D2R1_mean = np.sum(D2R1_jack, 0)
        D2R1_mean = D2R1_mean[1:]-D2R1_mean[:-1]
        R1R2_jack = np.delete(R1R2,index_ran[i],0)
        R1R2_mean = np.sum(R1R2_jack, 0)
        R1R2_mean = R1R2_mean[1:]-R1R2_mean[:-1]
        cl[i] = len(t_photo)/area/r1**2*((D1D2_mean-D1R2_mean-D2R1_mean)/R1R2_mean+1)
    np.save('./color_cut_red/%s_%.1f_%.1f_results_jack.npy' % (zs[k],ms[spec]+0.1,ml[photo]),cl)
    np.save('./color_cut_red/%s_%.1f_%.1f_param_jack.npy' % (zs[k],ms[spec]+0.1,ml[photo]),np.array([len(t_spec)]))
    
    
print("finish!")
