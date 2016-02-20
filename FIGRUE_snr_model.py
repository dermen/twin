from scipy.interpolate import interp1d
from pylab import *

import pickle

def norm_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma.min()
    cmax = cor_vec_ma.max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n / crange
    

# LOAD THE DETECTED PEAK POSITIONS
cpsi_vals, poss_peaks = load( 'cpsi_vals_peaks.npy')

# LOAD THE MODEL
f = open('pure_fcc_Dec2015_deca_100kAtoms_30kShots_model.pkl','r')
model_params = pickle.load(f)
I = interp1d( **model_params )
x = np.array(model_params['x'])
y = np.array(model_params['y'])
I_mirror = interp1d(x[x<=0]*-1, y[x<=0], bounds_error=0)

C = I(cpsi_vals[ cpsi_vals>=0] )
C_mirror= I_mirror( cpsi_vals[ cpsi_vals >=0])

    #############
C_model = sqrt( norm_corr(C)* norm_corr(C_mirror) )
    #############

# load the gaussian amplitudes, one for each of poss_peaks
pk_amps = np.load('amps_out.npy')
N_boot, N_nshot, N_pk = pk_amps.shape

pk_sig = sqrt(0.00019) # fixed

GAUSS = lambda x,A,SIG,MU : \
       A*exp( -(x-MU)**2 / (2*SIG**2) )

## LET THE SIGNAL BE THE GAUSSIANS EVALUATED
pk_amps_mn = pk_amps.mean(0)
sum_gauss = np.zeros_like(cpsi_vals)
for i_pk in xrange( N_pk):
    pk_amp = pk_amps_mn[-1,i_pk]
    pk_mu = cpsi_vals[poss_peaks[i_pk]]
    pk_gauss = GAUSS(cpsi_vals, pk_amp, 
                    pk_sig, pk_mu)
    sum_gauss += pk_gauss

# MIRROR THAT SHIT
I_mirror_dat = interp1d(-1* cpsi_vals[cpsi_vals<=0],
                        sum_gauss[cpsi_vals<=0],
                        bounds_error=0)
C_data = sum_gauss[ cpsi_vals >=0 ]
C_mirror_data = I_mirror_dat(cpsi_vals[cpsi_vals>=0])

    #############
C_data = sqrt( norm_corr(C_data)*norm_corr(C_mirror_data)  )
    #############


### FANCY PLOTTING BEGINS HERE
shot_series = [500*i for i in xrange(1,200)]
# calc variance in amp measurement
pk_amps_std = pk_amps.std(0)

# load noise and calc variance in noise estimate
noise_data = np.load('sigs_out.npy')
nz = noise_data.mean(0)
nz_std = noise_data.std(0)

# peak at 7/9 
i1_79 = 0
i2_79 = 13
pk_79 = sqrt(pk_amps_mn[:,i1_79] * pk_amps_mn[:,i2_79])
pk_79_nz = sqrt(pk_amps_std[:,i1_79] * pk_amps_std[:,i2_79])

# peak at 1/3
i1_13 = 5
i2_13 = 10
pk_13 = sqrt(pk_amps_mn[:,i1_13] * pk_amps_mn[:,i2_13])
pk_13_nz = sqrt(pk_amps_std[:,i1_13] * pk_amps_std[:,i2_13])

# peak at smaller 0.39
i1_39 = 4
i2_39 = 11
pk_39 = sqrt(pk_amps_mn[:,i1_39] * pk_amps_mn[:,i2_39])
pk_39_nz = sqrt(pk_amps_std[:,i1_39] * pk_amps_std[:,i2_39])

# spacing
import itertools

style.use('ggplot')
color_cycle = itertools.cycle( rcParams['axes.color_cycle'] )
red = color_cycle.next()
blue = color_cycle.next()
purple = color_cycle.next()
black = color_cycle.next()
yellow = color_cycle.next()
green = color_cycle.next()
pink = color_cycle.next()

red = 'r'
blue = 'Darkorange'
green = 'c'

fig = figure(1, figsize=(6,6), dpi=80)

# FIGURE 1
lab13 = r'$\cos\, \psi = 1/3$'
lab79 = r'$\cos\, \psi = 7/9$'
lab39 = r'$\cos\, \psi = 0.39$'

subplot(111)
fill_par = {'alpha':0.4}
plot_par = {'lw':0}



plot( shot_series, pk_13/nz ,
        color=blue,
        **plot_par)
plot( shot_series[1::2], (pk_13/nz) [1::2] ,'o',
        color=blue, ms=5, label=lab13)
err13 = sqrt( (pk_13_nz/pk_13)**2 + (nz_std/nz)**2)*pk_13/nz 
fill_between(x=shot_series,
        y1=pk_13/nz + err13,
        y2=pk_13/nz -err13, 
        color=blue,
        **fill_par ) 

##############


plot( shot_series, pk_79/nz ,
        color=red,
        **plot_par)
plot( shot_series[1::2], (pk_79/nz)[1::2] ,'s',
        color=red, ms=5, label=lab79)
err79 = sqrt( (pk_79_nz/pk_79)**2 + (nz_std/nz)**2)*pk_79/nz 
fill_between(x=shot_series,
        y1=pk_79/nz + err79,
        y2=pk_79/nz -err79, 
        color=red,
        **fill_par )


##############
plot( shot_series, pk_39/nz ,
        color=green,
        **plot_par)
plot( shot_series[1::2], (pk_39/nz)[1::2] ,'^',
        color=green, ms=5, label=lab39)
err39 = sqrt( (pk_39_nz/pk_39)**2 + (nz_std/nz)**2)*pk_39/nz 
fill_between(x=shot_series,
        y1=pk_39/nz + err39,
        y2=pk_39/nz -err39, 
        color=green,
        **fill_par )


plot( shot_series, ones(len(shot_series))*2, '--', lw=2 , color=black)
xlim(500,100000)
ax2 = gca()
ax2.yaxis.tick_left()
ax2.set_yscale('log')
fs = 12
tick_params(axis='both', which='major',length=6,width=1,labelsize=fs)
tick_params(axis='both', which='minor',length=4,width=1,labelsize=fs)
grid(b=1, which='both', color=black, linestyle='--', alpha=0.5, lw=1)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylim(1,20)
ax2.set_xlim(0,100000)
ax2.xaxis.tick_bottom()
ax2.set_ylabel(r'$Z$', fontsize=fs)
ax2.set_xlabel(r'$N_{\mathcal{P}}$', fontsize=fs) #, labelpad=1)
#ax2.yaxis.set_label_position("right")
ax2.set_axis_bgcolor('w')

leg = legend(numpoints=1, loc=2)
fr = leg.get_frame()
fr.set_facecolor('w')
fr.set_alpha(0.5)

######
# NEXT FIGURE
"""
subplot(121)


rng_79 = arange(460,495) #red
rng_13 = arange(180,230) #blue
rng_39 = arange(230,255) #green
xval = cpsi_vals[cpsi_vals >= 0]

bg_line_param = {'lw':5, 'alpha':1}
plot( xval[rng_79], C_data[rng_79], color=red,**bg_line_param)
plot( xval[rng_13], C_data[rng_13], color=blue, **bg_line_param)
plot( xval[rng_39], C_data[rng_39], color=green,**bg_line_param)

plot( xval, C_data, color=black, ls='-', 
        lw=0.5, marker='s', ms=1)
#plot( xval, C_model, lw=2,
#        color=purple, ls='--', marker='s', ms=5)

tick_params(axis='both', which='major',length=6,width=1,labelsize=fs)
ax1 = gca()
ax1.set_xticks( [x/10. for x in xrange(2,9,2)] )
ax1.set_xticklabels( [ '%.1f'%(x/10.) for x in xrange(2,9,2)] )
ax1.xaxis.tick_bottom()
ax1.yaxis.tick_left()
grid(b=1, which='both', color=black, linestyle='--', alpha=0.5, lw=1)

ax1.set_ylabel(r'$G_F\,(\cos(\psi))$',fontsize=fs)
ax1.set_xlabel(r'$|\cos(\psi)|$',fontsize=fs)


subplots_adjust(wpace=.1)
"""
