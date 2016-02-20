from pylab import *
from scipy.interpolate import interp1d
import pickle


import os
os.chdir('/data/work/mender/loki')

import postproc_helper as helper

#from postproc_helper import is_outlier, remove_peaks
os.chdir('/data/work/mender/FIGURE_snr_model')


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

fs = 12
# ring with peaks shaved off
d = np.load ('peaks_removed.npy')
coef = np.load('coef.npy')
mask = zeros_like(d)
mask[ d!=0] = 1
dm, mask, d_rm = helper.remove_peaks( d, mask, coef=coef, peak_thresh=2.5)

#w = is_outlier( d, 3 )

w = d_rm > median(d)

dm = ma.masked_equal( dm, 0 )

#dm = ma.masked_where( logical_or( w, d==0), d)
d_rm = ma.masked_where( ~w, d )

# NP sizes
s = np.load('NP_sizes_178802.npy')
s = s[ s < 500] 
s = s[ s > 10 ] 

a,b = np.histogram( s, bins=100 )
bins = [b[i]/2. + b[i+1]/2. for i in xrange(len(b)-1 ) ]



############
# NEW FIG
import matplotlib.gridspec as gridspec
fig = figure(1, figsize=(6.2,5), dpi=200) #, dpi=400)

gs = gridspec.GridSpec(2,3)
ax1 = subplot( gs[0,:-1])
ax2 = subplot(gs[1,:-1])

ax1.plot(d_rm, 's', color='Limegreen', ms=6, mec='Limegreen')
ax1.plot( d_rm, 'd',color=black,  ms=2)
ax1.plot( dm, 'd',color=black,  ms=2)
ax1.set_ylim(0,20000)
ax1.set_ylabel(r'$I_i(\phi)$ (counts)', fontsize=fs)

ax2.plot( dm, 'd', color=black, ms=3)
ax1.tick_params(axis='both', which='major',labeltop='on', 
        labelbottom='off', length=6, width=1, labelsize=9)
#ax1.set_xticklabels([])
ax1.set_xticks([2500/2, 2500*3/2])
ax1.set_xticklabels([r'$\pi/2$', r'$3\pi/2$'])

ax1.yaxis.tick_left()
ax1.xaxis.tick_top()
ax1.set_yticks([ x*1000 for x in xrange(2,17,2) ])
ax1.set_yticklabels([ '%dk'%x for x in xrange(2,17,2) ])
ax1.text(200, 14500, s='A', fontsize=12, color='r')

ax2.tick_params(axis='both', which='major', length=6, width=1, labelsize=9)
ax2.set_ylim(0,2700)
ax1.set_ylim(0,17500)

ax1.grid(1, ls='--', color='k', alpha=0.5) 
ax2.grid(1, ls='--', color='k', alpha=0.5) 
ax2.set_xlabel(r'$\phi \,(0-2\pi)$', fontsize=fs)
ax2.set_ylabel(r'$I_i(\phi)$ (counts)', fontsize=fs)

ax2.set_xticks([2500/2, 2500*3/2])
ax2.set_xticklabels([r'$\pi/2$', r'$3\pi/2$'])
ax2.yaxis.tick_left()
ax2.xaxis.tick_bottom()
ax2.text(200, 2200, s='B', fontsize=12, color='r')


ax2.set_yticks([ x*500 for x in xrange(1,6) ])
ax2.set_yticklabels([ '%.1fk'%(x/2.) for x in xrange(1,6) ])

ax3 = subplot(gs[:,-1])
ax3.grid(1, ls='--', color=black, alpha=0.3, 
        which='both', lw=1)

ax3.bar( bins[::2], 8*a[::2]/5000. , width=8, lw=1, color='Limegreen', 
            edgecolor=black)
ax3.text(170, 22, s='C', fontsize=12, color='r')
#bins = array( bins)
#a = array(a)
#bins = array(zip( bins, bins+.5))
#a = array(zip( a,a))
#ax3.plot( bins, 8*a/5000. , lw=2, color='Limegreen')

#ax3.hist( s , bins=100,normed=1, histtype='step', lw=2, color=black )
#ax3.hist( s , bins=100, histtype='stepfilled', 
#        color='Limegreen', lw=0, alpha=0.5 )

#ax3.set_xscale('log')
ax3.set_xlim(10,200)
ax3.set_xticks([50,100,150])

ax3.set_yticks([5,10,15,20])

ax3.tick_params(axis='x',right='on', which='major',
            length=6,width=1,labelsize=9)
ax3.yaxis.tick_right()
ax3.tick_params(axis='both', which='major', length=6, width=1)
ax3.tick_params(axis='both', which='minor', length=4, width=1)
ax3.yaxis.set_label_position("right")
ax3.set_ylabel(r'$N_\mathrm{spot}(d_L)$ per snapshot', fontsize=fs)
ax3.set_xlabel(r'$d_L$ (nm)', fontsize=fs)
ax3.xaxis.tick_bottom()

ax1.set_axis_bgcolor('w')
ax2.set_axis_bgcolor('w')
ax3.set_axis_bgcolor('w')


subplots_adjust(top=.9,bottom=.13, left=.15, right=.86, hspace=.04, wspace=.07 )
#savefig('size_estimate.png', dpi=300)


#estimate NPs in beam
def wd (d,vol='sphere', a=.4076 ):
    wg = 196967/6.022e23 # weight of gold atom
    if vol=='sphere':
        V = 4/3. * pi * (d/2.)**3
    V_unit = a**3
    return  wg*4* V / V_unit  
    
wL = sum( [ wd( bins[i] )*a[i] / 5000. / .8  
        for i in xrange(len(bins))  ])

c = 40 # mg/mL
Vml = 4.68e-10
wS = c*Vml - wL


"""
##########################################
lab1 = r'$\tilde {C}$'
lab2 = r'$C$'
plot_lines()
plot(xdata, difcors_flat, label=lab1, **par)
plot( xdata[::20], cors_flat[::20], 'd',
        label=lab2,  ms=2 ,color='k')#, alpha=0.4)
gca().set_yticks([])
# major xticks
tick_params(axis='both',labelbottom='off',
            labeltop='on',bottom='off', which='major',
            length=6,width=1.5,labelsize=9) #, pad=8)
# minor ticks
tick_params(axis='both',labelbottom='off',
            labeltop='on',bottom='off', which='minor',
            length=4,width=1.5,labelsize=9) #, pad=8)
gca().set_xticks(maj_ticks_top)
gca().set_xticklabels(maj_tlab_top)
gca().set_xticks(min_ticks_top, minor=1)
gca().set_xticklabels(min_tlab_top, minor=1)
gca().xaxis.set_tick_params(labeltop='on')
ylim(-0.1, 1.3)
grid(0)
leg = legend(loc=(0.37,0.56) , prop={'size':6},
        fontsize=9, numpoints=1,
        borderpad=0.57)
fr = leg.get_frame()
fr.set_edgecolor('black')

#####################################
subplot(312)
plot_lines()
plot( cpsi_model_full, C_model_full,**par)
ylim(-0.1,1.3)
gca().set_yticks([])
gca().set_xticks([])
grid(0)
ylabel(r'(photon counts)$\,^2$ (arb. scale)',fontsize=11)

####################################
subplot(313)
plot_lines()
plot( cpsi_vals, C_data, **par)
ylim(-0.1,1.3)
gca().set_yticks([])
tick_params(axis='both',top='off', which='major',
            length=6,width=1.5,labelsize=9)
tick_params(axis='both',top='off', which='minor',
            length=4,width=1.5,labelsize=9)

gca().set_xticks(min_ticks_bot, minor=1)
gca().set_xticklabels(min_tlab_bot, minor=1)
gca().set_xticks(maj_ticks_bot)
gca().set_xticklabels(maj_tlab_bot) 
xlabel(r'$\cos(\psi)$', fontsize=12)
grid(0)


subplots_adjust(bottom=.18, hspace=0.)
savefig('/home/mender/FIGURE2_results.png', dpi=300)

########
# OLD FIG
#subplot(311)
#plot_lines()
#plot(cpsi, cors, color=purple, lw=6, ls=':')
#plot(cpsi, cors, color=purple, lw=2)
#plot(cpsi, difcors, color='k', lw=2)
#ylim(-0.008, 0.008)
#xlim(-1,1)
#gca().set_yticks([])
#gca().set_xticks([])

#subplot(312)
#plot_lines()
#plot( cpsi, cors_m, color='k', lw=2)
#ylim(-3.72,1.76)
#xlim(-1,1)
#gca().set_yticks([])
#gca().set_xticks([])
#ylabel('CXS signal A.U.', fontsize=14)

#subplot(313)
#plot_lines()
#plot( cpsi_vals, C_model, color=green,lw=8, ls=':' )
#plot( cpsi_vals, C_model, color=green,lw=2 )
#plot( cpsi_vals, C_data, color='k', lw=2,alpha=0.5)
#xlim(-1,1)
#ylim(0,1.2)
#gca().set_yticks([])
#xlabel(r'$\cos(\psi)$',fontsize=18)

"""
