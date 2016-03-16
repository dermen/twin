import pandas
from pylab import *
import numpy as np
import itertools

#import  postproc_helper as helper

plt.style.use('ggplot')
color_cycle = itertools.cycle( plt.rcParams['axes.color_cycle']  )
colors = [ color_cycle.next() for i in xrange(10) ] 
red = colors[0]
blue = colors[1]
purple = colors[2]
black = colors[3]
yellow = colors[4]
green = colors[5]
pink = colors[6]
orange = 'Darkorange'
###############
# pR histogram
##############

fs=12

df = pandas.read_pickle('8.6keV_shots.pkl')

fig = figure(2, figsize=(6,4),dpi=100)
qval = lambda qinvang:(0.053 / 0.00005)*\
        tan(2*arcsin( qinvang*1.442/4/pi ))
qmin = qval( 2.6 )
hist( qmin + df.q111_rawindex.values, color=blue,
    bins=70, log=1, histtype='stepfilled', 
            lw=0, alpha=0.5)
#hist( qmin + df.q111_rawindex.values, color=blue,
#            bins=70, log=1, histtype='step', lw=3)
xlabel(r'$r^*_i\,(\mathrm{pixels\, units})}$',fontsize=fs, labelpad=9)
ylabel(r'$\mathrm {number\, of\, exposures}$', 
        fontsize=fs)
tick_params(axis='x', which='major',length=6,width=1,
            labelsize=fs)
tick_params(axis='y', which='major',length=0,width=0,
            labelsize=fs)
tick_params(axis='both', which='minor',length=4,width=1, 
            labelsize=fs)
#tick_params(axis='x', which='minor',length=4,width=1, 
#            labelsize=12, pad=10)
iqmin = 748
iqmax = 784
ax = gca()
mid_tick = 0#int(.5*iqmin +.5*iqmax)
ax.set_xticks([ iqmin, 760,772, iqmax] )
ax.set_xticklabels([r'$%d$'%iqmin, 
                    r'$760$', r'$772$',
                    r'$%d$'%iqmax])


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

new_yticks=ax.get_yticks()[1:-1]
ax.set_yticks( new_yticks)
ax.set_yticks([],minor=1)
ax.xaxis.tick_bottom()
ax.yaxis.tick_left()
ax.set_axis_bgcolor('w')
grid(lw=1, color=black, alpha=0.5, ls='--', axis='y')
#plot( ones( 2 )*iqmin, [ 0, 1e5], color=red, lw=2, 
#            ls='--', alpha=0.7)
#plot( ones( 2 )*iqmax, [ 0, 1e5], color=red, lw=2, 
#                    ls='--', alpha=0.7)
ylim(0.2,7e4)
xlim(iqmin-5, iqmax+5)
line_pos = arange(iqmin, iqmax+1,1)
for lp in line_pos:
    plot( ones( 2 )*lp, [ 0, 1e5], color=orange, 
                lw=1, ls='-')
subplots_adjust(left=.14, bottom=.14, right=.95,top=.97)

#############
# angular anisotropies #
##############
shot = np.load('polar_shot.npy')
pmask = np.load('polar_mask.npy')
shot_ma = np.ma.masked_equal(shot*pmask,0)
rad_pro = shot_ma.mean(1)



n,s,m = postproc_helper.get_ring( pdata=shot, pmask=pmask, 
        iq=rad_pro.argmax(), rm_peaks=0)

shot_norm, shot_ring, shot_mask = \
        postproc_helper.get_ring(pdata=shot, pmask=pmask, 
                iq=rad_pro.argmax() ) 
ofit_tmp, fit_tmp, vals_tmp = \
        postproc_helper.fit_periodic(shot_ring.copy(), 
            shot_mask.copy(), deg=15)

shot_ring2, shot_mask2,removed = \
        postproc_helper.remove_peaks( shot_ring.copy(),
        shot_mask.copy(), thick=10,
        coef=ofit_tmp, peak_thresh=2)

ofit_sub,fit_sub,vals_sub = \
    postproc_helper.fit_periodic(shot_ring2.copy(), 
    shot_mask2.copy(), deg=15)

norm = shot_ring2[ shot_mask2 > 0 ].mean()
ofit,fit,vals =\
    postproc_helper.fit_periodic(shot_ring2.copy()/norm, 
    shot_mask.copy(), deg=15)

y_tmp = np.ma.masked_equal(vals_tmp*shot_mask,0)
y = np.ma.masked_equal(vals_sub*shot_mask,0)

x = arange( len(y))

vals_tmp = np.polynomial.chebyshev.chebval( x, ofit_tmp)
vals_sub = np.polynomial.chebyshev.chebval( x, ofit_sub)

degs = x * 360 / len(x)  

fig = figure(1, figsize=(6,4))

semilogy( ma.masked_equal( s*m,0) , 'x', color=blue)
semilogy( ma.masked_equal( shot_ring2*shot_mask2,0) 
        , '^', color=red, ms=5, mec='r')
semilogy( vals_tmp, color=yellow, lw=7.5 )
semilogy( vals_sub,color='k', lw=10.5 , ls='--')  

tick_lab = [ '%.0f'%(tick*360 / len(y)) \
        for tick in gca().get_xticks() ]
gca().set_xticklabels(tick_lab)

tick_params(axis='both', which='major',
        length=0, width=1, labelsize=12)
tick_params(axis='both', which='minor',
        length=0, width=1, labelsize=12)
ylabel(r'$I_i(\phi)\,(\mathrm{counts})$', fontsize=12)#, labelpad=10)
xlabel(r'$\phi\,(0-2\pi)$', fontsize=12)#, labelpad=10)
xlim(-250,5250)
ylim(260,30944)
ax = gca()
grid(1, which='major', axis='both', lw=1, ls='--', 
        color='#777777')
ax.xaxis.tick_bottom()
ax.yaxis.tick_left()
ax.set_axis_bgcolor('w')

ax.set_xticks([0,int(5000/4.),2500,int(5000 * 3./4),5000])
ax.set_xticklabels([r'$0$',r'$\pi/2$',r'$\pi$',
        r'$3\pi/2$', r'$2\pi$'])

subplots_adjust(left=.15,bottom=.15)
savefig('Fig_angprof_supp.pdf')

