from pylab import *
from scipy.interpolate import interp1d

import pickle

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

def norm_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma[:-100].min()
    cmax = cor_vec_ma[:-100].max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n / crange
    
# LOAD THE DETECTED PEAK POSITIONS
cpsi_vals, poss_peaks = load( 'cpsi_vals_peaks.npy')

# LOAD THE MODEL
#f = open('pure_fcc_Dec2015_deca_100kAtoms_30kShots_model.pkl','r')
#model_params = pickle.load(f)
x,y = model = np.load( 'model.npy')

I = interp1d( x,y,bounds_error=0) #**model_params )

cpsi_model_full = x
C_model_full = norm_corr( y)
#x = np.array(model_params['x'])
#y = np.array(model_params['y'])

C = I(cpsi_vals )

C_model = norm_corr(C)

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

C_data = norm_corr( sum_gauss)

# load the correlation datas
cpsi, cors_m = np.load( 'all_cors_m2.npy')
_, cors, difcors = np.load('cor_vs_difcor.npy')

difcors = cors_m

Icors = interp1d( cpsi, cors, bounds_error=0)
Idifcors = interp1d( cpsi, difcors, bounds_error=0)

#xdata = np.linspace( cpsi.min(), 0.99, 2500)
xdata = np.linspace( cpsi.min(), 0.8, 2500)
Pcors = np.polyfit( xdata, Icors(xdata), deg=6 )
Pdifcors = np.polyfit( xdata, Idifcors(xdata), deg=6 )

cors_flat = Icors(xdata) - polyval(Pcors, xdata )
difcors_flat = Idifcors(xdata) #- polyval(Pdifcors, xdata )

cors_flat = norm_corr( cors_flat)
difcors_flat = norm_corr( difcors_flat)

from mpl_toolkits.axes_grid.inset_locator import inset_axes

cpsi_min = xdata.min()
def plot_lines():
    line_par = {'lw':2, 'ls':'-', 'color':red , 'alpha':0.6}
    plot( ones(2)*7/9., [-10,10], **line_par )
    plot( ones(2)*1/3., [-10,10], **line_par )
    plot( ones(2)*5/9., [-10,10], **line_par ) 
    plot( ones(2)*-7/9., [-10,10], **line_par )
    plot( ones(2)*-1/3., [-10,10], **line_par )
    plot( ones(2)*-5/9., [-10,10], **line_par )
    plot( ones(2)*cpsi_min, [-10,10], lw=3, 
                alpha=0.25, color='k', ls='--')
    plot( ones(2)*-1*cpsi_min, [-10,10], lw=3, 
                alpha=0.25, color='k', ls='--')

############
# NEW FIG
par = {'lw':3 }
fs = 12
#######################
# TOP TICKS
min_ticks_top =[-1, -7/9., -5/9., 
        -1/3., 1/3., 5/9., 7/9. ,1] 
maj_ticks_top =[ cpsi_min ,0, -cpsi_min ] 
min_tlab_top = ['' for x in min_ticks_top]
maj_tlab_top = [  r'$\cos\, \psi_{\max}$',r'$\cos\, \psi = \pi/2$',  r'$-\cos\, \psi_{\max}$' ]
####################
# BOTTOM TICKS
maj_ticks_bot = [-1,  -7/9., -5/9.,
        -1/3., 1/3., 5/9., 7/9. ,1] 
maj_tlab_bot = [ r'$-1$', r'$-7/9$',  r'$-5/9$', r'$-1/3$', 
          r'$1/3$', r'$5/9$', r'$7/9$', r'$1$']
min_ticks_bot =[ cpsi_min ,0,  -cpsi_min ] 
min_tlab_bot =[ '', '', ''] 


tx = 0.65
ty = 0.8
fig = figure(1, figsize=(6.2,5), dpi=100)
##########################################
subplot(311)
#plot_lines()
gca().set_yticks([])
# major xticks
tick_params(axis='both',labelbottom='off',
            labeltop='on',bottom='off', which='major',
            length=6,width=1,labelsize=fs) #, pad=8)
# minor ticks
tick_params(axis='both',labelbottom='off',
            labeltop='on',bottom='off', which='minor',
            length=0,width=0,labelsize=fs) #, pad=8)

gca().set_xticks(maj_ticks_top)
gca().set_xticklabels(maj_tlab_top)
gca().set_xticks(min_ticks_top, minor=1)
gca().set_xticklabels(min_tlab_top, minor=1)
gca().xaxis.set_tick_params(labeltop='on')
grid(1, ls='--', lw=1, alpha=0.5, color=black, which='both')
lab0=r'$C_{\mathrm{deca}}$'#(\cos\,(\psi))$'
plot0 = plot( cpsi_model_full, C_model_full,color='Limegreen', label=lab0, **par)
plot( cpsi_model_full, C_model_full,color=black, marker='s', mec=black, ms=1, lw=0)
ylim(-0.1,1.3)
xlim(-.95,.95)
ax = gca()
ax.set_axis_bgcolor('white')
#leg =legend(loc=9)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_facecolor('white')
#fr.set_alpha(0.5)
ax1 = gca()


ax1.text(tx, ty, s='A', fontsize=12, color='r')
ax1.spines['bottom'].set_visible(0)

#####################################
subplot(312)
lab1 = r'$D - P_D$'
lab2 = r'$P_D$'
#lab2 = r'$C - P_C$'
plot1 = plot(xdata, difcors_flat, label=lab1, color=blue , **par)
#plot11 = plot(xdata, norm_corr( polyval(Pdifcors, xdata )), label=lab2, color=black) 
#plot2 = plot( xdata[::20], cors_flat[::20], 'd',
#        label=lab2,  ms=3 ,color=black)#, alpha=0.4)
#plot_lines()
#gca().set_yticks([])
#gca().set_xticks([])
#grid(0)
ylabel(r'counts$\,^2$ A.U.',fontsize=fs)
#leg = legend(loc=(0.37,0.56) , # , prop={'size':6},
       # fontsize=fs, numpoints=1,
      # borderpad=0.57)

ylim(-0.1, 1.3)
#leg =legend(loc=9)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_linewidth(0.0)
#fr.set_facecolor('white')
#fr.set_alpha(0.8)


ax2 = gca()
ax2.set_xticks(min_ticks_bot + min_ticks_top)
ax2.set_yticks([])
#ax2.get_xaxis().set_visible(0)
#ax2.get_yaxis().set_visible(0)
ax2.tick_params(axis='both',top='off', which='both',
            length=0,width=0,labelsize=0,pad=8)
ax2.grid(1, ls='--', lw=1, alpha=0.5, color=black, which='both')
ax2.spines['top'].set_visible(0)
ax2.set_axis_bgcolor('white')
ax2.grid(1, ls='--', lw=1, alpha=0.5, color=black, which='both')

ax2.set_xlim(-.95,.95)
ax2.text(tx, ty, s='B', fontsize=12, color='r')

####################################
subplot(313)
#plot_lines()
lab3 = r'$G$'
plot3 = plot( cpsi_vals, C_data, color='#CD00CD',alpha=0.5,  label=lab3, **par)

ylim(-0.1,1.3)
gca().set_yticks([])
tick_params(axis='both',top='off', which='major',
            length=6,width=1,labelsize=fs,pad=8)
tick_params(axis='both',top='off', which='minor',
            length=0,width=0,labelsize=fs)

gca().set_xticks(min_ticks_bot, minor=1)
gca().set_xticklabels(min_tlab_bot, minor=1)
gca().set_xticks(maj_ticks_bot)
gca().set_xticklabels(maj_tlab_bot) 
xlabel(r'$\cos\, \psi $', fontsize=fs)
#grid(0)#, lw=2, color=black)
xlim(-.95,.95)

grid(1, ls='--', lw=1, alpha=0.5, color=black, which='both')
ax3 = gca()
ax3.set_axis_bgcolor('white')
ax3.spines['top'].set_visible(0)


#leg =legend(loc=9)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_facecolor('white')
#fr.set_alpha(0.5)

#leg = figlegend(handles = plots, labels=labels, loc=9) 
        #(0.37,0.56) , # , prop={'size':6},
       #fontsize=fs, numpoints=1,
      # borderpad=0.57)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_facecolor('white')
ax3.text(tx, ty, s='C', fontsize=12, color='r')

subplots_adjust(bottom=.18, hspace=0)
#savefig('/home/mender/FIGURE2_results.png', dpi=300)

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
#xlabel(r'$\cos\,(\psi)$',fontsize=18)


