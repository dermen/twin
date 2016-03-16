from pylab import *
from scipy.interpolate import interp1d

import pickle

# spacing
import itertools

style.use('ggplot')
color_cycle = itertools.cycle( \
            rcParams['axes.color_cycle'] )
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

def norm_corr_( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma[:-100].min()
    cmax = cor_vec_ma[:-100].max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n / crange, crange


def positive_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma[:-100].min()
    cmax = cor_vec_ma[:-100].max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n 


def friedel_corr(cx,cy ):
    I=interp1d(cx,cy)
    cx_sym = linspace( cx.min(), -cx.min(), len(cx) )
    cy_sym = np.zeros_like( cx_sym )
    for i,cospsi in enumerate(cx_sym):
        cy_sym[i] =  .5*I(cospsi) + .5*I(-cospsi)
    return cx_sym, cy_sym

def make_friedel(cx,cy ):
    I=interp1d(cx,cy)
    I2=interp1d(-cx,cy)
    c_lower = .5*(I( cx[ cx <0]) + I(cx[ cx <0]))
    c_upper = .5*(I2(cx[cx >= 0]) + I2(cx[ cx >=0]))

    return np.concatenate((c_lower, c_upper))

xx, difcors = np.load( 'sym-mean-corsFeb29-0pix.npy')
difcors_flat = norm_corr(difcors)
xx, C_datafit = np.load( 'data_fit.npy')
C_datafit, norm = norm_corr_( C_datafit)
intercor_full = np.load('inter_mean_full_87943.npy')
intercor_full/= norm
noise =  intercor_full.std()


cpsi_vals = xx
xdata = xx
# LOAD THE DETECTED PEAK POSITIONS
#cpsi_vals, poss_peaks = load( 'cpsi_vals_peaks.npy')

cpsi_vals = xx

# LOAD THE MODEL
#f = open('pure_fcc_Dec2015_deca_100kAtoms_30kShots_model.pkl','r')
#model_params = pickle.load(f)
x,y = model = np.load( 'model.npy')
#I = interp1d( x,y)#,bounds_error=0) #**model_params )
cpsi_model_full = x
C_model_full = norm_corr( y)
where_x = logical_and( x > cpsi_vals.min(), x < -cpsi_vals.min() )
cpsi_model = x[where_x]
y_inds = np.searchsorted( x, cpsi_model)
C = y[y_inds] #I(cpsi_model )
C_model = norm_corr(C)
#C_model = make_friedel(cpsi_model, norm_corr(C))

#=========================================================
# load the gaussian amplitudes, one for each of poss_peaks
#=========================================================
#pk_amps = np.load('amps_out.npy')
#N_boot, N_nshot, N_pk = pk_amps.shape

#pk_sig = sqrt(0.00019) # fixed

#GAUSS = lambda x,A,SIG,MU : \
#       A*exp( -(x-MU)**2 / (2*SIG**2) )

## LET THE SIGNAL BE THE GAUSSIANS EVALUATED
#pk_amps_mn = pk_amps.mean(0)
#sum_gauss = np.zeros_like(cpsi_vals)
#for i_pk in xrange( N_pk):
#    pk_amp = pk_amps_mn[-1,i_pk]
#    pk_mu = cpsi_vals[poss_peaks[i_pk]]
#    pk_gauss = GAUSS(cpsi_vals, pk_amp, 
#                    pk_sig, pk_mu)
#    sum_gauss += pk_gauss

#C_data = np.load( 'data_fit.npy')
#C_data = norm_corr( sum_gauss)

# load the correlation datas
#cpsi, cors_m = np.load( 'all_cors_m2.npy')
#_, cors, difcors = np.load('cor_vs_difcor.npy')
#difcors = cors_m

#cpsi, difcors = np.load( 'sym-mean-corsFeb29-0pix.npy')
#difcors_flat = norm_corr( difcors)
#Icors = interp1d( cpsi, cors, bounds_error=0)
#xdata = np.linspace( cpsi.min(), 0.99, 2500)
#xdata = np.linspace( cpsi.min(), -cpsi.min(), 2500)
#Pcors = np.polyfit( xdata, Icors(xdata), deg=6 )
#Pdifcors = np.polyfit( xdata, Idifcors(xdata), deg=6 )
#cors_flat = Icors(xdata) - polyval(Pcors, xdata )
#difcors_flat = Idifcors(xdata) #- polyval(Pdifcors, xdata )
#cors_flat = norm_corr( cors_flat)
#difcors_flat = norm_corr( difcors_flat)


cpsi_min = xdata.min()

fs = 12

top_ticklabels = [r'$\cos\, \psi_{\max}$',
        r'$\cos\, \psi = \pi/2$',
        r'$-\cos\, \psi_{\max}$' ]
top_ticks = [ cpsi_min ,0,  -cpsi_min ] 

major_ticklabels = [ r'$-7/9$',  
            r'$-5/9$', r'$-1/3$',
            r'$1/3$', r'$5/9$', 
            r'$7/9$']
major_ticks =[ -7/9., 
            -5/9., -1/3.,
            1/3.,  5/9.,
            7/9.] 
  
minor_ticks = [-.39, -0.037,0.037,.39]

tx = 0.65
ty = 0.8

fig = figure(1, figsize=(6,5), dpi=100)
fs = 12
ax = gca()
ax.set_xlim(-.95,.95)
ax.set_axis_bgcolor('white')
ax.set_ylabel(r'$\mathrm{counts}\,^2\, \mathrm{A.U.}$',
            fontsize=fs)
ax.set_xlabel(r'$\cos \, \psi$', fontsize=fs,
                 labelpad=-7)
ax.set_xticks( major_ticks )
ax.set_xticks(minor_ticks, minor=1)
ax.set_xticklabels( major_ticklabels)
ax.tick_params(which='major', length=7, 
            labelsize=fs,pad=0)
ax.tick_params(which='minor', length=0, 
            labelsize=fs,pad=0)

ax_ = ax.twiny()
ax_.set_xlim(ax.get_xlim())
ax_.set_xticks( top_ticks )# minor=1)
ax_.set_xticklabels( top_ticklabels )# minor=1)
ax_.xaxis.tick_top()
ax_.tick_params(which='major', length=7, 
            labelsize=fs,pad=-3)
ax_.set_yticks([])
ax.set_yticks([])
#ax_.set_ylim(-0.1,3.45)

ax.grid(lw=1, color=red, ls='--', alpha=0.5, which='both')
ax_.grid(lw=2, color='#777777', ls='-', alpha=0.5)

off_model = 2.4
off_data = 1.2
##################
# PLOT THE MODEL #
##################
cpsi_model, C_model = friedel_corr( cpsi_model_full,
                                    C_model_full)
ax.plot( cpsi_model, 
            C_model+off_model,color='Limegreen', 
            lw=3)
ax.text(tx, ty+off_model, s='a', fontsize=fs, color=red)
#################
# PLOT THE DATA #
#################
ax.plot(xdata, 
        difcors_flat+off_data, 
        color=blue , lw=1)
ax.text(tx, ty+off_data, s='b', fontsize=fs, color='r')
################
# PLOT THE FIT #
################
cpsi_fit , C_fit = friedel_corr( cpsi_vals, C_datafit )
ax.plot( cpsi_fit, C_fit, 
            color='#CD00CD',alpha=0.5,  
            lw=3)
ax.text(tx, ty, s='c', fontsize=fs, color='r')

ax.plot( [-.88,.88], ones(2)*noise*2.5, lw=2, 
        alpha=0.7, color=yellow  )
ax.text(0.9, noise*1.5, r'$Z=2.5$', fontsize=12 , 
            color='#777777')
ax.fill_between( cpsi_fit, zeros_like( cpsi_fit), 
        ones_like(cpsi_fit)*noise*2.5, 
        color=yellow, alpha=0.3)

ax.set_ylim(-0.1,3.45)
ax_.tick_params(pad=-8)

#intercor_full = np.load('inter_mean_full_87943.npy')
#phis = arange( 5000 ) * 2 * pi / 5000.
#th = arcsin( 2.668 * 1.442 / 4 / pi )
#cpsi = cos(phis) * cos(th)**2 + sin(th)**2

#noise = intercor_full[ logical_and( 
#            cpsi >cpsi_vals.min(),
#            cpsi < -cpsi_vals.min() )]

#cpsi_ = cpsi[ logical_and( 
#            cpsi >cpsi_vals.min(),
#            cpsi < -cpsi_vals.min())]

#ax.plot( cpsi_, noise, 's',ms=2 , 
#    color='#89E894', mew=0)




##########################################
# major xticks
#tick_params(axis='both',labelbottom='off',
#            labeltop='on',bottom='off', which='major',
#            length=6,width=1,labelsize=fs) #, pad=8)
# minor ticks
#tick_params(axis='both',labelbottom='off',
#            labeltop='on',bottom='off', which='minor',
#            length=0,width=0,labelsize=fs) #, pad=8)

#gca().set_xticks(maj_ticks_top)
#gca().set_xticklabels(maj_tlab_top)
#gca().set_xticks(min_ticks_top, minor=1)
#gca().set_xticklabels(min_tlab_top, minor=1)
#gca().xaxis.set_tick_params(labeltop='on')
#grid(1, ls='--', lw=1, alpha=0.5, color=black, which='major')
#grid(1, ls='-', lw=1, alpha=0.5, color=red, which='minor')
#lab0=r'$C_{\mathrm{deca}}$'#(\cos\,(\psi))$'
#plot0 = plot( cpsi_vals, C_model,color='Limegreen', label=lab0, **par)
#plot( cpsi_vals, C_model,color=black, marker='s', mec=black, ms=1, lw=0)
#ylim(-0.1,1.3)
#xlim(0,.95)
#leg =legend(loc=9)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_facecolor('white')
#fr.set_alpha(0.5)


#ax1.spines['bottom'].set_visible(0)

#####################################
#subplot(312)
#lab1 = r'$D - P_D$'
#lab2 = r'$P_D$'
#lab2 = r'$C - P_C$'
#plot1 = plot(xdata, difcors_flat,  color=black , marker='s', 
#                    mec=black, ms=1, lw=0)
#plot11 = plot(xdata, norm_corr( polyval(Pdifcors, xdata )), label=lab2, color=black) 
#plot2 = plot( xdata[::20], cors_flat[::20], 'd',
#        label=lab2,  ms=3 ,color=black)#, alpha=0.4)
#plot_lines()
#gca().set_yticks([])
#gca().set_xticks([])
#grid(0)
#leg = legend(loc=(0.37,0.56) , # , prop={'size':6},
       # fontsize=fs, numpoints=1,
      # borderpad=0.57)

#ylim(-0.1, 1.3)
#leg =legend(loc=9)
#fr = leg.get_frame()
#fr.set_edgecolor('black')
#fr.set_linewidth(0.0)
#fr.set_facecolor('white')
#fr.set_alpha(0.8)

#ax2 = gca()
#ax2.set_xticks(min_ticks_bot, minor=1)
#ax2.set_xticks(min_ticks_top)

#ax2.set_yticks([])

#ax2.get_xaxis().set_visible(0)
#ax2.get_yaxis().set_visible(0)
#ax2.tick_params(axis='both',top='off', which='both',
#            length=0,width=0,labelsize=0,pad=8)
#ax2.grid(1, ls='--', lw=1, alpha=0.5,
#                color=black, which='minor')
#ax2.grid(1, ls='-', lw=1, alpha=0.5, 
#            color=red, which='major')
#ax2.spines['top'].set_visible(0)
#ax2.set_axis_bgcolor('white')

#ax2.set_xlim(-.95,.95)
#ax2.set_xlim(0,.95)
#ax2.text(tx, ty, s='B', fontsize=12, color='r')

####################################
#subplot(313)
#plot_lines()
#lab3 = r'$G$'


#ylim(-0.1,1.3)
#gca().set_yticks([])
#tick_params(axis='both',top='off', which='major',
#            length=6,width=1,labelsize=fs,pad=8)
#tick_params(axis='both',top='off', which='minor',
#            length=0,width=0,labelsize=fs)


#intercor_full = np.load('inter_mean_full_87943.npy')
#phis = arange( 5000 ) * 2 * pi / 5000.
#th = arcsin( 2.668 * 1.442 / 4 / pi )
#cpsi = cos(phis) * cos(th)**2 + sin(th)**2
#I_noise = interp1d( cpsi, intercor_full )
#noise = intercor_full[ cpsi_vals.min() < cpsi_vals < -cpsi_vals.min() ] 

#noise = intercor_full[ logical_and( cpsi >  cpsi_vals.min(),cpsi < -cpsi_vals.min() )]
#cpsi_ = cpsi[ logical_and( cpsi >  cpsi_vals.min(),cpsi < -cpsi_vals.min() )]
#I_noise( cpsi_vals)

#ylim(-0.1,1.3)
#gca().set_yticks([])
#tick_params(axis='both',top='off', which='major',
#            length=6,width=1,labelsize=fs,pad=8)
#tick_params(axis='both',top='off', which='minor',
#            length=0,width=0,labelsize=fs)

#ax3 = gca()
#ax3.set_xticks(min_ticks_bot, minor=1)
#ax3.set_xticklabels(min_tlab_bot, minor=1)
#ax3.set_xticks(maj_ticks_bot)
#ax3.set_xticklabels(maj_tlab_bot) 
#xlabel(r'$\cos\, \psi $', fontsize=fs)
#xlim(-.95,.95)

#ax3.set_axis_bgcolor('white')
#ax3.spines['top'].set_visible(0)

#grid(1, ls='-', lw=1, alpha=0.5, color=red, which='major')
#grid(1, ls='--', lw=1, alpha=0.5, color=black, which='minor')

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
#ax3.text(tx, ty, s='C', fontsize=12, color='r')

#subplots_adjust(bottom=.18, hspace=0)

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


