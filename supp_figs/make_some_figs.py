from pylab import *
import pandas
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
import itertools

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


tableau20 = [(31, 119, 180), (174, 199, 232), 
            (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138),
             (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213),
             (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210),
             (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141),
             (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)  

####################################
#
#  FIGURE: BRAGG RING POSITION     #
#
####################################
#style.use('fivethirtyeight')
style.use('bmh')
df = pandas.read_pickle( 'interped_178802.pkl' )
rp = df['radial_profile'].values[1000].tolist()

Gauss = lambda x,amp, var, mu,offset: \
        amp*np.exp(-(x-mu)**2/(2.*var)) + offset
A,W,C,O = df.ix[1000, ['pk_amp', 'pk_width',\
                'pk_pos', 'background_offset']] 
rp_fit = Gauss( arange(len(rp)), A,W**2,C,O)


im = np.load( '178802_1000.npy')

qval = lambda qinvang : (0.053 / 0.00005)*\
            tan(2* arcsin( qinvang*1.442/4/pi )  )
qmin = int( qval( 2.6 ) )

fig = figure(1, figsize=(6,4))
fs = 12
gs = gridspec.GridSpec(nrows=1, ncols=3)
ax1 = plt.subplot(gs[0, :-1])
ax2 = plt.subplot(gs[0,-1])
qticks = arange( 10,71,20) 
qlabels = map(lambda x: r'$%d$'%x,qticks+qmin)

ax1.set_xticks([int(5000 / 4.), int(5000 * 3. / 4)])
ax1.set_xticklabels([r'$\pi/2$', r'$3\pi/2$'])
ax1.set_yticks(qticks)
ax1.set_yticklabels(qlabels)
ax1.imshow(im, aspect='auto', cmap='hot',
            vmin=500, vmax=2000, origin='lower')
ax1.set_xlabel(r'$\phi\,\, (0-2\pi)$', fontsize=fs)
ax1.tick_params(axis='y', length=6,width=1,
        labelsize=fs, color='#777777')
ax1.tick_params(axis='x')
ax1.yaxis.tick_left()
ax1.set_ylabel(r'$r$', fontsize=fs, rotation=0)
ax1.yaxis.set_label_position("left")
ax1.set_axis_bgcolor('w')



ax2.plot( rp, arange(len(rp)), 's', color='c',
                ms=5, alpha=0.7, label='raw data')
ax2.plot( rp_fit, arange(len(rp)),color='Darkorange',
                lw=2, label='Gaussian fit')
ax2.yaxis.tick_right()
ax2.set_yticks(qticks)
ax2.set_yticklabels(qlabels)
ax2.set_xticks([700,1100])

ax2.set_xticklabels( map(lambda x: r'$%d$'%x, \
                [700,1100])  )
ax2.yaxis.set_label_position("right")
ax2.set_ylabel(r'$r$', fontsize=fs, rotation=0)
ax2.tick_params(axis='both', length=0,width=1,
        labelsize=fs)
ax2.set_xlabel(r'$\langle I_i(r, \phi) \rangle_\phi$',\
            fontsize=fs)

leg = ax2.legend(prop={'size':10}, numpoints=1, loc=1)
fr = leg.get_frame()
fr.set_facecolor('w')
ax2.set_axis_bgcolor('w')

ax1.grid(1, lw=1, color='w', alpha=0.5, ls='--')
ax2.grid(1, lw=1, color='#777777', alpha=0.5, ls='--')

ax1.set_ylim(qticks[0]-8, qticks[-1]+8)
ax2.set_ylim(qticks[0]-8, qticks[-1]+8)


ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

subplots_adjust(bottom=.22,left=.12, 
    right=.89, top=.98,wspace=.08)

savefig('bragg_peak_position.pdf',
            facecolor='w')


#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#sdlfklfdsldfLK#LKL#K$L#KL$K#K#
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#dkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t




####################################
#
#  FIGURE: SHOT MEAN HISTOGRAM     #
#
####################################

df = pandas.read_pickle('8.6keV_shots.pkl')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig = figure( 2, figsize=(5.5,4))

ax = gca()

# make logarithmic bins
nbins= 200
power_ten_start=1
power_ten_stop=4
bins = np.logspace(power_ten_start, 
            power_ten_stop , nbins)

ax.hist( df.shot_mean.values, bins=bins, log=1,  
        histtype='stepfilled' , alpha=0.5,
        color='#348ABD')

ax.plot( ones( 2 )*300, [ 0, 1e5], color=orange, 
                lw=2, ls='-')
ax.plot( ones( 2 )*3000, [ 0, 1e5], color=orange, 
                lw=2, ls='-')

ax.set_xscale('log')
ax.set_xlim(80,7.9e3)
ax.set_ylim(0.2, 8e3)

ax.tick_params(axis='y',which='major', length=6,width=1,
        labelsize=fs, color='#777777')
ax.tick_params(axis='both',which='minor', length=0,width=1,
        labelsize=fs, color='#777777')

ax.set_xlabel(r'$\mathrm{mean\, counts}\, \bar{I_i}$',
            fontsize=fs)
ax.set_ylabel(r'$\mathrm{number\, of\, exposures}$',
            fontsize=fs)

ax.set_axis_bgcolor('w')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)


ax.grid(1, axis='x',ls='--', color='#777777', lw=1, 
        alpha=0.5, which='both')
ax.grid(1, axis='y',ls='--', color='#777777', lw=1, 
        alpha=0.5, which='major')
ax.xaxis.tick_bottom()

subplots_adjust(left=.12,bottom=.15,right=.91,top=.95)
savefig('shot_mean_histogram.pdf')


#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#eg9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t
#s;kdjfklsdfjljsdnflsdf
#sldkjfkdsjfksdnfjksdnksdfvlserv
#ewirjtwp98gnpwegnpw5g9ut893t



##############################
#
#  FIGURE: COR VS DIFCOR     #
#
##############################

style.use('ggplot')

x, cors, difcors = np.load('cor_vs_difcor.npy')
xdata = np.linspace( x.min(), 0.99, 2500)
Icors = interp1d( x, cors, bounds_error=0)
Idifcors = interp1d( x, difcors, bounds_error=0)
Pcors = np.polyfit( xdata, Icors(xdata), deg=6 )
Pdifcors = np.polyfit( xdata, Idifcors(xdata), deg=6 )
cors_flat = Icors(xdata) - polyval(Pcors, xdata )
difcors_flat = Idifcors(xdata) - polyval(Pdifcors, xdata )
cors_norm = norm_corr( cors_flat)
difcors_norm = norm_corr( difcors_flat)

fs = 12
fig,ax = subplots(nrows=2, ncols=2, num=3,figsize=(6,4))
ax00 = ax[0,0]
ax01 = ax[0,1]
ax10 = ax[1,0]
ax11 = ax[1,1]


ticks = [ -7/9., -5/9., -1/3., 1/3. , 5/9., 7/9.  ]
ticklabels = [ '-7/9', '-5/9', '-1/3', 
                '1/3' , '5/9', '7/9'  ]
ticklabels = map( lambda lab: r'$%s$'%lab, 
                    ticklabels )



ax00.plot(  xdata, Icors(xdata), lw=8 , color='b')
ax00.plot(  xdata, polyval(Pcors, xdata), 
                lw=3 , color='y', ls='--')
ax00.set_xlim(-0.84,0.84)
ax00.set_ylim(-0.008,0.006)
ax00.set_yticks([])
ax00.set_xticklabels([])
ax00.set_xticks(ticks)
ax00.grid(1,lw=1,ls='--',color='#777777',alpha=.5)
ax00.tick_params(axis='x', length=0)
ax00.legend()
ax00.set_axis_bgcolor('w')
ax00.set_ylabel(r'$C(\cos \,\psi)$', fontsize=12)
#lab00 = ( r'$C(\cos \,\psi)$', r'$\mathrm{fit}$'  )
#leg=ax00.legend( lab00, loc=9 )
#fr = leg.get_frame()
#fr.set_facecolor('w')
#fr.set_alpha(0.5)

ax01.plot(  xdata, Idifcors(xdata), lw=5 , color='m')
ax01.plot(  xdata, polyval(Pdifcors,xdata), 
            lw=4 , color='y', ls='--')
ax01.set_xlim(-0.84,0.84)
ax01.set_ylim(-0.00018,0.00005)
ax01.set_xticklabels([])
ax01.tick_params(axis='x', length=0)
ax01.set_yticks([])
ax01.set_xticks(ticks)
ax01.grid(1,lw=1,ls='--',color='#777777', alpha=.5)
ax01.set_axis_bgcolor('w')
ax01.set_ylabel(r'$D(\cos \, \psi)$',fontsize=fs)
            #rotation=180)
#ax01.yaxis.set_label_position('right')


ax10.plot( xdata, cors_norm, color='b', lw=2 )
ax10.set_xlim(-0.84,0.84)
ax10.set_ylim(-.1,1.1)
ax10.xaxis.tick_bottom()
ax10.set_yticks([])
ax10.tick_params(axis='x', length=6, labelsize=fs)
ax10.grid(1,lw=1,ls='--',color='#777777',alpha=.5)
ax10.set_xticks(ticks)
ax10.set_xticklabels(ticklabels, rotation=45)
ax10.set_xlabel(r'$\cos\, \psi$', fontsize=fs,labelpad=-5)
ax10.set_axis_bgcolor('w')
ax10.set_ylabel(r'$\mathrm{residual}$',fontsize=fs)

ax11.plot( xdata, difcors_norm, color='m' , lw=2)
ax11.set_xlim(-0.84,0.84)
ax11.set_ylim(-.1,1.1)
ax11.xaxis.tick_bottom()
ax11.set_yticks([])
ax11.tick_params(axis='x', length=6, labelsize=fs)
ax11.set_xticks(ticks)
ax11.set_xticklabels(ticklabels, rotation=45)
ax11.set_xlabel(r'$\cos\, \psi$', fontsize=fs,labelpad=-5)
ax11.grid(1,lw=1,ls='--',color='#777777',alpha=.5)
ax11.set_axis_bgcolor('w')
ax11.set_ylabel(r'$\mathrm{residual}$',fontsize=fs)

subplots_adjust(left=.09,bottom=.16,right=.97,top=.97,
        wspace=.29, hspace=.2)

#savefig('compare_correlations.pdf')
##################################
#  FIGURE: FRIEDEL COMPARE       #
# (uses stuff from previous fig) # 
##################################
ydata = Idifcors(xdata)
x_friedel, difcors_friedel = friedel_corr(xdata, 
            ydata)

ydata_norm = norm_corr(ydata)

ticks = [pi/4, pi/2, 3*pi/4]
ticklabels = [ r'$\pi/4$', r'$\pi/2$',
                r'$3\pi/4$']

fig, ax = subplots(nrows=2,ncols=1,num=4,figsize=(6,4) )
ax0 = ax[0]
ax1 = ax[1]

#x_m = xdata.min()

xlm = (0.62302099794755916, 2.4980604661828156)

xpsi = arccos( xdata)
order = argsort( xpsi )

ax0.plot( xpsi[order], ydata_norm[order],
        color=purple, lw=2 ) 
#ax0.set_xlim(arccos(x_m),arccos(-x_m))
ax0.set_xlim(xlm)
ax0.set_yticks([])
ax0.set_ylim(.12,1.05)
ax0.set_xticks(ticks)
#ax0.set_xticklabels(ticklabels)
ax0.set_xticklabels([])
ax0.tick_params(length=0)
ax0.set_ylabel(r'$D(\psi)$',fontsize=fs)
ax0.grid(1,lw=1,ls='--',color='#777777',alpha=.5)
ax0.set_axis_bgcolor('w')


xpsi_friedel = arccos( x_friedel)
order = argsort( xpsi_friedel )

ax1.plot( xpsi_friedel[order],
        norm_corr(difcors_friedel)[order],
            color=green, lw=2)

ax1.set_xticks(ticks)
ax1.set_xticklabels(ticklabels)

ax1.set_yticks([])
ax1.tick_params(length=6, labelsize=fs, axis='x')

#ax1.set_xlim(arccos(x_m),arccos(-x_m))
ax1.set_xlim(xlm)
ax1.set_ylim(-0.3,1.3)
ax1.set_ylabel(r'$D_F( \psi)$',fontsize=fs)
ax1.grid(1,lw=1,ls='--',color='#777777',alpha=.5)
ax1.xaxis.tick_bottom()
ax1.set_xlabel(r'$\mathrm{angle}\, \psi\,(\mathrm{radians})$',
            fontsize=fs,labelpad=2)
ax1.set_axis_bgcolor('w')

subplots_adjust(left=.07,bottom=.14,right=.96,top=.95,
        wspace=.2, hspace=.09)

savefig('compare_friedel.pdf')




####################
# GAUSSIAN FITTING #
####################
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from scipy import optimize
from scipy.special import wofz
from itertools import cycle
import sys
from lmfit import minimize, Parameters, Parameter, report_fit

sys.path = ['/Users/mender/loki']+ sys.path
import postproc_helper
sm = postproc_helper.smooth


xx, C = np.load( '/Users/mender/Documents/twin/sym-corsFeb29-0pix.npy')
xx, cn = np.load( '/Users/mender/Documents/twin/sym-mean-corsFeb29-0pix.npy')
cnn = positive_corr( cn).data
fit_data = savgol_filter( cnn, 31, 2 )
fit_data = sm( fit_data, 30, 50 )

local_max = argrelextrema(fit_data, np.greater, order=1)[0]
N_peaks = len( local_max)

local_mins = argrelextrema(fit_data, np.less, order=1)[0]
local_mins = np.append( [0], local_mins )
local_mins = np.append( local_mins, [cn.shape[0]-1 ] )
while local_max[0] < local_mins[0]:
    local_max = local_max[1:]
while local_max[-1] > local_mins[-1]:
    local_max = local_max[:-1] 

fig = figure(1, figsize=(6,6))
ax = gca()

labels = (r'$\mathrm{Savitzski-Golay}$',
        r'$\cos \,\psi_\gamma $',
        r'$\mathrm{local\, min}$',
            r'$\mathrm{raw\, data}$' )

ax.plot( xx, cnn, 'o', ms=5,color=tableau20[18], 
            label=labels[3])
ax.plot( xx, fit_data, color=tableau20[5], lw=3, 
            label=labels[0] )
ax.plot( xx[local_max], fit_data[local_max],
        's',color=tableau20[3], ms=9, 
        label=labels[1], alpha=0.8)
#ax.plot( xx[local_mins], fit_data[local_mins], 
#        '^',color=tableau20[11], ms=10,
#            label=labels[2], alpha=0.8)


leg = ax.legend(numpoints=1, loc=2)
fr = leg.get_frame()
fr.set_alpha(0.3)
fr.set_facecolor('w')
ax.set_axis_bgcolor('w')
ax.set_xlim(-0.8,-0.2)

ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
ax.tick_params(axis='both', length=6, labelsize=fs )

ax.set_xticklabels(ticklabels)

ax.set_xlabel(r'$\cos \, \psi$',fontsize=fs, 
            labelpad=10)

ax.set_yticks([])
ax.set_ylabel(r'$D_F(\cos \, \psi)$',
            fontsize=12, labelpad=10)
subplot_adjust(top=.97)
#savefig('peak_detectction.pdf')


###########################
# GAUSSIAN PARTIAL FITTING #
###########################
def sum_gauss( params, xdata, ydata , Ngauss, ctrs):
    model = zeros_like( xdata)
    for i in xrange( Ngauss):
        amp = params['amp%d'%i ].value
        wid = params['wid%d'%i ].value
        off = params['off%d'%i ].value
        model = model+  amp * np.exp( \
                    -((xdata - ctrs[i])/wid)**2) + off
    return model - ydata

colors_fit = cycle(tableau20[0::2])
colors = cycle(tableau20[1::2])
markers = cycle(  ['^', 's', 'o', 'd', '<']  )

local_max_split = np.array_split( local_max, 10)
args_data = []
kwargs_data = []
args_fit = []
kwargs_fit = []

mecs =cycle([blue, 'Darkorange', green, red, purple  ])

for local_max_ in local_max_split:
    lowest_max = local_max_[0]
    highest_max = local_max_[-1]
    lower_bound = local_mins[ local_mins < lowest_max][-1]
    upper_bound = local_mins[ local_mins > highest_max][0]

    ctrs = xx[ local_max_]
    xdata = xx[ np.arange(lower_bound, upper_bound +1)]
    ydata = cnn[ np.arange(lower_bound, upper_bound +1)]
   

    args = [xdata, ydata, markers.next()]
    kwargs = {'ms':8, 'mec':mecs.next() ,'mew':1,
            'color':colors.next(),
                'alpha':0.5}
       
    args_data.append( args)
    kwargs_data.append(kwargs)
    
    W = 0.01252791440815894
    N_gauss = len(local_max_)
    params = Parameters()
    for i_ in xrange(N_gauss):
        params.add('off%d'%i_, value= 0.1, min=0, max=2)
        params.add('wid%d'%i_, value= W, 
                    min=W-W*0.7, max=W+W*0.7)
        params.add('amp%d'%i_, value= 1 , min=0, max=2)
    
    result_gauss = minimize(sum_gauss, params, 
                args=(xdata, ydata, N_gauss, ctrs ))

    rp = result_gauss.params
    offs = zeros_like(ctrs)
    amps = zeros_like(ctrs)
    wids = zeros_like(ctrs)
    for i_ in xrange( N_gauss):
        offs[i_] = rp['off%d'%i_].value
        amps[i_] = rp['amp%d'%i_].value
        wids[i_] = rp['wid%d'%i_].value

    fit_partial = sum_gauss( result_gauss.params, 
                xdata, zeros_like( xdata),
                                N_gauss, ctrs )
    c = colors_fit.next()

    args = [xdata, fit_partial-sum(offs)]
    kwargs = {'color':c,'lw':3}
    args_fit.append( args)
    kwargs_fit.append(kwargs)


fig = figure( 5, figsize=(6,6) )
fs=12
ax = gca()
for arg_data, arg_fit,\
    kwarg_data, kwarg_fit in zip( args_data, \
            args_fit, kwargs_data, kwargs_fit ):
    ax.plot( *arg_data, **kwarg_data)
    ax.plot( *arg_fit, **kwarg_fit)

ax.set_xlim(-.82,.01)
ax.set_ylim(-0.05,2.05)
ax.set_xlabel(r'$\cos \, \psi$', fontsize=fs,
        labelpad=7)

ticks = map( lambda x: r'$%.1f$'%x, ax.get_xticks())
ax.set_xticklabels(ticks)

ax.grid(1, lw=1, ls='--', color='#777777', alpha=.5)
ax.set_axis_bgcolor('w')
ax.tick_params(length=6, labelsize=fs)
ax.set_yticks([])
ax.xaxis.tick_bottom()
subplots_adjust(top=.98, right=.96, left=.04)
savefig('gaussian_partial.pdf')




##################
#
#  WIDTH OF PEAK #
#
###################

style.use('ggplot')

cpsi_model_full, C_model_full = np.load( '../model.npy')
C_model_full = norm_corr( C_model_full)
cpsi_model, C_model = friedel_corr( cpsi_model_full,
                        C_model_full)

xsim = arccos(cpsi_model[390:450])
ysim = C_model[390:450]

#xsim = sort( xsim )
#ysim = ysim[ argsort(xsim)]

gauss = lambda x,ctr,width,amp,off: off + \
                amp*exp(- ((x-ctr)/width)**2 )

fit_sim = optimize.curve_fit(gauss, xdata=xsim,
                ydata=ysim, p0=(arccos(5/9.),
                        .06, 0.6, 0 ) )

width_sim = fit_sim[0][1]
beta_sim = 2.3458 * width_sim / sqrt(2)
wavelen = 1.442
q = 2.668
th = arcsin( q* wavelen / 4 / pi )
s_sim = .89 * sqrt(2) * wavelen / beta_sim / cos(th)
#s_sim = 41.6 inverse angstroms

cpsi_data, C_data = \
        np.load( '../sym-mean-corsFeb29-0pix.npy')
C_data = norm_corr(C_data)


xdat = arccos(cpsi_data[2050:2170])
ydat = C_data[2050:2170]

#xdat = sort(xdat)
#ydat = ydat[argsort(xdat)]

fit_dat = optimize.curve_fit(gauss, xdata=xdat,
                ydata=ydat, p0=(arccos(5/9.),0.03, 
                        0.7, 0.15 ) )

width_dat = fit_dat[0][1]
beta_dat = 2.3458 * width_dat / sqrt(2)
s_dat = .89 * sqrt(2) * wavelen / beta_dat / cos(th)


ctr_dat = fit_dat[0][0]
amp_dat = fit_dat[0][2]

amp_sim = fit_sim[0][2]
ctr_sim = fit_sim[0][0]


xoffset = -ctr_sim + ctr_dat

fig,ax = subplots(nrows=2, ncols=1, num=7, figsize=(5,5) )
ctr = arccos(5/9.)
ticks = ctr + linspace( -0.03, 0.03, 7 )

ticklabels = map( lambda x: (r'%.3f'%x).strip('0'),
                ticks)
ticklabels = map( lambda x: (r'$%s$'%x),
                ticklabels)
ax1 = ax[0]
ax2 = ax[1]

# PLOT THE SIMULATION
ctr_simoff = ctr_sim + xoffset
ax1.plot(xsim+xoffset, ysim, '^',
        label=r'$\mathrm{simulation}$',
        ms=9, mec='b',alpha=.5, color=blue)
ax1.plot(  xsim+xoffset, gauss(xsim, *fit_sim[0] ), 
                lw=3, ls='--', color=blue)

ax1.annotate(
    '', xy=(ctr_simoff-beta_sim/2.,
            amp_sim/2.), xycoords='data',
    xytext=(ctr_simoff+beta_sim/2., 
            amp_sim/2.), textcoords='data',
    arrowprops={'arrowstyle': '<->', 'color':'k', 'lw':1})
ax1.text(x=.986, y=amp_sim/2 + 0.05, 
        s=r'$\beta^*_{\mathrm{sim}}$', fontsize=fs)

#ax1.set_xlim(.94344, 1.02043)
ax1.set_ylim(0,.65)

ax1.grid(1, lw=1, ls='--', color='#777777', alpha=.5)
ax1.set_axis_bgcolor('w')
ax1.tick_params(length=0, labelsize=fs)
ax1.set_yticks([])
ax1.xaxis.tick_bottom()
ax1.set_xticks(ticks)
ax1.set_xticklabels([])
ax1.set_xlim(0.9417653565786227, 1.0217653565786227 )

# PLOT THE DATA
ax2.plot(xdat, ydat - fit_dat[0][3], 'o',
        label=r'$\mathrm{simulation}$' ,
        ms=8,mec='r', alpha=.15, color=red)
ax2.plot(  xdat, gauss(xdat, *fit_dat[0] ) - fit_dat[0][3], 
                lw=3, ls='-', color=red)
ax2.annotate(
    '', xy=(ctr_dat-beta_dat/2.,
            amp_dat/2.), xycoords='data',
    xytext=(ctr_dat+beta_dat/2., 
            amp_dat/2.), textcoords='data',
    arrowprops={'arrowstyle': '<->', 'color':'k', 'lw':1})

ax2.text(x=.983, y=amp_dat/2 + 0.05, 
        s=r'$\beta^*_{\mathrm{data}}$', fontsize=fs)
ax2.set_xlim(0.9417653565786227, 1.0217653565786227 )
#ax2.set_xlim(.94344, 1.02043)
ax2.set_ylim(0,.85)

ax2.grid(1, lw=1, ls='--', color='#777777', alpha=.5)
ax2.set_axis_bgcolor('w')
ax2.tick_params(length=0, labelsize=fs)
ax2.set_yticks([])
ax2.xaxis.tick_bottom()

ax2.set_xlabel(r'$\psi\,(\mathrm{radians})$', fontsize=fs,
        labelpad=9)

ax2.set_xticks(ticks)
ax2.set_xticklabels(ticklabels)
subplots_adjust(top=.96, right=.96, left=.1, 
            bottom=.12,hspace=0)
fig.text(0.04, 0.5, r'$D_F(\psi)$', va='center',
        rotation='vertical', fontsize=fs)

savefig('compare_widths.pdf')
#
