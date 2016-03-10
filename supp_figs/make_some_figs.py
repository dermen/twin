from pylab import *
import pandas
import matplotlib.gridspec as gridspec
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

fig = figure(1, figsize=(5,3))
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
ax1.tick_params(axis='both', length=0,width=1,
        labelsize=fs)
ax1.tick_params(axis='x')
ax1.set_ylabel(r'$p_r$', fontsize=fs, rotation=0)
ax1.yaxis.set_label_position("left")
ax1.set_axis_bgcolor('w')



ax2.plot( rp, arange(len(rp)), 's', color='c',
                ms=5, alpha=0.7, label='raw data')
ax2.plot( rp_fit, arange(len(rp)),color='Darkorange',
                lw=2, label='Gaussian fit') #,alpha=0.6)
ax2.yaxis.tick_right()
ax2.set_yticks(qticks)
ax2.set_yticklabels(qlabels)
ax2.set_xticks([700,1100])

ax2.set_xticklabels( map(lambda x: r'$%d$'%x, \
                [700,1100])  )
ax2.yaxis.set_label_position("right")
ax2.set_ylabel(r'$p_r$', fontsize=fs, rotation=0)
ax2.tick_params(axis='both', length=0,width=1,
        labelsize=fs)
ax2.set_xlabel(r'$\langle I_i(p_r, \phi) \rangle_\phi$',\
            fontsize=fs)
subplots_adjust(bottom=.22,left=.12, 
    right=.89, top=.98,wspace=.08)

ax2.legend(prop={'size':8}, numpoints=1, loc=1)
ax2.set_axis_bgcolor('w')


savefig('bragg_peak_position.png', dpi=150, 
            facecolor='w')
