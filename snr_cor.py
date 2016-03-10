from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from scipy import optimize
from scipy.special import wofz
from itertools import cycle

from lmfit import minimize, Parameters, Parameter, report_fit
import postproc_helper
sm = postproc_helper.smooth

def norm_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma.min()
    cmax = cor_vec_ma.max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n #/ crange

xx, C = np.load( '/Users/mender/Documents/twin/sym-corsFeb29-0pix.npy')
xx, cn = np.load( '/Users/mender/Documents/twin/sym-mean-corsFeb29-0pix.npy')
cnn = norm_corr( cn).data
cnn_smooth = savgol_filter( cnn, 31, 2 )
cnn_smooth = sm( cnn_smooth, 30, 50 )

local_max = argrelextrema(cnn_smooth, np.greater, order=1)[0]

local_mins = argrelextrema(cnn_smooth, np.less, order=1)[0]
local_mins = np.append( [0], local_mins )
local_mins = np.append( local_mins, [cn.shape[0]-1 ] )
while local_max[0] < local_mins[0]:
    local_max = local_max[1:]
while local_max[-1] > local_mins[-1]:
    local_max = local_max[:-1] 

plot( xx, cnn, 'kx')
plot( xx, cnn_smooth, 'b', lw=3 )
plot( xx[local_max], cnn_smooth[local_max], 's',color='Limegreen', ms=8)
plot( xx[local_mins], cnn_smooth[local_mins], 'd',color='r', ms=8)

#shot_series , cnn_series_sets = np.load( 'cnns.npy')


shot_series = array([   25,  3150,  6275,  9400, 12525, 15650, 18775])

#shot_series = array([ 200*i for i in xrange(1,100) ])
#cnn_series_ = [ [norm_corr( C[np.random.randint(0,C.shape[0], 
#                s ) ].mean(0) ).data for s in shot_series]
#                for jj in range(400) ]
#cnn_series_sets_ = array( cnn_series)
#np.save( 'cnns_200', [shot_series, cnn_series_sets] )


N = C.shape[0]
#N_sets_all = map( int, [ float(N/s) for s in shot_series  ] )

N_sets_all = [ 100 for s in shot_series]
N_sets_all[0] = 1
#N_sets = cnn_series_sets.shape[0]
N_nshots = len( shot_series)
N_peaks = len( local_max)


#aksjdtkjahfl
############################
# GAUSSIAN SUM FIT CODE    #
############################
def sum_gauss( params, xdata, ydata , Ngauss):#, ctrs):
    model = zeros_like( xdata)
    for i in xrange( Ngauss):
        amp = params['amp%d'%i ].value
        wid = params['wid%d'%i ].value
        off = params['off%d'%i ].value
        ctr = params['ctr%d'%i ].value
        #model = model+  amp * np.exp( -((xdata - ctrs[i])/wid)**2) + off
        model = model+  amp * np.exp( -((xdata - ctr)/wid)**2) + off
    return model - ydata


local_max_split = np.array_split( local_max, 10)
ctrs = xx[ local_max]
fit_data4 = {}
fit_data4['amps'] = {}
fit_data4['offs'] = {}
fit_data4['wids'] = {}
fit_data4['ctrs'] = {}
for i_nshot in xrange( 0, len(shot_series)):
    #np.random.shuffle(C)
    N_sets = N_sets_all[i_nshot]
    fit_data4['amps'][i_nshot] = zeros( (N_sets, N_peaks))
    fit_data4['offs'][i_nshot] = zeros( (N_sets, N_peaks))
    fit_data4['wids'][i_nshot] = zeros( (N_sets, N_peaks))
    fit_data4['ctrs'][i_nshot] = zeros( (N_sets, N_peaks))

    for i_max, local_max_ in enumerate( local_max_split):
        lowest_max = local_max_[0]
        highest_max = local_max_[-1]
        lower_bound = local_mins[ local_mins < lowest_max][-1]
        upper_bound = local_mins[ local_mins > highest_max][0]

        xdata = xx[ np.arange(lower_bound, upper_bound +1)]
        local_max_indices = searchsorted( local_max, local_max_ )
        ctrs = xx[ local_max_]
        
        W = 0.01252791440815894
        N_gauss = len(local_max_)
        params = Parameters()
        
        C_partial = C[:,lower_bound: upper_bound+1]
        #C_partial_chunks = np.array_split( C_partial, N_sets)
        inds = [ np.random.randint(0,N , shot_series[i_nshot])
                    for _ in xrange(N_sets)]
        C_partial_chunks = array([ C_partial[ ind ].mean(0)
                            for ind in inds ])

        for i_, ctr in enumerate( ctrs): # in xrange(N_gauss):
            params.add('off%d'%i_, value=0.1, min=0, max=40)
            params.add('wid%d'%i_, value=W, min=W-W*0.7, max=W+W*0.7)
            params.add('amp%d'%i_, value=1 , min=0, max=40)
            #ctr = ctrs [local_max_[i_] ]
            params.add('ctr%d'%i_, value=ctr , min=ctr-0.01*ctr,\
                                            max=ctr+0.01*ctr)
        
        #for i_nshot in xrange( N_nshots):
        #N_sets = N_sets_all[i_nshot]

        for i_set, cnn_set in enumerate( C_partial_chunks):
            ydata = norm_corr( cnn_set )
            #ydata = norm_corr( cnn_set.mean(0) )
            #ydata = cnn_mean[ np.arange(lower_bound, upper_bound +1)]
            result_gauss = minimize(sum_gauss, params, 
                        args=(xdata, ydata, N_gauss ))
            
            rp = result_gauss.params
            for i_, i_peak in enumerate(local_max_indices):
                fit_data4['offs'][i_nshot][i_set, i_peak] =\
                            rp['off%d'%i_].value
                fit_data4['amps'][i_nshot][i_set, i_peak] =\
                            rp['amp%d'%i_].value
                fit_data4['wids'][i_nshot][i_set, i_peak] =\
                            rp['wid%d'%i_].value
                fit_data4['ctrs'][i_nshot][i_set, i_peak] =\
                            rp['ctr%d'%i_].value



#===========================
import h5py
fo = h5py.File('fit_data4.h5py', 'w')

for param, param_dict in fit_data4.iteritems():
    for param_number, param_val in param_dict.iteritems():
        h5_key = param +'/'+ param+ str(param_number)
        fo.create_dataset(h5_key, data = param_val )

fo.close()

#==========================================================
#
# FIND THE CONVERGED SIGNAL AMPLITUDES
#
#=========================================================

f2 = h5py.File('fit_data2.h5py', 'r')
f3 = h5py.File('fit_data3.h5py', 'r')
f4 = h5py.File('fit_data4.h5py', 'r')

amps2 = f2['amps']
amps3 = f3['amps']
amps4 = f4['amps']

wids2 = f2['wids']
wids3 = f3['wids']
wids4 = f4['wids']

ctrs2 = f2['ctrs']
ctrs3 = f3['ctrs']
ctrs4 = f4['ctrs']

#amps_all = {}
A_means = []
W_means = []
C_means = []
Num_samples = zeros_like(shot_series, dtype=float)
for i_shot,nshot in enumerate(shot_series):
    keyA = 'amps%d'%i_shot
    keyW = 'wids%d'%i_shot
    keyC = 'ctrs%d'%i_shot
    amps = ( amps2[keyA], amps3[keyA],amps4[keyA] )
    wids = ( wids2[keyW], wids3[keyW], wids4[keyW])
    ctrs = ( ctrs2[keyC], ctrs3[keyC], ctrs4[keyC])
    
    A = concatenate(amps,0)
    W = concatenate(wids,0)
    C = concatenate(ctrs,0)
    #amps_all[i_shot] = A
    
    N_ = np.arange( 1,A.shape[0]+1 , dtype=float)
    Num_samples[i_shot] = A.shape[0]

    A_ = np.cumsum( A, 0)
    A_means.append(  A_.T / N_ ) 
    
    W_ = np.cumsum( W, 0)
    W_means.append(  W_.T / N_ ) 

    C_ = np.cumsum( C, 0)
    C_means.append(  C_.T / N_ ) 


A_mean_after_nshots = [cm.mean(1)
                    for cm in A_means]
A_stdev_after_nshots = [cm.std(1)
                    for cm in A_means]
W_mean_after_nshots = [cm.mean(1)
                    for cm in W_means]
W_stdev_after_nshots = [cm.std(1)
                    for cm in W_means]
C_mean_after_nshots = [cm.mean(1)
                    for cm in C_means]
C_stdev_after_nshots = [cm.std(1)
                    for cm in C_means]

amp_79 = [ a_[0] for a_ in A_mean_after_nshots]
ctr_79 = [ c_[0] for c_ in C_mean_after_nshots]
wid_79 = [ w_[0] for w_ in W_mean_after_nshots]

amp_13 = [ a_[13] for a_ in A_mean_after_nshots]
ctr_13 = [ c_[13] for c_ in C_mean_after_nshots]
wid_13 = [ w_[13] for w_ in W_mean_after_nshots]

amp_4 = [ a_[12] for a_ in A_mean_after_nshots]
ctr_4 = [ c_[12] for c_ in C_mean_after_nshots]
wid_4 = [ w_[12] for w_ in W_mean_after_nshots]

amp_59 = [ a_[8] for a_ in A_mean_after_nshots]
ctr_59 = [ c_[8] for c_ in C_mean_after_nshots]
wid_59 = [ w_[8] for w_ in W_mean_after_nshots]


CRV = optimize.curve_fit
#fit amplitudes (first point is really noisy)
amp_79_fit = CRV(fitfunc,
            ydata=amp_79[3:],
            xdata=2*shot_series[3:],
            p0=(0.5,0.1))

amp_59_fit = CRV(fitfunc,
            ydata=amp_59[3:],
            xdata=2*shot_series[3:],
            p0=(0.5,0.1))

amp_13_fit = CRV(fitfunc,
            ydata=amp_13[3:],
            xdata=2*shot_series[3:],
            p0=(0.5,0.1))

amp_4_fit = CRV(fitfunc,
            ydata=amp_4[3:],
            xdata=2*shot_series[3:],
            p0=(0.5,0.1))

amp_79_eval = fitfunc(2*shot_series, *amp_79_fit[0])
amp_59_eval = fitfunc(2*shot_series, *amp_59_fit[0])
amp_13_eval = fitfunc(2*shot_series, *amp_13_fit[0])
amp_4_eval = fitfunc(2*shot_series, *amp_4_fit[0])

# inter_shot_series
shot_series_ic = array([ 25*i for i in xrange(1,1000) ])
ic = np.load('inter_means.npy')
ic_inds = np.searchsorted( shot_series_ic, shot_series)
noise = ic[ ic_inds]

noise_fit = CRV(fitfunc,
            ydata=noise.std(1), 
            xdata=2*shot_series, 
            p0=(0.5,-0.4) )

noise_eval = fitfunc( 2*shot_series, *noise_fit[0])

#noise_fitfull = CRV(fitfunc,
#            ydata=ic.std(1), 
#            xdata=2*shot_series_ic, 
#            p0=(0.5,-0.4) )
#
#noise_evalfull = fitfunc( 2*shot_series_ic, *noise_fitfull[0])


snr_79_fit = CRV(fitfunc,
                ydata=amp_79_eval/noise_eval,
                xdata=2*shot_series,
                p0=(0.5,0.5))

snr_59_fit = CRV(fitfunc,
                ydata=amp_59_eval/noise_eval,
                xdata=2*shot_series,
                p0=(0.5,0.5))


snr_13_fit = CRV(fitfunc,
                ydata=amp_13_eval/noise_eval,
                xdata=2*shot_series,
                p0=(0.5,0.5))

snr_4_fit = CRV(fitfunc,
                ydata=amp_4_eval/noise_eval,
                xdata=2*shot_series,
                p0=(0.5,0.5))

snr_79_eval = fitfunc(2*shot_series, *snr_79_fit[0])
snr_59_eval = fitfunc(2*shot_series, *snr_59_fit[0])
snr_13_eval = fitfunc(2*shot_series, *snr_13_fit[0])
snr_4_eval = fitfunc(2*shot_series, *snr_4_fit[0])

#======================
# GET ERROR BAR DATA
std_79 = [a_[0] for a_ in A_stdev_after_nshots ]
std_59 = [a_[8] for a_ in A_stdev_after_nshots ]
std_13 = [a_[13] for a_ in A_stdev_after_nshots ]
std_4 = [a_[12] for a_ in A_stdev_after_nshots ]

#err_79 = std_79 / sqrt(Num_samples) / noise.std(1)
#err_59 = std_59 / sqrt(Num_samples) / noise.std(1)
#err_13 = std_13 / sqrt(Num_samples) / noise.std(1)
#err_4 = std_4 / sqrt(Num_samples) / noise.std(1)

err_79 = std_79 /  noise.std(1)
err_59 = std_59 /  noise.std(1)
err_13 = std_13 / noise.std(1)
err_4 = std_4 / noise.std(1)


#======================
style.use('ggplot')
color_cycle = cycle( rcParams['axes.color_cycle'] )
red = color_cycle.next()
blue = color_cycle.next()
purple = color_cycle.next()
black = color_cycle.next()
yellow = color_cycle.next()
green = color_cycle.next()
pink = color_cycle.next()

red = red
orange = 'Darkorange'
blue = 'c'
green = '#89E894'

ms = 10
mew = 1
mec='k'#black
ec='k' #black
elw=2
lw=2

fig = figure(1, figsize=(6,6), dpi=80)

lab13 = r'$\cos \, \psi=1/3$; $k=%.2f$'%np.round(snr_13_fit[0][1],2)
lab59 = r'$\cos \, \psi=5/9$; $k=%.2f$'%np.round(snr_59_fit[0][1],2)
lab79 = r'$\cos \, \psi=7/9$; $k=%.2f$'%np.round(snr_79_fit[0][1],2)
lab4 = r'$\cos \, \psi=0.4$; $k=%.2f$'%np.round(snr_4_fit[0][1],2)

errorbar(x=2*shot_series, 
        y=snr_13_eval,
        yerr=err_13, marker='s', 
        color=red , mew=mew, ms=ms, mec=mec,
        lw=lw,elinewidth=elw,ecolor=ec, label=lab13)

errorbar(x=2*shot_series, 
        y=snr_59_eval,
        yerr=err_59, marker='^', 
        color=blue , mew=mew, ms=ms,mec=mec,
        lw=lw,elinewidth=elw,ecolor=ec,label=lab59)
errorbar(x=2*shot_series, 
        y=snr_79_eval,
        yerr=err_79, marker='o', 
        color=orange , mew=mew, ms=ms,mec=mec,
        lw=lw,elinewidth=elw,ecolor=ec, label=lab79)
errorbar(x=2*shot_series, 
        y=snr_4_eval,
        yerr=err_4, marker='D', 
        color=green , mew=mew, ms=ms,mec=mec,
        lw=lw,elinewidth=elw,ecolor=ec, label=lab4)


ax = gca()
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
leg = legend(handles,labels,numpoints=1, 
            loc=2, borderpad=0)
fr = leg.get_frame()
fr.set_facecolor('w')
fr.set_alpha(0.5)


ax.set_yscale('log')
ax.set_xscale('log')
ax.set_axis_bgcolor('white')
ax.yaxis.tick_left()

fs = 12
tick_params(axis='both', which='major',length=6,width=1,labelsize=fs)
tick_params(axis='both', which='minor',length=0,width=0,labelsize=fs)
grid(b=1, which='both', color=black, linestyle='--', alpha=0.3, lw=1)

xlabel(r'$N$',fontsize=fs)
ylabel(r'$Z(\cos\, \psi)\,\, \propto\,\, N^{k}$', fontsize=fs)
ylim(0.5, 20)
xlim(6000,40000)
ax.xaxis.tick_bottom()

ax.set_xticks([7000,10000,20000,30000 ])
#ax.set_xticklabels([r'$\mathdefault{7} $\cdot$ $\mathdefault{10^{3}}$',r'$\mathdefault{1 \cdot 10^{4}}$',
#                    r'$\mathdefault{2 \cdot 10^{4}}$',r'$\mathdefault{3 \cdot 10^{4}}$'])
ax.set_xticklabels([r'$7 \cdot 10^{3}$',r'$1 \cdot 10^{4}$',
                    r'$2 \cdot 10^{4}$',r'$3 \cdot 10^{4}$'])

ax.set_yticks([1,5,10])
ax.set_yticklabels([r'$1$', r'$5$', 
                    r'$10$'])

subplots_adjust(left=0.15,bottom=.1,right=.91,top=1)



## ===========
# Z-score

#
shots[ argmin( np.abs( fitfunc( shots, *snr_4_fit[0]) - 2.5 ) ) ]
shots[ argmin( np.abs( fitfunc( shots, *snr_13_fit[0]) - 2.5 ) ) ]
shots[ argmin( np.abs( fitfunc( shots, *snr_59_fit[0]) - 2.5 ) ) ]
shots[ argmin( np.abs( fitfunc( shots, *snr_79_fit[0]) - 2.5 ) ) ]
"""

#================================
# PLOT CONVERGENCES OF AMPLITUDES
#================================
pk_lab = ['7/9', '5/9', '0.4', '1/3']
for i_, ipk in enumerate( [0,8,12,13 ] ):

    for i_shot in xrange( 1,7):
        plot( A_means[i_shot][ipk],lw=2 )

    labs =[ r'$N=%d$'%(int(s))  for s in  2*shot_series]
    leg = legend(labs, loc=4)
    fr = leg.get_frame()
    fr.set_alpha(0.5)

    suptitle('Convergence of CXS signal $(S)$ for\npeak at $\cos \,\psi = %s$'%pk_lab[i_],fontsize=12)
    ylabel(r'$S(N; N_{G})$',fontsize=fs)
    xlabel(r'$N_G$',fontsize=fs)
    savefig('/Users/mender/Desktop/SNR_%d.png'%i_,dpi=150 )
    clf()

#plot( shot_series , snr_4_eval, 's' )
#f = h5py.File( 'stats_test2.hdf5', 'w')
#f.create_dataset('off', data=offs_out)
#f.create_dataset('amp', data=amps_out)
#f.create_dataset('wid', data=wids_out)
#f.create_dataset('ctrs', data=xx[local_max])
#f.close()


#f = h5py.File('stats_test2.hdf5', 'r')
#amps = f['amp'].value
#offs = f['off'].value
#wids = f['wid'].value
#ctrs = f['ctrs'].value
#f.close()

def fitfunc(x,p1,p2):
    y = p1*(np.power(x,p2))
    return y

nz_coefs = []
snr_coefs = []

SNR_mean = zeros( N_nshots  )
#for i_set, cnn_series in enumerate(cnn_series_sets):

#nz = array([ std(c-savgol_filter(c, 31,5) )
#        for c in cnn_series])
nz_fit = optimize.curve_fit( fitfunc,
            ydata=nz, xdata=shot_series, 
            p0=(0.5,-0.4) )

#nz_coef = nz_fit[0][1]    
#nz_coefs.append( nz_coef)

SNR = amps[0] / nz[:,None]

#SNR_mean += SNR[:,0]

snr_fit = optimize.curve_fit(fitfunc, ydata=SNR[:,13], 
                xdata=2*shot_series, p0=(8,0.4))

#snr_coef = snr_fit[0][1]    
#snr_coefs.append( snr_coef)
    
snr_fit = optimize.curve_fit(fitfunc, ydata=SNR_mean, 
                xdata=2*shot_series, p0=(8,0.4))

optimize.curve_fit(fitfunc, ydata=nz, 
                xdata=2*shot_series, p0=(8,0.4))

plot( 2*shot_series, SNR[:,13], 's' )
plot( 2*shot_series, fitfunc( 2*shot_series, *snr_fit[0]  ), lw=2 )

def G(ipeak, n ,iset=0):
    return gauss_wid( xx, amps[iset,n,ipeak], wids[iset,n,ipeak], 
                    xx[local_max[ipeak]], offs[iset,n,ipeak])
def A( ipeak, n, iset=0):
    return amps[iset, n, ipeak]
def O( ipeak, n, iset=0):
    return offs[iset, n, ipeak]

plot( xx, G(0,600)-O(0,600) )
plot( xx, G(13,600)-O(13,600) )
plot( xx[ local_max[13] ], A(13,600), 's' )
plot( xx[ local_max[0] ], A(0,600) , 's')

#cnn1 = cnn_series_sets[0,0]
#cnn25 = cnn_series_sets[0,24]
#cnn50 = cnn_series_sets[0,-1]

#cnns = [cnn1, cnn25, cnn50]

#nn = ['1k', '25k', '50k']
#for i,c in enumerate(cnns):
#    subplot( 3,1,i+1)
#    c_smooth = savgol_filter(c, 31, 2)
#    plot( xx, c, 'sk', ms=5, label='raw data N=%s'%nn[i])
#    plot( xx, c_smooth, 'Limegreen', lw=2, label='Savitzki-Golay filtered')
#    plot( xx, c-c_smooth, 'rs', ms=5, label='residual (noise)')
#    xlim(-0.8,0.8)
#    legend(loc=9, numpoints=1)
"""

