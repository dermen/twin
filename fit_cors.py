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

plot( xx, cnn, 'kx')
plot( xx, fit_data, 'b', lw=3 )
plot( xx[local_max], fit_data[local_max], 's',color='Limegreen', ms=8)
plot( xx[local_mins], fit_data[local_mins], 'd',color='r', ms=8)




#aksjdtkjahfl
############################
# GAUSSIAN SUM FIT CODE    #
############################
def sum_gauss( params, xdata, ydata , Ngauss, ctrs):
    model = zeros_like( xdata)
    for i in xrange( Ngauss):
        amp = params['amp%d'%i ].value
        wid = params['wid%d'%i ].value
        off = params['off%d'%i ].value
        model = model+  amp * np.exp( -((xdata - ctrs[i])/wid)**2) + off
    return model - ydata

#W = 0.01252791440815894
#params = Parameters()
#for i,pk in enumerate(local_max):
#    params.add('off%d'%i, value= 0.1, min=0, max=2)
#    params.add('wid%d'%i, value= W, min=W-W*0.7, max=W+W*0.7)
#    params.add('amp%d'%i, value= 1 , min=0, max=2)
#result_gauss = minimize(sum_gauss, params, 
#            args=(xx, cnn, N_peaks, xx[local_max] ))

####
#alslsakdlaskd
colors = cycle(  ['r', 'b', 'k'] )
colors_fit = cycle(  ['Limegreen', 'Darkorange', 'DeepPink'] )
local_max_split = np.array_split( local_max, 10)

all_results_gauss = []
fit = zeros_like( xx)

#for ii,cnn in enumerate(cnn_series):
#subplot(3,2,ii+1)
#xlim(-0.79, 0)
#title('N=%d'%shot_series[ii])
for local_max_ in local_max_split:
    lowest_max = local_max_[0]
    highest_max = local_max_[-1]
    lower_bound = local_mins[ local_mins < lowest_max][-1]
    upper_bound = local_mins[ local_mins > highest_max][0]

    ctrs = xx[ local_max_]
    xdata = xx[ np.arange(lower_bound, upper_bound +1)]
    ydata = cnn[ np.arange(lower_bound, upper_bound +1)]
   
    #figure(1)
    #plot( xdata, ydata , 'd',ms=9, color=colors.next())
    
    W = 0.01252791440815894
    N_gauss = len(local_max_)
    params = Parameters()
    for i_ in xrange(N_gauss):
        params.add('off%d'%i_, value= 0.1, min=0, max=2)
        params.add('wid%d'%i_, value= W, min=W-W*0.7, max=W+W*0.7)
        params.add('amp%d'%i_, value= 1 , min=0, max=2)
    
    result_gauss = minimize(sum_gauss, params, 
                args=(xdata, ydata, N_gauss, ctrs ))

    all_results_gauss.append( [result_gauss, ctrs])
    
    rp = result_gauss.params
    offs = zeros_like(ctrs)
    amps = zeros_like(ctrs)
    wids = zeros_like(ctrs)
    for i_ in xrange( N_gauss):
        offs[i_] = rp['off%d'%i_].value
        amps[i_] = rp['amp%d'%i_].value
        wids[i_] = rp['wid%d'%i_].value

    fit_partial = sum_gauss( result_gauss.params, xdata, zeros_like( xdata),
                                N_gauss, ctrs )
    
    fit_full = sum_gauss( result_gauss.params, xx, zeros_like( xx),
                             N_gauss, ctrs )
    
    fit += fit_full - sum(offs)

    #c = colors_fit.next()
    #plot( xdata, ydata, 'bs', ms=3, mec='b' )
    #plot( xdata, fit_partial , color='Limegreen', lw=2)
    #plot( xdata, ydata - fit_partial , 's', color='k', mec='k',ms=3)
    #plot( xdata, fit_partial-sum(offs), color='Darkorange', lw=2)

#np.save('data_fit', [xx,fit] )

shot_series = array([ 50*i for i in xrange(1,500) ])
cnn_series = [ [norm_corr( C[np.random.randint(0,C.shape[0], 
                s ) ].mean(0) ).data for s in shot_series]
                for jj in range(3) ]
np.save( 'cnns_50', [shot_series, cnns] )


