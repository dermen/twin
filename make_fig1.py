import PIL
# spacing
from pylab import *


import os
import itertools
from scipy.interpolate import interp1d

import matplotlib as mpl

def norm_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma.min()
    cmax = cor_vec_ma.max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n / crange


style.use('ggplot')
color_cycle = itertools.cycle( rcParams['axes.color_cycle'] )
red = color_cycle.next()
blue = color_cycle.next()
purple = color_cycle.next()
black = color_cycle.next()
yellow = color_cycle.next()
green = color_cycle.next()
pink = color_cycle.next()


#import matplotlib.gridspec as gridspec
#gs = gridspec.GridSpec(3, 2)
#ax = plt.subplot(gs[0, :])

##########
# LOAD THE DATA
cpsi_cub, cm_cub = np.load('cor_cuboctahedron.npy')
cpsi_cub = cpsi_cub[10:180]
cm_cub = cm_cub[10:180]

cpsi_tet, cm_tet = np.load('cor_twins.npy')
cpsi_tet = cpsi_tet[10:180]
cm_tet = cm_tet[10:180]
#############
###


from PIL import ImageOps

fig1 = figure(figsize=(6,6), frameon=False)

im1 =  PIL.Image.open('fig1_top_v7_filled.png')
#im1 = im1.resize( (im1.size[0]+100, im2.size[1] ) )
#im2 =  PIL.Image.open('fig1q_insertV3.png')


#delta = 20
lx = im1.size[0] + 500
ly = im1.size[1]
layer = PIL.Image.new('RGB', (lx,ly), (255,255,255))
layer.paste(im1, (0,0))

resize = 0.9
new_shape_im2= tuple( int(x*resize) for x in im2.size )
im2_new = im2.resize( new_shape_im2) 

lx2 = im1.size[0] + (im2.size[0] - im2_new.size[0])/2
ly2 = (im2.size[1] - im2_new.size[1])/2
layer.paste( im2_new, (lx2,ly2+delta) )

line_pos = (im1.size[0] + lx2 )/2
plot( ones(2)*(line_pos-delta/2), [0,layer.size[1]], 'k-', lw=2 )
xlim(0,layer.size[0] )
ylim(layer.size[1],0 )
imshow(layer)


ax = axes( [  ] )


lx = 1.4
ly = 2.2
#width of figure
font_size =12
# fig=plt.figure( 1, figsize=(fig_x,fig_y), dpi=DPI,frameon=False)

ax.axis('off')
purple = "#73006c"
yellow = "#fcff78"
orange = "#f7631a"
blue = "#8bbcf2"
gold ="#dfc63c"

lwf = 0.5
x = linspace(0,1,100)
# Ewald sphere
plot( x, sqrt(1-x**2), lw=2*lwf, alpha=0.5,color='k' , ls='--')
plot( x, -sqrt(1-x**2), lw=2*lwf,alpha=0.5, color='k', ls='--' )




arrowprops=dict(edgecolor='Black', arrowstyle = '<|-,head_length=0.4,head_width=0.3', 
            shrinkA = 0, shrinkB = 0, lw=lwf*2,color='Black')
annotate('', xy=(1,0), xytext=(.6,.8),arrowprops=arrowprops)
annotate('', xy=(1,0), xytext=(.6,-.8),arrowprops=arrowprops)

alphaVal=0.6
# original laser beam
arrowprops2=dict(edgecolor='k', linestyle='dashed',arrowstyle = '-', 
             shrinkA = 0, shrinkB = 0, lw=lwf*2,alpha = alphaVal)
annotate('', xy=(0,0), xytext=(1.0,0),arrowprops=arrowprops2)
## extend the original beam
#arrowprops2=dict(edgecolor='k', linestyle='dashed',arrowstyle = '-', 
#             shrinkA = 0, shrinkB = 0, lw=lwf*2,alpha = alphaVal)
#annotate('', xy=(1.0,0), xytext=(lx,0),arrowprops=arrowprops2)

# scattered laser beam
arrowprops2=dict(edgecolor='k', linestyle='dashed',arrowstyle = '-', 
             shrinkA = 0, shrinkB = 0, lw=lwf*2,alpha = alphaVal)
annotate('', xy=(0,0), xytext=(.6,.8),arrowprops=arrowprops2)


# label the two q's
s1=r"$q_1$"
t1 = text(.55,.6,s=s1, fontsize=font_size,  weight='bold')
s2 = r"$q_2$"
t2 = text(.55,-.6,s=s2, fontsize=font_size,  weight='bold')

# label psi
s3 = r"$\psi$"
t3 = text(.75,-0.04,s=s3, fontsize=font_size,  weight='bold')
t3.set_bbox({'color':yellow, 'alpha':0.6, 'edgecolor':"Black"})


# label theta{111}
s7 = r"$2\theta_{111}$"
t7 = text(0.25,0.1,s=s7, fontsize=font_size,  weight='bold')


#draw beam circle
circle_beam=Circle((0,0),1,color=yellow,alpha=0.5)
gca().add_artist(circle_beam)

#draw sample circle
sample_radius = .1
circle_sample=Circle((0,0),sample_radius,color=blue)
gca().add_artist(circle_sample)

#draw sample molecules
mol_radius = 0.01
positions = (1.-2.0*np.random.random_sample(size=[50,2]))*sample_radius*0.7
this_mol=Circle([0.0,0.0],mol_radius,color=gold,alpha=0.5,ec='k')
gca().add_artist(this_mol)
for this_pos in positions:
    this_mol=Circle(this_pos,mol_radius,color=gold,alpha=1,ec='k',linewidth=0.2)
    gca().add_artist(this_mol)
    
#savefig('test.png', dpi=DPI
#         ,transparent=True
#        ,bbox_inches='tight', pad_inches=0)





"""

ax = axes( [ 0, 0, .9, .9 ])
ax.imshow(im, 
        aspect='auto')
ax = gca()
ax.set_yticks([])
ax.set_xticks([])
ax.text(858,437,s=r'$2\theta$', fontsize=12)



fig2 = figure(figsize=( 5,5))

t = array( [-7/9., -5/9., -1/3., 0, 1/3., 5/9., 7/9.] )
tlab = [r'$-\frac{7}{9}$', r'$-\frac{5}{9}$', r'$-\frac{1}{3}$', '$0$', r'$\frac{1}{3}$', r'$\frac{5}{9}$', r'$\frac{7}{9}$']

a1 = axes([.2, .175, .29, .13 ])
a1.plot( cpsi_cub, norm_corr(cm_cub), '-', color='Limegreen', lw=2 )
a1.plot( cpsi_cub, norm_corr(cm_cub), 'ms', color=black, ms=.3 )
a1.set_ylim(-0.25,1.25)
a1.grid(1, ls='-', color=black, alpha=0.5)
a1.set_yticks([])
a1.set_xticks(t)
a1.set_xticklabels(tlab, fontsize=8)
a1.set_xlim(-0.1,1)
a1.tick_params(pad=0.1, length=0)
a1.set_xlabel(r'$\cos(\psi)$',fontsize=8, labelpad=2)

a2 = axes([.6, .175, .27, .13 ])
a2.plot( cpsi_tet, norm_corr(cm_tet), '-', color='Limegreen', lw=2 )
a2.plot( cpsi_tet, norm_corr(cm_tet), 'ms', color=black, ms=.3 )
a2.set_ylim(-0.25,1.25)
a2.grid(1, ls='-', color=black, alpha=0.5)
a2.set_yticks([])
a2.set_xticks(t)
a2.set_xticklabels(tlab, fontsize=8)
a2.set_xlim(-0.1,1)
a2.tick_params(pad=0.1, length=0)
a2.set_xlabel(r'$\cos(\psi)$',fontsize=8,labelpad=2)


#savefig('fig1_v7_5x5_real.pdf', format='pdf')

#os.chdir('/data/work/mender')
#pol = vstack( [h5py.File('cub%d.hdf5'%x,'r')['polar_intensities'].value[:,0]
#            for x in range( 10) ])
#f = h5py.File('cub%d.hdf5'%x,'r')
#k= f['k'].value
#wavelen = 2*pi / k
#th = arcsin( 2.668 * wavelen / (4*pi) )
#cpsi = cos(phis)*cos(th)**2 + sin(th)**2
#from popi.corr import correlate as Corr
#cor = zeros_like(cor)
#for i in xrange( pol.shape[0]):
#    cor[i] = Corr(pol[i], pol[i])

#os.chdir('/data/work/mender')
#pol = vstack( [h5py.File('tet%d.hdf5'%x,'r')['polar_intensities'].value[:,0]
#            for x in range( 10) ])
#f = h5py.File('tet%d.hdf5'%x,'r')
#k= f['k'].value
#wavelen = 2*pi / k
#th = arcsin( 2.668 * wavelen / (4*pi) )
#cpsi = cos(phis)*cos(th)**2 + sin(th)**2
#from popi.corr import correlate as Corr
#cor = zeros_like(cor)
#for i in xrange( pol.shape[0]):
#    cor[i] = Corr(pol[i], pol[i])
"""

