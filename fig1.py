from pylab import *
import PIL
import itertools

def norm_corr( cor_vec):
    cor_vec_ma = ma.masked_invalid(cor_vec)
    cmin = cor_vec_ma.min()
    cmax = cor_vec_ma.max()
    crange = cmax-cmin
    cor_vec_n = cor_vec_ma - cmin
    return cor_vec_n / crange


purple = "#73006c"
yellow = "#fcff78"
orange = "#f7631a"
blue = "#8bbcf2"
gold ="#dfc63c"

lwf = 0.5
x = linspace(0,1,100)

font_size=12

extent = 600
gap=50.
im1 =  PIL.Image.open('tmp_fig1A_top.png')
im1_old_size = (1839, 863)
im1 = im1.resize( im1_old_size)
lx = im1.size[0] + extent
ly = im1.size[1]
layer = PIL.Image.new('RGB', (lx,ly), (255,255,255))
layer.paste(im1, (0,0))

fig = figure(figsize=(6,6), frameon=False)

ax1 = axes([0,0.6,1,0.4  ])


ax1.axis('off')
ax1.imshow( layer, aspect='auto')
ax1.text(212, 116, 'A', color='r', fontsize=font_size)


### MAKE axis 2 (insert axis)
val = (layer.size[0]-extent+gap)/layer.size[0]
ax2 = axes([ val  ,0.6,1-val, 0.4],frameon=False)
ax2.set_ylim(-1.1,1.1)
ax2.set_xlim(0,1.35)
ax2.axis('off')


# Ewald sphere
#ax2.axis('off')
ax2.plot( x, sqrt(1-x**2), lw=2*lwf, alpha=0.5,color='k' , ls='--')
ax2.plot( x, -sqrt(1-x**2), lw=2*lwf,alpha=0.5, color='k', ls='--' )
#draw beam circle
circle_beam=Circle((0,0),1,color=yellow,alpha=0.2)
ax2.add_artist(circle_beam)

arrowprops=dict(edgecolor='Black',
            arrowstyle = '<|-,head_length=0.4,head_width=0.3',
            shrinkA = 0, shrinkB = 0, 
            lw=lwf*2,
            color='Black')

ax2.annotate('', xy=(1,0), xytext=(.6,.8),
            arrowprops=arrowprops)
ax2.annotate('', xy=(1,0), xytext=(.6,-.8),
            arrowprops=arrowprops)

alphaVal=0.6
# original laser beam
arrowprops2=dict(edgecolor=orange, 
                linestyle='solid',arrowstyle = '-', 
             shrinkA = 0, shrinkB = 0, 
                lw=lwf*2,alpha = 0.3) 
ax2.annotate('', xy=(0,0), xytext=(1.0,0),
        arrowprops=arrowprops2)
## extend the original beam
plot( linspace(1,1.25,1000), zeros(1000), 'k--', alpha=0.5)

# scattered laser beam
ax2.annotate('', xy=(0,0), xytext=(.6,.8),
            arrowprops=arrowprops2)


# label the two q's
s1=r"$q_1$"
t1 = ax2.text(.55,.45,s=s1, fontsize=font_size,  
        weight='bold')
s2 = r"$q_2$"
t2 = ax2.text(.5,-.61,s=s2, fontsize=font_size,  
        weight='bold')

# label psi
t_psi = linspace( pi-arctan2( 2,1.), pi+arctan2( 2.,1.),1000 )
x_psi = 1+ 0.1* cos(t_psi)
y_psi = 0.1* sin(t_psi)
ax2.plot( x_psi, y_psi, 'k')

s3 = r"$\psi_{\max}$"#< \pi$"
ax2.annotate(s3, xy=(1-.1,0), xytext=(0,-.42),
            arrowprops=dict(edgecolor='k',
                arrowstyle = '->'), weight='bold')

# label 2 theta{111}
t_psi = linspace(0, arctan2( 4,3.), 1000 )
x_psi = 0.1* cos(t_psi)
y_psi = 0.1* sin(t_psi)
ax2.plot( x_psi, y_psi, 'k')
s7 = r"$2\theta_{111}$"
t7 = text(0.15,0.06,s=s7, fontsize=font_size,  weight='bold')

t_psi = linspace(0, pi - arctan2( 2.,1), 1000, endpoint=False )
#x_psi = 1+0.14* cos(t_psi)
#y_psi = 0.14* sin(t_psi)
#ax2.plot( x_psi, y_psi, 'k')
x_psi = 1+0.16* cos(t_psi)
y_psi = 0.16* sin(t_psi)
ax2.plot( x_psi, y_psi, 'k')

ax2.annotate(r'$\frac{\pi}{2} + \theta_{111}$', xytext=(.8,.8), 
                xy=(1.07,.17),
                arrowprops=dict(edgecolor='k',
                arrowstyle = '->'), weight='bold')

#label phi1
ax1.text(1508,280,s=r"$\phi_1$", fontsize=font_size,  
        weight='bold', color='white')
# label phi2
ax1.text(1410,200+layer.size[1]/2.,s=r"$\phi_2$", fontsize=font_size,  
        weight='bold', color='white')

ax1.text(1520,420, s=r'$\Delta=\pi$', color='white', fontsize=font_size)

det_center = (1500, 412)
#phi1 = (1676,323)
#phi2 = (1686,624)
phi1 = (1492,161)
phi2 = (1505,682)
ax1.annotate('',xy=det_center, xytext=phi1,
                arrowprops=dict(edgecolor='white',
                arrowstyle = '<-', lw=1))
ax1.annotate('',xy=det_center, xytext=phi2,
                arrowprops=dict(edgecolor='white',
                arrowstyle = '<-', lw=1))


# draw 2theta
samp_cent = (644,480)
ax1.annotate('',xy=samp_cent, xytext=(853,315),
                arrowprops=dict(edgecolor='black',
                arrowstyle = '-'))

ax1.annotate('',xy=samp_cent, xytext=(868,448),
                arrowprops=dict(edgecolor='black',
                arrowstyle = '-'))
# label 2theta
ax1.text(772,438,s=r"$2\theta$", fontsize=font_size,  
        weight='bold', color='black')

#ax1.add_artist( Circle(phi1,10,color='white')  )
#ax1.add_artist( Circle(phi2,10,color='white')  )
ax1.add_artist( Circle(det_center,8,color='white')  )

ax2.text(1.1, -.8, 'B', color='r', fontsize=font_size)

im3 = PIL.Image.open( 'pymol_cubo_highest.png')

ax3 = axes( [0.05,0.1,0.45,0.5] )
#ax3.axis('off')
#ax3.imshow( im3, aspect='auto')


im4 = PIL.Image.open( 'pymol_deca_highres_wline.png')
ax4 = axes( [.55,0.1,0.45,0.5] )
#ax4.axis('off')
#ax4.imshow( im4, aspect='auto')


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



t = array( [-7/9., -5/9., -1/3., 0, 1/3., 5/9., 7/9.] )
tlab = [r'$-\frac{7}{9}$', r'$-\frac{5}{9}$', r'$-\frac{1}{3}$', '$0$', r'$\frac{1}{3}$', r'$\frac{5}{9}$', r'$\frac{7}{9}$']

ax3.plot( cpsi_cub, norm_corr(cm_cub), '-', color='Limegreen', lw=2 )
ax3.plot( cpsi_cub, norm_corr(cm_cub), 'ms', color='k', ms=1 )
ax3.set_ylim(-0.25,1.25)
ax3.grid(1, ls='-', color='k', alpha=0.5)
ax3.set_yticks([])
ax3.set_xticks(t)
ax3.set_xticklabels(tlab, fontsize=font_size)
ax3.set_xlim(-0.1,1)
ax3.tick_params(pad=7, length=0)
#ax3.set_xlabel(r'$\cos\, \psi$',fontsize=font_size,labelpad=0)

ax4.plot( cpsi_tet, norm_corr(cm_tet), '-', color='Limegreen', lw=2 )
ax4.plot( cpsi_tet, norm_corr(cm_tet), 'ms', color='k', ms=1 )
ax4.set_ylim(-0.25,1.25)
ax4.grid(1, ls='-', color='k', alpha=0.5)
ax4.set_yticks([])
ax4.set_xticks(t)
ax4.set_xticklabels(tlab, fontsize=font_size)
ax4.set_xlim(-0.1,1)
ax4.tick_params(pad=7, length=0 )
#ax4.set_xlabel(r'$\cos\,\psi$',fontsize=font_size,labelpad=0)



ax5 = axes( [.27, .37, .2, .2] )
ax5.set_xticks([])
ax5.set_yticks([])
ax5.imshow( im3,aspect='auto')

ax6 = axes( [.27+.5, .37, .2, .2]) 
ax6.set_xticks([])
ax6.set_yticks([])
ax6.imshow(im4,aspect='auto')


ax3.spines['top'].set_visible(0)
ax3.spines['right'].set_visible(0)
ax3.spines['left'].set_visible(0)
ax3.spines['bottom'].set_visible(0)

ax4.spines['top'].set_visible(0)
ax4.spines['right'].set_visible(0)
ax4.spines['left'].set_visible(0)
ax4.spines['bottom'].set_visible(0)

ax3.set_xlim(.3333-0.22,1)
ax4.set_xlim(.3333-0.22,1)

fig.text(.49, .02, s=r'$\cos\, \psi$',fontsize=font_size, weight='bold' )

ax3.set_ylabel(r'$C(\cos\,\psi)$', fontsize=font_size,labelpad=0 )
ax3.text( 0.2, 1.1 ,'A',  color='r', fontsize=font_size)
ax4.text( 0.2, 1.1 ,'B',  color='r', fontsize=font_size)

ax4.annotate('', xy=(-0.05, 0), xycoords='axes fraction', 
            xytext=(-0.05, 1), 
            arrowprops=dict(arrowstyle="-", ls='dashed', 
                color='dimgrey'))

ax2.annotate('', xy=(-0.09, -1), xycoords='axes fraction', 
            xytext=(-0.09, 1), 
            arrowprops=dict(arrowstyle="-", ls='solid', 
                color='dimgrey'))


savefig('Fig1Fig2.png', dpi=500, facecolor='w')
#ax4.set_ylabel(r'$C(\cos(\psi))$', fontsize=font_size ,labelpad=0)
# FORGO the sample circle for now
#draw sample circle
#sample_radius = .1
#circle_sample=Circle((0,0),sample_radius,color=blue)
#gca().add_artist(circle_sample)

#draw sample molecules
#mol_radius = 0.01
#positions = (1.-2.0*np.random.random_sample(size=[50,2]))*\
#                sample_radius*0.7
#this_mol=Circle([0.0,0.0],mol_radius,
#            color=gold,alpha=0.5,ec='k')
#gca().add_artist(this_mol)
#for this_pos in positions:
#    this_mol=Circle(this_pos,mol_radius,
#        color=gold,alpha=1,ec='k',
#        linewidth=0.2)
#    gca().add_artist(this_mol)
    


#savefig('tmp.png', dpi=200)
#ax2.imshow()

#####
# MAKE axis 3,4 for particle models



