import numpy as np
import matplotlib.pyplot as plt
import pickle
plt.ion()

plt.rcParams['font.size']=20

# filename = 'Plots_Fig3/SFA_data.txt'
# px, pz, M2_neg_v, M2_pos_v, M2_neg_l, M2_pos_l = np.loadtxt(filename, unpack=True)

DumpFileName = 'OAM_Lin_DataDump'

def loadDataDump(DataKey, DumpFileName):
    '''Function to retrieve data from pickle dump'''
    DataDumpFile = open(DumpFileName,'rb')
    DataOut = pickle.load(DataDumpFile)[DataKey]
    DataDumpFile.close()
    return DataOut

[Fig2InsetData_v, Fig2Data_v] = loadDataDump('Fig_Paper_Velocity', DumpFileName)
[Fig2InsetData_l, Fig2Data_l] = loadDataDump('Fig_Paper_Length', DumpFileName)

dp=0.01
pzList = np.arange(-1.7, 1.7, dp)
pxList = np.arange(0, 1., dp)
pzGrid, pxGrid = np.meshgrid(pzList, pxList)

def plot(ZQs, norm=1, vmin=None, vmax=None, diff=False):
    fig,ax = plt.subplots(figsize=(4,6))
    if diff:
        mesh = ax.pcolormesh(pxGrid, pzGrid, ZQs, cmap='seismic', \
                             rasterized=True, shading='gouraud', vmin=vmin, vmax=vmax)
    else:
        mesh = ax.pcolormesh(pxGrid, pzGrid, ZQs/norm, cmap='inferno', \
                             rasterized=True, shading='gouraud', vmax=vmax)
    ax.set_aspect('equal')
    ax.set_ylim(-1.2,1.2)
    ax.set_xlim(0, 0.7)
    ax.set_xticks([0, 0.35, 0.7])
    ax.set_xticklabels(['0.0', '0.35', '0.7'])
    ax.set_yticks([-1.2, -0.6, 0, 0.6, 1.2])
    fig.colorbar(mesh, aspect=50)
    fig.tight_layout(pad=0.5)
    fig.tight_layout(pad=0.5)
    return fig


def gaussian(x,mu,sigma):
    return np.exp(-(x-mu)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2)
                  
def blur(f,sigma):
    l = len(f)
    f_ = np.zeros(l)
    for i in range(l):
        f_[i] = (gaussian(np.arange(l), i, sigma) * f).sum()
    return f_

# Load data
M2_length_plus1  = Fig2InsetData_l[0][0] 
M2_length_CD = -Fig2InsetData_l[1][0]
M2_length_minus1 = -(M2_length_CD - 2)/(M2_length_CD + 2) * M2_length_plus1
# M2_length_minus1 =  M2_length_plus1[:,::-1]
M2_length_tot = M2_length_plus1 + M2_length_minus1
# M2_length_CD = 2*(M2_length_plus1 - M2_length_minus1)/M2_length_tot
vmax_length = max(M2_length_tot.flatten())
M2_length_CD = np.where(M2_length_tot/2>0.02*vmax_length, M2_length_CD, np.nan)

M2_vel_plus1  = Fig2InsetData_v[0][0] 
M2_vel_CD = -Fig2InsetData_v[1][0]
M2_vel_minus1 = -(M2_vel_CD - 2)/(M2_vel_CD + 2) * M2_vel_plus1
# M2_vel_minus1 =  M2_vel_plus1[:,::-1]
M2_vel_tot = M2_vel_plus1 + M2_vel_minus1
# M2_vel_CD = 2*(M2_vel_plus1 - M2_vel_minus1)/M2_vel_tot
vmax_vel = max(M2_vel_tot.flatten())
M2_vel_CD = np.where(M2_vel_tot/2>0.02*vmax_vel, M2_vel_CD, np.nan)

# Plot figures
norm = max(M2_length_tot.flatten())
vmax = None
fig_length_tot = plot(M2_length_tot, norm=norm)
fig_length_plus = plot(M2_length_plus1, norm=norm)#, vmax=0.5)
fig_length_minus = plot(M2_length_minus1, norm=norm)#, vmax=0.5)
fig_length_CD = plot(M2_length_CD, diff=True)

norm = max(M2_vel_tot.flatten())
fig_vel_tot = plot(M2_vel_tot, norm=norm)
fig_vel_plus = plot(M2_vel_plus1, norm=norm)#, vmax=0.9)
fig_vel_minus = plot(M2_vel_minus1, norm=norm)#, vmax=0.9)
fig_vel_CD = plot(M2_vel_CD, diff=True, vmin=-2, vmax=2)

# # Export figures
fig_length_tot.savefig('Supplemental_Fig2(a)_length.png')
fig_length_plus.savefig('Supplemental_Fig2(b)_length_plus.png')
fig_length_minus.savefig('Supplemental_Fig2(b)_length_minus.png')
fig_length_CD.savefig('Supplemental_Fig2(c)_length.png')

fig_vel_tot.savefig('Supplemental_Fig2(a)_vel.png')
fig_vel_plus.savefig('Supplemental_Fig2(b)_vel_plus.png')
fig_vel_minus.savefig('Supplemental_Fig2(b)_vel_minus.png')
fig_vel_CD.savefig('Supplemental_Fig2(c)_vel.png')

# fig,ax = plt.subplots()
# ax.plot(M2_length_plus1[])

# norm = max(M2_length_tot.flatten())
# i = 25
# M2_length_plus1_blur = blur(M2_length_plus1[i],10)
# M2_length_minus1_blur = blur(M2_length_minus1[i],10)
# M2_length_tot_blur = M2_length_plus1_blur + M2_length_minus1_blur
# M2_length_CD_blur = 2*(M2_length_plus1_blur - M2_length_minus1_blur)/M2_length_tot_blur
# M2_length_CD_blur = np.where(M2_length_tot_blur/2>0.02*vmax_length, M2_length_CD_blur, np.nan)

# fig,ax = plt.subplots()
# l, = ax.plot(pzGrid[25], M2_length_plus1[25]/norm, \
#         label=r'$|M_{+1}(p_{\parallel},p_{\perp}=0.25)|^2$')
# ax.plot(pzGrid[25], M2_length_plus1_blur/norm, '--', color=l.get_color())
# l, = ax.plot(pzGrid[25], M2_length_minus1[25]/norm,\
#         label=r'$|M_{-1}(p_{\parallel},p_{\perp}=0.25)|^2$')
# ax.plot(pzGrid[25], M2_length_minus1_blur/norm, '--', color=l.get_color())
# l, = ax.plot(pzGrid[25], M2_length_CD[25],\
#         label=r'$\mathrm{PEVD}_{1}(p_{\parallel},p_{\perp}=0.25)$')
# ax.plot(pzGrid[25], M2_length_CD_blur, '--', color=l.get_color())

# ax.set_xlabel(r'$p_{\perp}$ (a.u.)')
# ax.legend(loc='best')
# fig.tight_layout(pad=0.5)
