import matplotlib.pyplot as plt


charm_energies = [0.3879999998751531, 0.8420076252175448, 1.0371822051751995]
bottom_energies = [0.09086700922681423, 0.47327644877899033, 0.58994433097599]
charm_aligned = [0.04326747162164409, 3.130187115216843e-28, 0.03013194434398317]
charm_antialigned = [-0.12980241486493227, -9.390561345650528e-28, -0.09039583303194951]
bottom_aligned = [0.016984157921053088, 2.4605924185003373e-28, 0.009343895479365883]
bottom_antialigned = [-0.05095247376315926, -7.381777255501012e-28, -0.02803168643809765]

charm_full = []
bottom_full = []

for i in range(0,3):
    charmlow = charm_energies[i] + charm_antialigned[i] 
    charm_full.append(charmlow)
    charmhigh = charm_energies[i] + charm_aligned[i]
    charm_full.append(charmhigh)
    
    bottomlow = bottom_energies[i] + bottom_antialigned[i]
    bottom_full.append(bottomlow)
    bottomhigh = bottom_energies[i] + bottom_aligned[i]
    bottom_full.append(bottomhigh)
    
print(charm_full)
print(bottom_full)
plt.figure(1)
xsc = [1]*len(charm_full)
xsb = [1.6]*len(bottom_full)
plt.scatter(xsc, charm_full, s=9000, marker="_", linewidth=2, zorder=3, color = 'red', label = 'charmonium')
plt.scatter(xsb, bottom_full, s=9000, marker="_", linewidth=2, zorder=3, color = 'blue', label = 'bottomonium')
plt.xlim(0.8,1.8)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) 
plt.ylabel('Binding energy (GeV)',fontsize = 13)
#plt.legend(loc = 'lower left', s = 10)
labels = ['$1^{1}S_{0}$', '$1^{3}S_{1}$','$l=1$','','$2^{1}S_{0}$','$2^{3}S_{1}$']
for i in range(0,len(labels)-1):
    plt.text(1.16,charm_full[i],labels[i],color = 'red')
    plt.text(1.35,bottom_full[i],labels[i],color = 'blue')
plt.text(1.16,charm_full[-1],labels[-1], color = 'red')
plt.text(1.35,bottom_full[-1]+0.03,labels[-1],color = 'blue')
plt.text(0.87,-0.1,'Charmonium',color = 'red',fontsize = 14)
plt.text(1.47,-0.1,'Bottomonium',color = 'blue',fontsize = 14)
plt.savefig('energy levels.pdf',bbox_inches="tight")
plt.show()

labels_p = ['$1^{1}P_{1}$','$1^{3}P_{0}$','$1^{3}P_{1}$','$1^{3}P_{2}$']
fig, (ax1, ax2) = plt.subplots(1,2)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) 
Charm_lsplits = [0.842007637,
0.845899197,
0.838116077,
0.849790757]
Bottom_lsplits = [0.473276449,0.477118242,0.469434656,0.480960035]
xsc =[1]*len(Charm_lsplits)
xsb = [1]*len(Bottom_lsplits)
ax1.scatter(xsc, Charm_lsplits, s=9000, marker="_", linewidth=2, zorder=3, color = 'red', label = 'charmonium')
ax1.set(ylabel = 'Binding Energy (GeV)', xlabel = 'Charmonium')
ax1.xaxis.label.set_color('red')
ax2.scatter(xsb, Bottom_lsplits, s=9000, marker="_", linewidth=2, zorder=3, color = 'blue', label = 'charmonium')
ax2.yaxis.tick_right()
ax2.set(ylabel = 'Binding Energy (GeV)', xlabel = 'Bottomonium')
ax2.xaxis.label.set_color('blue')
for i in range(0,len(labels_p)):
    ax1.text(1.036,Charm_lsplits[i]-0.0001,labels_p[i],color = 'red')
    ax2.text(1.036,Bottom_lsplits[i]-0.0001,labels_p[i],color = 'blue')
for ax in (ax1,ax2):
    ax.set_xticks([])
plt.savefig('l splits.pdf',bbox_inches="tight")