import numpy as np
from scipy.integrate import odeint, simpson
import matplotlib.pyplot as plt
#from scipy.special import sph_harm
#from matplotlib import cm, colors
#from mpl_toolkits import mplot3d

alpha = 1/137
def coulomb(r):
        return -alpha/r

def schrodinger(y,r, potential, E, l, mu):
    u,v = y
    dydr = [v, (((l*(l+1))/r**2)-2*mu*(E-potential(r)))*u]
    return dydr

#Natural units in MeV

rs = np.linspace(0.00000001,14, 1500)
E = -1.304569*10**-5
l = 0
mu = 0.51099895000
y0 = [0,1]
a0 = 1/(mu*alpha)



# solution = odeint(schrodinger, y0,rs, (coulomb, E, l, mu)) #Gives un-normalised wavefunction
# u = solution[:,0]
# plt.plot(rs/a0,u)
# plt.show()
# norm = simpson(u*u, x = rs)
# unorm = u/np.sqrt(abs(norm))

# pdist = unorm*unorm
#print(unorm[-1])
# plt.plot(rs/a0,pdist)
# plt.title('ground state')
# plt.show()
'Now have a method for solving the schrodinger equation for a given energy,'
'just need to implement a modified bisection method'#


def getEnergy(E0, E1, potential, tolerance, l, mu):
    # E0 = E0*10**-6
    # E1 = E1*10**-6
    
    t1 = odeint(schrodinger, y0, rs, (potential, E1, l, mu))[:,0]
    t0 = odeint(schrodinger, y0, rs, (potential, E0, l, mu))[:,0]
    #plt.plot(rs, t0, label = 'T0')
    #plt.plot(rs, t1, label = 'T1')
    #plt.legend()
    #plt.show()
    norm1 = np.abs(simpson(t1*t1, x = rs))
    norm0 = np.abs(simpson(t0*t0, x = rs))
    t0 = t0/np.sqrt(norm0)
    t1 = t1/np.sqrt(norm1)
    test_0 = t0[-1]
    test_1 = t1[-1]
    error = 'No sign change on interval'
    #print('test_0: ' + str(test_0))
    #print('test_1 : ' + str(test_1))
    if test_0*test_1 >0:
        return error
    
    if test_0 < 0:
        Elow = E0
        Ehigh = E1
    else:
        Ehigh  = E0
        Elow = E1
    Emid = (Ehigh+Elow)*0.5
    #overflow_check = np.max(t1)
    #print(abs(test_0 - test_1))
    count = 0
    while abs(Elow - Ehigh) > tolerance or count<60:
        
        count = count+1
        tlow = odeint(schrodinger, y0, rs, (potential, Elow, l, mu))[:,0]
        normlow = np.abs(simpson(tlow*tlow, x = rs[::-1]))
        tlow = tlow/np.sqrt(normlow)
        test_0 = tlow[-1]

        thigh = odeint(schrodinger, y0, rs, (potential, Ehigh, l, mu))[:,0]
        normhigh = np.abs(simpson(thigh*thigh, x = rs[::-1]))
        #print(Emid)
        thigh = thigh/np.sqrt(normhigh)
        test_1 = thigh[-1]
        
       
        tmid = odeint(schrodinger, y0, rs,(potential, Emid, l, mu))[:,0]
        normmid = np.abs(simpson(tmid*tmid, x = rs[::-1]))
        tmid = tmid/np.sqrt(normmid)
        
        testmid = tmid[-1]
        #plt.plot(rs/a0, tmid)
        #plt.title('Get smeared, idiot')
        
        if testmid < 0:
            Elow = Emid
        
        else:
            Ehigh = Emid
        Emid = (Ehigh + Elow)*0.5
    #overflow_check = np.max(tmid)
    #print(overflow_check)
    #plt.show()
    #if overflow_check > 1e7:
     #   return 'likely overflow'
    #else:
    return Emid

    
def milestone(ns, potential, energy_start, ls, mu):
    Es = []
    i = 0
    for n in ns:
        E_init = energy_start[i]
        E0 = E_init + 0.1
        E1 = E_init - 0.1
        #print(E0, E1)
        l = ls[i]
        i = i+1
        
        Enew = getEnergy(E0, E1, potential, 0.0001, l, mu)
        #print('Enew = ' + str(Enew))
        Es.append(Enew)
        #error = (E_anal- Enew)/E_anal
        #print('Energy error: ' + str(error) + '%')
        # u = odeint(schrodinger, y0, rs, (coulomb, E, l, mu))[:,0]
        # unorm = np.abs(simpson(u,x = rs))
        # u = u/np.sqrt(unorm)
        
    return Es
            
ns = [1,2]

# E_list = milestone(ns)           
# Energies = E_list[0]
# ns = E_list[1]
# ls = E_list[2]
# colours = ['#2f8dff', 'blue', '#682860']
'With a method to obtain energy eigenvalues, can now plot the probability distributions'

# plt.figure(figsize = (11,7))

# for i in range(0,3):
#     E = Energies[i]*10**-6
#     l = ls[i]
#     n = ns[i]
#     U = odeint(schrodinger, y0, rs, (coulomb, E,l,mu))[:,0]
#     norm = np.abs(simpson(U*U,x = rs))
#     U = U/np.sqrt(norm)
#     Prob = U*U
#     plt.plot(rs/a0, Prob, label = '(n,l) = ' + str((n,l)), color = colours[i])


# plt.xlim(0,max(rs/a0)+0.5)
# plt.xlabel(r'$\frac{r}{a_{0}}$', fontsize = '25')
# plt.xticks(fontsize = '20')

# plt.ylim(0, 0.00215)
# plt.ylabel('Radial Probability density', fontsize = '22')
# plt.yticks(fontsize = '17')

# plt.legend(fontsize = '17')
# plt.show()

alpha_s_charm = 0.40
alpha_s_bottom = 0.28

'First will use Charmonium measurements to iterate to a value of beta'

E0_charm = 0.388 #ground state charmonium energy, from this point all measured in GeV
mu_charm = (1.34)/(2)
mu_bottom = 4.70/2
def get_beta(B0, B1, alpha, tolerance, l, mu, E0):
    def cornell0(r):
        return (-4*alpha)/(3*r) + B0*r 
    
    def cornell1(r):
        return (-4*alpha)/(3*r) + B1*r 
    
    #print((((l*(l+1))/rs**2)-2*mu*(E-cornell1(rs))))
    t1 = odeint(schrodinger, y0, rs, (cornell1, E0, l, mu))[:,0]
    t0 = odeint(schrodinger, y0, rs, (cornell0, E0, l, mu))[:,0]
    norm1 = np.abs(simpson(t1*t1, x = rs))
    norm0 = np.abs(simpson(t0*t0, x = rs))
    t0 = t0/np.sqrt(norm0)
    t1 = t1/np.sqrt(norm1)
    test_0 = t0[-1]
    test_1 = t1[-1]
    error = 'No sign change on interval'
    if test_0*test_1 >0:
        return error
    
    if test_0 < 0:
        Blow = B0
        Bhigh = B1
    else:
        Bhigh  = B0
        Blow = B1
    Bmid = (Bhigh+Blow)*0.5
    runs = 0
    while abs(test_0-test_1) > tolerance:
        def cornell_low(r):
            return (-4*alpha)/(3*r) + Blow*r
        
        def cornell_mid(r):
            return (-4*alpha)/(3*r) + Bmid*r
        
        def cornell_high(r):
            return (-4*alpha)/(3*r) + Bhigh*r
        
        tlow = odeint(schrodinger, y0, rs, (cornell_low, E0_charm, l, mu))[:,0]
        normlow = np.abs(simpson(tlow*tlow, x = rs[::-1]))
        tlow = tlow/np.sqrt(normlow)
        test_0 = tlow[-1]

        thigh = odeint(schrodinger, y0, rs, (cornell_high, E0_charm, l, mu))[:,0]
        normhigh = np.abs(simpson(thigh*thigh, x = rs[::-1]))
        #print(Bmid)
        thigh = thigh/np.sqrt(normhigh)
        test_1 = thigh[-1]
        
        #print(test_1)
        tmid = odeint(schrodinger, y0, rs,(cornell_mid, E0_charm, l, mu))[:,0]
        normmid = np.abs(simpson(tmid*tmid, x = rs[::-1]))
    
        tmid = tmid/np.sqrt(normmid)
        
        testmid = tmid[-1]
        
        if testmid < 0:
            Blow = Bmid
        
        else:
            Bhigh = Bmid
        Bmid = (Bhigh + Blow)*0.5
        #print(Bmid)
        runs = runs+1
    print("iterated " + str(runs) +" " +'times')
    plt.figure(1)
    plt.plot(rs,tmid)
    plt.xlabel('r $(GeV^{-1})$')
    plt.ylabel('u(r)')
    plt.title('Radial wavefunction for returned Beta')
    plt.show()
    return Bmid

B_lit_charm = 0.195

B0 = 0.185
B1 = 0.205

B_lit_bottom = 0.168
beta = get_beta(B0,B1,alpha_s_charm, 0.001, 0, mu_charm, E0_charm)
print('beta: ' + str(beta))

def cornell_charm(r):
    return (-4*alpha_s_charm)/(3*r) + beta*r 

def cornell_bottom(r):
    return (-4*alpha_s_bottom)/(3*r) + beta*r

#charm_ground = odeint(schrodinger, y0, rs, (cornell_charm, E0_charm, 0, mu))
#charm_ground_U = charm_ground[:,0]
#norm = np.sqrt(simpson(charm_ground_U, x = rs))

#charm_gr'#ound_deriv = charm_ground[:,1]/np.sqrt(norm)
#R_zero = charm_ground_deriv[0]

#spin_coupling = ((8*alpha_s_charm)/(9*1.34**2))*np.abs(R_zero)**2

V0_charm = np.sqrt((4*alpha_s_charm)/(3*beta))

accessible_cols = ['red','blue','black']
Energies = milestone([1,1,2], cornell_charm, [E0_charm,0.9, 1.1], [0,1,0], mu_charm)  
Energies_bottom = milestone([1,1,2], cornell_bottom, [0.01, 0.4,0.6], [0,1,0], mu_bottom)   
ls = [0,1,0]
ns = [1,1,2]

print('charm energies: ' + str(Energies))
print('Bottom energies: ' + str(Energies_bottom))

plt.figure(2, figsize = (11,7), frameon=False)
plt.axhline(0,0,14, linestyle = 'dashdot', color = 'grey')
for i in range(0,3):
    n = ns[i]
    l = ls[i]
    u = odeint(schrodinger, y0, rs, (cornell_charm, Energies[i], l, mu_charm))[:,0]
    unorm = np.abs(simpson(u*u,x = rs[::-1]))
    u = u/np.sqrt(unorm)
    prob = u*u
    plt.plot(rs, u,label = '(n,l) = ' + str((n,l)), color = accessible_cols[i])    
plt.legend(fontsize = '17')

plt.xlim(0,14)
plt.xlabel('Quark seperation $(GeV^{-1})$', fontsize = '22')
plt.xticks(fontsize = '20')

plt.ylabel('u(r) $(Gev^{1/2})$', fontsize = '22')
#plt.ylim(0,0.45)
plt.yticks(fontsize = '17')
#for pos in ['right', 'top']: 
    #plt.gca().spines[pos].set_visible(False) 
#plt.show() 

#plt.vlines(V0_charm, 0,0.68,)
plt.savefig('Charmonium wavefunctions.pdf',bbox_inches="tight")


'''if given odeint solution, returns the derivative at 0'''
def getderiv_0(solution): 
    #print("make sure you have input a raw odeint solution or this wont work")
    u = solution[:,0]
    unorm = np.abs(simpson(u*u,x = rs[::-1]))
    solution = solution/np.sqrt(unorm)
    deriv = solution[:,1]
    return deriv[0]

E_big = E0_charm + 0.1
E_small = E0_charm - 0.1

u0_big = odeint(schrodinger, y0, rs, (cornell_charm, E_big, 0, mu_charm))[:,0]
u0_big = u0_big/np.sqrt(np.abs(simpson(u0_big*u0_big, x = rs)))
u0_actual = odeint(schrodinger, y0, rs, (cornell_charm, E0_charm, 0, mu_charm))[:,0]
u0_actual = u0_actual/np.sqrt(np.abs(simpson(u0_actual*u0_actual, x = rs)))
u0_small = odeint(schrodinger, y0, rs, (cornell_charm, E_small, 0, mu_charm))[:,0]
u0_small = u0_small/np.sqrt(np.abs(simpson(u0_small*u0_small, x = rs)))

plt.figure(3, figsize = (10,7))
plt.plot(rs,u0_big, color = 'blue', label = 'Energy = 0.588 GeV')
plt.plot(rs,u0_actual, color = 'black', label = 'Energy = 0.388 GeV')
plt.plot(rs,u0_small, color = 'r', label = 'Energy = 0.188 GeV')
plt.xlim(0,14)
plt.legend(loc = 'upper left', fontsize = '14')
plt.xlabel('Quark seperation $(GeV^{-1})$', fontsize = '22')
plt.ylabel('u(r) $(Gev^{1/2})$',fontsize = '22')
plt.xticks(fontsize = '16')
plt.yticks(fontsize = '16')
plt.savefig('Bisection illustration.pdf',bbox_inches="tight")
#for pos in ['right', 'top']: 
#    plt.gca().spines[pos].set_visible(False) 
plt.show()


def spin_split(u,mass,alpha):
    R0 = getderiv_0(u)
    delta = ((8*alpha)/(9*mass*mass))*R0**2
    return delta

'''deltas give prefactor to spin dot product'''
charm_aligned = []
charm_antialigned = []
bottom_aligned = []
bottom_antialigned = []



for i in range(0,3):
    Ec = Energies[i]
    Eb = Energies_bottom[i]
    uc = odeint(schrodinger, y0, rs, (cornell_charm, Ec, [0,1,0][i], mu_charm))
    ub = odeint(schrodinger, y0, rs, (cornell_bottom, Eb, [0,1,0][i], mu_bottom))
    split_c = spin_split(uc,2*mu_charm,alpha_s_charm)
    split_b = spin_split(ub,2*mu_bottom,alpha_s_bottom)
    charm_aligned.append(split_c*0.25)
    charm_antialigned.append(split_c*-0.75)
    bottom_aligned.append(split_b*0.25)
    bottom_antialigned.append(split_b*-0.75)
    
print('charm aligned: ' + str(charm_aligned))
print('charm anti aligned:' + str(charm_antialigned))
print('bottom aligned: ' + str(bottom_aligned))
print('bottom antialigned: ' + str(bottom_antialigned))


'''Extra code for finding the LS splitting'''
corrections_charm = []
corrections_bottom = []

LS = [0,-2,-1,1]
DV_charm = (2*alpha_s_charm)/(2*1.34**2*rs**3) - beta/(2*1.34**2*rs)
DV_bottom = (2*alpha_s_bottom)/(2*1.34**2*rs**3) - beta/(2*1.34**2*rs)
U_state_charm = odeint(schrodinger, y0, rs, (cornell_charm, Energies[1], l, mu_charm))[:,0]
unorm = np.abs(simpson(U_state_charm**2,x = rs[::-1]))
U_state_charm = U_state_charm/np.sqrt(unorm)
U_state_bottom = odeint(schrodinger, y0, rs, (cornell_charm, Energies_bottom[1], l, mu_bottom))[:,0]
unorm = np.abs(simpson(U_state_bottom**2,x = rs[::-1]))
U_state_bottom = U_state_bottom/np.sqrt(unorm)

for i in LS:
    LS_state = LS[i]
    DV_c = DV_charm*LS_state
    Correction_charm  = simpson(U_state_charm**2*DV_c, x = rs)
    corrections_charm.append(Correction_charm)
    
    DV_b = DV_bottom*LS_state
    Correction_bottom = simpson(U_state_bottom**2*DV_c, x = rs)
    corrections_bottom.append(Correction_bottom)
print("Charm P corrections: " + str(corrections_charm))
print("Bottom P corrections: " + str(corrections_bottom))


energies_check = []

for i in range(0,100):
    E_it = milestone([1,1,2], cornell_charm, [E0_charm,0.9, 1.1], [0,1,0], mu_charm)  
    energies_check.append(E_it)
    
energies_check = np.array(energies_check)
ground = energies_check[:,0]
first = energies_check[:,1]
second = energies_check[:,2]
xs = [1]*len(ground)
checkg = np.std(ground)
checkf = np.std(first)
checks = np.std(second)

print('Ground state stdev '+str(checkg))
print('(1,1) state stdev ' +str(checkf))
print('(2,0) state stdev ' +str(checks))
