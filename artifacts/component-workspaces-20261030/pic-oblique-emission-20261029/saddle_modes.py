"""Derived orbital local-saddle diagnostics; never a native particle source."""
import math

qe=1.602176634e-19
me=9.1093837139e-31
kb=1.380649e-23
h=6.62607015e-34
hbar=h/(2*math.pi)
eps0=8.8541878188e-12
c=qe/(16*math.pi*eps0)
C=4*math.pi*me*qe**3/h**3
A0=4*math.pi*me*qe*kb**2/h**3

def matrix(rows,cols):return [[0.]*cols for _ in range(rows)]
def transpose(a):return [list(x) for x in zip(*a)]
def multiply(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def maximum(a):return max(abs(x) for row in a for x in row)
def difference(a,b):return [[x-y for x,y in zip(ar,br)] for ar,br in zip(a,b)]

J=matrix(6,6)
for i in range(3):J[i][i+3]=1;J[i+3][i]=-1

def modes(bn,bt):
    """Dimensionless B components omega_cn/omega0, omega_ct/omega0."""
    bc2=bn*bn+bt*bt
    a=1-bc2
    disc=math.hypot(a,2*bn)
    u2=(a+disc)/2 if a>=0 else 2*bn*bn/(disc-a)
    if bn==0:
        return {'unstable':math.sqrt(max(0.,u2)), 'stable':math.sqrt(max(0.,bc2-1)), 'basis':None}
    if u2<=0:raise ValueError('No numerically resolved unstable mode')
    s2=bn*bn/u2
    u,s=math.sqrt(u2),math.sqrt(s2)
    y=bn*bt/(u2+bn*bn)
    xs=-bt*s2/(bn*(1+s2))
    au,ass=math.hypot(1,y),math.hypot(1,xs)
    qu=[1/au,y/au,0,0,0,0]
    pu=[0,0,-bt/(au*(u2+bn*bn)),1/au,y/au,0]
    qs=[xs/ass,1/ass,0,0,0,0]
    ps=[0,0,-1/(bn*ass),xs/ass,1/ass,0]
    root=math.sqrt(abs(bn))
    gq=[0,1/root,0,0,0,-bn/root]
    gp=[0,0,math.copysign(1,bn)/root,0,0,0]
    S=transpose([qu,qs,gq,pu,ps,gp])
    H=matrix(6,6)
    H[0][0]=bt*bt-1;H[1][1]=bn*bn
    H[0][1]=H[1][0]=-bn*bt
    H[0][5]=H[5][0]=-bt
    H[1][5]=H[5][1]=bn
    for i in range(3,6):H[i][i]=1
    target=matrix(6,6)
    for i,x in enumerate([-u2,s2,0,1,1,0]):target[i][i]=x
    symplectic=maximum(difference(multiply(transpose(S),multiply(J,S)),J))
    diagonal=maximum(difference(multiply(transpose(S),multiply(H,S)),target))/max(1,maximum(target))
    # Independent mechanical-coordinate Lorentz generator, q=-e.
    T=matrix(6,6);Ti=matrix(6,6)
    for i in range(6):T[i][i]=Ti[i][i]=1
    T[5][0]=-bt;T[5][1]=bn;Ti[5][0]=bt;Ti[5][1]=-bn
    L=matrix(6,6)
    for i in range(3):L[i][i+3]=1
    L[3][0]=1;L[3][5]=bt;L[4][5]=-bn;L[5][3]=-bt;L[5][4]=bn
    lorentz=maximum(difference(multiply(T,multiply(J,multiply(H,Ti))),L))/max(1,maximum(L))
    return {'unstable':u,'stable':s,'basis':S,'symplectic_residual':symplectic,
            'diagonal_energy_residual':diagonal,'independent_Lorentz_residual':lorentz}

def logsinhc(x):
    if x==0:return 0.
    if abs(x)<1e-4:return -x*x/6+x**4/180
    if x>40:return math.log(2*x)-x-math.log1p(-math.exp(-2*x))
    return math.log(x/math.sinh(x))

def logadd(*values):
    largest=max(values)
    if largest==-math.inf:return largest
    return largest+math.log(sum(math.exp(x-largest) for x in values))

def upperexp(x):return math.exp(x) if x>=-744 else math.ulp(0.)

def spectrum(T,F,Bn,Bt,phi=4.3):
    if not (T>0 and F>0):raise ValueError('Local SN saddle needs positive F and thermal T')
    t=kb*T/qe;xs=math.sqrt(c/F);Us=phi-2*math.sqrt(c*F)
    omega0=math.sqrt(2*qe*c/(me*xs**3))
    bn,bt=qe*Bn/(me*omega0),qe*Bt/(me*omega0)
    m=modes(bn,bt);lam,stable=m['unstable'],m['stable']
    if lam<=0:return {'supported_local_unstable_channel':False,'physical_current':None,'reason':'No unstable saddle channel; not a zero gross thermal surface supply'}
    gamma=hbar*omega0*lam/(2*math.pi*qe)
    delta=hbar*omega0*stable/qe
    au=math.pi*gamma/t;ass=delta/(2*t)
    if au>=math.pi:return {'supported_local_unstable_channel':False,'physical_current':None,'reason':'Activated MB saddle integral is outside convergence domain'}
    quantum=1 if au==0 else au/math.sin(au)
    logJ=math.log(A0)+2*math.log(T)-Us/t+math.log(lam)+math.log(quantum)+logsinhc(ass)
    Jmb=math.exp(logJ)
    unstable_mean=t*(1-au*au/3 if au<1e-5 else au/math.tan(au))
    stable_mean=t*(1+ass*ass/3 if ass<1e-5 else ass/math.tanh(ass))
    mean=Us+unstable_mean+stable_mean
    # Full FD lies below MB, with a controlled split-energy bound even when
    # the naive second Boltzmann moment diverges. No total current is set zero.
    cut=Us*(1-gamma/t)
    if cut<=0:raise ValueError('No positive split for the declared activated-model bound')
    beta=1/gamma-1/t
    logPref=math.log(C)+math.log(lam)+math.log(gamma)+logsinhc(delta/(2*gamma))
    lowlog=logPref-math.log(beta)+(cut-Us)/gamma-cut/t
    errlog=logadd(lowlog,logJ-cut/t)
    neglog=logPref-Us/gamma-2*math.log(beta)
    abs_energy_log=logadd(logJ+math.log(abs(mean)),math.log(2)+neglog)
    lowenergy_log=logadd(logPref+(cut-Us)/gamma-cut/t+math.log(cut/beta+1/(beta*beta)),math.log(2)+neglog)
    heaterrlog=logadd(lowenergy_log,abs_energy_log-cut/t)
    return {'supported_local_unstable_channel':True,'T_K':T,'F_V_m':F,'B_normal_T':Bn,'B_tangent_T':Bt,
            'saddle_x_m':xs,'saddle_energy_eV':Us,'omega0_s-1':omega0,'omega_c_over_omega0':math.hypot(bn,bt),
            'unstable_frequency_over_omega0':lam,'stable_frequency_over_omega0':stable,
            'activated_MB_model_current_A_m2':Jmb,'log_model_current_A_m2':logJ,
            'activated_MB_model_mean_Fermi_energy_eV':mean,'activated_MB_complete_heat_W_m2':-Jmb*mean,
            'FD_current_difference_upper_A_m2':upperexp(errlog),'log_FD_current_difference_upper_A_m2':errlog,
            'FD_absolute_heat_difference_upper_W_m2':upperexp(heaterrlog),
            'FD_relative_current_bound':upperexp(errlog-logJ),
            'physical_current':None,'particle_phase_mapping_qualified':False,'full_SN_anharmonic_error_bound':None,
            'actual_magnetized_band_error_bound':None,'physical_foil_model_admitted':False}
