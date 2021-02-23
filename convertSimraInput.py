import pandas as pd
# Typical parm.dat file
#1, 0, 250, 0        ! isolve, kont, iout, istrat
#8000, 0.125         ! nstep, dt
#100, 0.01, 0.25     ! niter, eps, supg  
#0.000015, 1.3       ! visc, rho0
#2000., 20.          ! hscal, uinf   
#0, 1                ! itract, irelax
#0.0                 ! mesh rotation
#0.25, 1.0           ! cr_min, cr_max

# Define default parameters
# Missing parameters: isolve = 1, supg = 0.25, itract = 0
parms = {}
parms['nstep'] = 4000    # NSTEP
parms['dt'] = 1.0        # DT
parms['istrat'] = 1      # STRATIFICATION
parms['iout'] = 50       # SAVE_RES
parms['kont'] = 0        # RESTART
parms['niter'] = 200     # MAXIT
parms['eps'] = 0.01      # CC
parms['visc'] = 1.48e-5  # VIS
parms['rho0'] = 1.225    # RHOREF
parms['uinf'] = 10.      # UREF
parms['hscal'] = 2000    # LENREF
parms['tau_su'] = 0.25   # TAU_SU
parms['theta'] = 1.0     # PSI
parms['itemp'] = 0       # ROTBOUND
parms['cr_min'] = 0.5    # COURANT_MIN
parms['cr_max'] = 3.0    # COURANT_MAX
parms['irelax'] = 1.0    # RELAX


for line in open("parm.dat"):
    fields=line.replace('!',',').replace('\n','').replace(' ','').split(',')
    n = len(fields)
    for i in range(0,n//2):
        parms[fields[n//2 + i]] = float(fields[i])
    

print(parms)
        
with open('simra.in', 'w') as f:
    f.write('&PARAM_DATA\n')
    f.write(' NSTEP=%d,        ! Number of time steps\n' % parms['nstep'])
    f.write(' DT=%g,          ! Initial time step\n' % parms['dt'])
    f.write(' STRATIFICATION=%d,  ! = 1: stratification included\n' % parms['istrat'])
    f.write(' SAVE_RES=%d,      ! time-steps between output\n' % parms['iout'])
    f.write(' RESTART=%d,         ! Toggle startup from previous run \n' % parms['kont'])
    f.write(' MAXIT=%d,         ! max number of iterations\n' % parms['niter'])
    f.write(' CC=%g,           ! Tolerance level for linear solver\n' % parms['eps'])
    f.write(' VIS=%g,       ! kinematic viscosity\n' % parms['visc'])
    f.write(' RHOREF=%g,        ! reference density\n' % parms['rho0'])
    f.write(' UREF=%g,           ! reference velocity\n' % parms['uinf'])
    f.write(' LENREF=%g,       ! reference length\n' % parms['hscal'])
    f.write(' TAU_SU=%g,       ! SUPG parameter (0, 10, 0.25)\n' % parms['tau_su'])
    f.write(' PSI=%g,             ! Time stepping scheme (Theta method), PSI=1 (backward Euler), PSI=0.5 (Crank-Nicolson)\n' % parms['theta'])
    f.write(' ROTBOUND=0,        ! variable outer side boundary BCs\n')
    f.write(' COURANT_MIN=%g,  ! min acceptable Cr-number, for time-step control\n' % parms['cr_min'])
    f.write(' COURANT_MAX=%g,     ! max acceptable Cr-number, for time-step control\n' % parms['cr_max'])
    f.write(' RELAX=%g/           ! outer BC relaxation\n' % parms['irelax'])
    f.write('\n')
    f.write('&EXTRA_DATA\n')
    f.write(' TURB=1,            ! = 1: turbulence model included\n')
    f.write(' SAVE_HIST=20,      ! time-steps between saved output\n')
    f.write(' WRITE_NORMS=.TRUE.,! print absolute norm of updates\n')
    f.write(' SOURCE=0/          ! Source function flag\n')

