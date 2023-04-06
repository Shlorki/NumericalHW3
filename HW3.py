import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg

###### Define Backward Time Centered Space method
def BTCS(u0,a,b,T,Nx,Nt):
    # Setup discretization
    [t,ht] = np.linspace(0,T,Nt,retstep=True); [x,hx] = np.linspace(a,b,Nx,retstep=True)

    # Initialize solution matrix
    U = np.zeros((Nx,Nt)); U[Nx-1,0] = u0(x[0])
    for k in range(0,Nx-1):
        U[k,0] = u0(x[k])

    # Construct Backward Time, Centered Space difference equation matrix
    L = ht/hx
    Ld = L/2*np.ones(Nx-2); D = np.ones(Nx-1)
    A = sp.diags([Ld,D,-Ld],[-1,0,1]); A = A.tolil()
    A[0,Nx-2] = L/2; A[Nx-2,0] = -L/2; A = A.tocsr()

    # iterate time steps
    for n in range(1,Nt):
        b = U[range(Nx-1),n-1]
        U[range(Nx-1),n] = sp.linalg.spsolve(A,b); U[Nx-1,n] = U[0,n]

    return U,x,t

###### Test BTCS and plot numerical solution
def u0(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

# Run the Backward Time Centered Space method
Nx = 201; Nt = 225
U,x,t = BTCS(u0,0,1,1,Nx,Nt)
x1 = x+1

# Setup grid lines
N_int = 3000
X = np.linspace(0,1,3000)
grid1 = np.ones(N_int); grid2 = -0.1*np.ones(N_int); grid3 = np.zeros(N_int); grid4 = np.linspace(0,1,N_int)
grid5 = 0.5*np.ones(N_int); grid6 = 1.5*np.ones(N_int)

# Plot grid lines and numerical solution
plt.figure(1)
plt.plot(X,grid1,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid1,color='k',linestyle='--',linewidth=.75)
plt.plot(X,grid3,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid3,color='k',linestyle='--',linewidth=.75)
plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(grid5,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid6,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(x,U[:,0],color='b', label='t=0'); plt.plot(x+1,U[:,0],color='b')
plt.plot(x,U[:,int((Nt-1)/2)],color='g', label='t=0.5')
plt.plot(x+1,U[:,int((Nt-1)/2)],color='g')
plt.plot(x,U[:,Nt-1],color='r', label='t=1')
plt.plot(x+1,U[:,Nt-1],color='r')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Backward Time Centered Space')
plt.legend(loc='center left', fontsize = 'small')
plt.savefig('BTCS.png')
plt.show()


###### Define the Crank-Nicholson Method
def CNCS(u0,a,b,T,Nx,Nt):
    # Setup discretization
    [t,ht] = np.linspace(0,T,Nt,retstep=True); [x,hx] = np.linspace(a,b,Nx,retstep=True)

    # Initialize solution matrix
    U = np.zeros((Nx,Nt)); U[Nx-1,0] = u0(x[0])
    for k in range(Nx-1):
        U[k,0] = u0(x[k])

    # Construct Crank-Nicholson Centered Space difference equation matrix
    L = ht/hx
    Ld = L/4*np.ones(Nx-2); D = np.ones(Nx-1)
    A = sp.diags([Ld,D,-Ld],[-1,0,1]); A = A.tolil()
    A[0,Nx-2] = L/4; A[Nx-2,0] = -L/4; A = A.tocsr()

    # iterate time steps
    for n in range(1,Nt):
        b = np.zeros(Nx-1)
        b[0] = U[0,n-1]+L/4*(U[1,n-1]-U[-2,n-1]); b[Nx-2] = U[Nx-2,n-1]+L/4*(U[0,n-1]-U[Nx-3,n-1])
        b[range(1,Nx-2)] = U[range(1,Nx-2),n-1] + L/4*(U[range(2,Nx-1),n-1]-U[range(Nx-3),n-1])
        U[range(Nx-1),n] = sp.linalg.spsolve(A,b); U[Nx-1,n] = U[0,n]

    return U,x,t

###### Test CNCS and plot numerical solution
def u0(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

# Run the Crank-Nicholson method with central difference in space
Nx = 201; Nt = 225
U,x,t = CNCS(u0,0,1,1,Nx,Nt)
x1 = x+1

# Setup grid lines
N_int = 3000
X = np.linspace(0,1,3000)
grid1 = np.ones(N_int); grid2 = -0.1*np.ones(N_int); grid3 = np.zeros(N_int); grid4 = np.linspace(0,1,N_int)
grid5 = 0.5*np.ones(N_int); grid6 = 1.5*np.ones(N_int)

# Plot grid lines and numerical solution
plt.figure(1)
plt.plot(X,grid1,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid1,color='k',linestyle='--',linewidth=.75)
plt.plot(X,grid3,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid3,color='k',linestyle='--',linewidth=.75)
plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(grid5,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid6,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(x,U[:,0],color='b', label = 't=0'); plt.plot(x+1,U[:,0],color='b')
plt.plot(x,U[:,int((Nt-1)/2)],color='g', label = 't=0.5')
plt.plot(x+1,U[:,int((Nt-1)/2)],color='g')
plt.plot(x,U[:,Nt-1],color='r', label = 't=1')
plt.plot(x+1,U[:,Nt-1],color='r')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Crank-Nicholson Time Centered Space')
plt.legend(loc='center left', fontsize = 'small')
plt.savefig('CNCS.png')
plt.show()


###### Define the Lax-Friedrichs method
def LaxFried(u0,a,b,T,Nx,Nt):
    # Setup discretization
    [t,ht] = np.linspace(0,T,Nt,retstep=True); [x,hx] = np.linspace(a,b,Nx,retstep=True)

    # Initialize solution matrix
    U = np.zeros((Nx,Nt)); U[Nx-1,0] = u0(x[0])
    for k in range(Nx-1):
        U[k,0] = u0(x[k])

    # Iterate the solution through time using Lax-Friedrichs
    L = ht/hx
    for n in range(1,Nt):
        U[0,n] = (U[1,n-1]+U[-1,n-1])/2 + L/2*(U[1,n-1]-U[-1,n-1]); U[-1,n] = U[0,n]
        U[Nx-2,n] = (U[0,n-1]+U[Nx-3,n-1])/2 + L/2*(U[0,n-1]-U[Nx-3,n-1])
        U[range(1,Nx-2),n] = (U[range(2,Nx-1),n-1]+U[range(Nx-3),n-1])/2 + L/2*(U[range(2,Nx-1),n-1]-U[range(Nx-3),n-1])

    return U,x,t

###### Test LaxFried and plot numerical solution
def u0(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

# Run the Lax-Friedrichs method
Nx = 201; Nt = 225
U,x,t = LaxFried(u0,0,1,1,Nx,Nt)
x1 = x+1

# Setup grid lines
N_int = 3000
X = np.linspace(0,1,3000)
grid1 = np.ones(N_int); grid2 = -0.1*np.ones(N_int); grid3 = np.zeros(N_int); grid4 = np.linspace(0,1,N_int)
grid5 = 0.5*np.ones(N_int); grid6 = 1.5*np.ones(N_int)

# Plot grid lines and numerical solution
plt.figure(1)
plt.plot(X,grid1,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid1,color='k',linestyle='--',linewidth=.75)
plt.plot(X,grid3,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid3,color='k',linestyle='--',linewidth=.75)
plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(grid5,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid6,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(x,U[:,0],color='b', label = 't=0'); plt.plot(x+1,U[:,0],color='b')
plt.plot(x,U[:,int((Nt-1)/2)],color='g', label = 't=0.5')
plt.plot(x+1,U[:,int((Nt-1)/2)],color='g')
plt.plot(x,U[:,Nt-1],color='r', label = 't=1')
plt.plot(x+1,U[:,Nt-1],color='r')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Lax-Friedrichs')
plt.legend(loc='center left', fontsize = 'small')
plt.savefig('LaxFried.png')
plt.show()


###### Define the Lax-Wendroff method
def LaxWend(u0,a,b,T,Nx,Nt):
    # Setup discretization
    [t,ht] = np.linspace(0,T,Nt,retstep=True); [x,hx] = np.linspace(a,b,Nx,retstep=True)

    # Initialize solution matrix
    U = np.zeros((Nx,Nt)); U[Nx-1,0] = u0(x[0])
    for k in range(Nx-1):
        U[k,0] = u0(x[k])

    # Iterate the solution through time using Lax-Wendroff
    L = ht/hx
    for n in range(1,Nt):
        U[0,n] = U[0,n-1] + L/2*(U[1,n-1]-U[-1,n-1]) + L**2/2*(U[1,n-1]-2*U[0,n-1]+U[-1,n-1]); U[-1,n] = U[0,n]
        U[Nx-2,n] = U[Nx-2,n-1] + L/2*(U[0,n-1]-U[Nx-3,n-1]) + L**2/2*(U[0,n-1]-2*U[Nx-2,n-1]+U[Nx-3,n-1])
        U[range(1,Nx-2),n] = U[range(1,Nx-2),n-1] + L/2*(U[range(2,Nx-1),n-1]-U[range(Nx-3),n-1]) + L**2/2*(U[range(2,Nx-1),n-1]-2*U[range(1,Nx-2),n-1]+U[range(Nx-3),n-1])

    return U,x,t

###### Test LaxWend and plot the numerical solution
def u0(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

# Run the Lax-Wendroff method
Nx = 201; Nt = 225
U,x,t = LaxWend(u0,0,1,1,Nx,Nt)
x1 = x+1

# Setup grid lines
N_int = 3000
X = np.linspace(0,1,3000)
grid1 = np.ones(N_int); grid2 = -0.1*np.ones(N_int); grid3 = np.zeros(N_int); grid4 = np.linspace(0,1,N_int)
grid5 = 0.5*np.ones(N_int); grid6 = 1.5*np.ones(N_int)

# Plot grid lines and numerical solution
plt.figure(1)
plt.plot(X,grid1,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid1,color='k',linestyle='--',linewidth=.75)
plt.plot(X,grid3,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid3,color='k',linestyle='--',linewidth=.75)
plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(grid5,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid6,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(x,U[:,0],color='b', label = 't=0'); plt.plot(x+1,U[:,0],color='b')
plt.plot(x,U[:,int((Nt-1)/2)],color='g', label = 't=0.5')
plt.plot(x+1,U[:,int((Nt-1)/2)],color='g')
plt.plot(x,U[:,Nt-1],color='r', label = 't=1')
plt.plot(x+1,U[:,Nt-1],color='r')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Lax-Wendroff')
plt.legend(loc='center left', fontsize = 'small')
plt.savefig('LaxWend.png')
plt.show()


###### Define the 4th order Runge-Kutta time and 4th Order Compact Differences in Space
def RK4CD4(u0,a,b,T,Nx,Nt):
    # Setup discretization
    [t,ht] = np.linspace(0,T,Nt,retstep=True); [x,hx] = np.linspace(a,b,Nx,retstep=True)

    # Initialize solution matrix
    U = np.zeros((Nx,Nt)); U[-1,0] = u0(x[0])
    for k in range(Nx-1):
        U[k,0] = u0(x[k])

    # Fourth-order compact difference scheme for spatial derivative
    def CD4(u):
        L = 1/4*np.ones(Nx-2); D = np.ones(Nx-1)
        A = sp.diags([L,D,L],[-1,0,1]); A = A.tolil()
        A[0,Nx-2] = 1/4; A[Nx-2,0] = 1/4; A = A.tocsr()

        b = np.zeros(Nx-1)
        b[0] = 3/(4*hx)*(u[1]-u[-1]); b[Nx-2] = 3/(4*hx)*(u[0]-u[Nx-3])
        b[range(1,Nx-2)] = 3/(4*hx)*(u[range(2,Nx-1)]-u[range(Nx-3)])

        du = sp.linalg.spsolve(A,b)
        return du

    # RK4 for time stepping
    for n in range(1,Nt):
        K1 = CD4(U[range(Nx-1),n-1]); K2 = CD4(U[range(Nx-1),n-1]+ht/2*K1)
        K3 = CD4(U[range(Nx-1),n-1]+ht/2*K2); K4 = CD4(U[range(Nx-1),n-1]+ht*K3)

        U[range(Nx-1),n] = U[range(Nx-1),n-1] + ht*(K1+2*K2+2*K3+K4)/6; U[-1,n] = U[0,n]

    return U,x,t

###### Test RK4CD4 and plot the numerical solution
def u0(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

# Run the Backward Time Centered Space method
Nx = 201; Nt = 225
U,x,t = RK4CD4(u0,0,1,1,Nx,Nt)
x1 = x+1

# Setup grid lines
N_int = 3000
X = np.linspace(0,1,3000)
grid1 = np.ones(N_int); grid2 = -0.1*np.ones(N_int); grid3 = np.zeros(N_int); grid4 = np.linspace(0,1,N_int)
grid5 = 0.5*np.ones(N_int); grid6 = 1.5*np.ones(N_int)

# Plot grid lines and numerical solution
plt.figure(1)
plt.plot(X,grid1,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid1,color='k',linestyle='--',linewidth=.75)
plt.plot(X,grid3,color='k',linestyle='--',linewidth=.75); plt.plot(X+1,grid3,color='k',linestyle='--',linewidth=.75)
plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid1,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(grid5,grid4,color='k',linestyle='--',linewidth=.75); plt.plot(grid6,grid4,color='k',linestyle='--',linewidth=.75)
plt.plot(x,U[:,0],color='b', label = 't=0'); plt.plot(x+1,U[:,0],color='b')
plt.plot(x,U[:,int((Nt-1)/2)],color='g', label = 't=0.5')
plt.plot(x+1,U[:,int((Nt-1)/2)],color='g')
plt.plot(x,U[:,Nt-1],color='r', label = 't=1')
plt.plot(x+1,U[:,Nt-1],color='r')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('RK4 Time CD4 Space')
plt.legend(loc='center left', fontsize = 'small')
plt.savefig('RK4CD4.png')
plt.show()