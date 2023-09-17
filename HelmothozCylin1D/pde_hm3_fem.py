"""
Basic function: (Helmothoz Equation)
    1/r *d/dr(r*du/dr) = 0
Functional: 
    1/2 * r * ((du/dr)^2 + u^2)
BC:
    Diriclet r = 1, u = 1, (set as right side colume)
    Natural r = 0, du/dr = 0, (don't need to do anything)
True solution:
    u(r) = bessel(0,r*i)/bessel(0,i)
"""
" =============================Finite Element Analysis ======================="
" =============================Import"
import numpy as np
from scipy import special
from matplotlib import pyplot
" =============================Basis parameters"
" The number of element"
def FEM(element_number):
    " The number of the node. Two nodes for each element"
    node_number = element_number + 1
    " For python convenience"
    N = node_number-1
    " The step size"
    h = 1/element_number
    " =============================Generate the grids"
    " Uniform grid is utilized here"
    grid = np.linspace(0,1,node_number)
    " =============================Finite element method"
    "%%%"
    " A is used to record the coefficient matrix"
    A = np.zeros([node_number,node_number])
    " Part B"
    A[0,0] = A[0,0] + (grid[1]**2-grid[0]**2)/(2*h**2)
    A[N,N] = A[N,N] + (grid[N]**2-grid[N-1]**2)/(2*h**2)
    for i in range(1,node_number-1):
        A[i,i] = A[i,i] + (grid[i]**2-grid[i-1]**2)/(2*h**2)
        A[i,i] = A[i,i] + (grid[i+1]**2-grid[i]**2)/(2*h**2)
    A[0,1] = A[0,1] - (grid[1]**2-grid[0]**2)/(2*h**2)
    A[N,N-1] = A[N,N-1] - (grid[N]**2-grid[N-1]**2)/(2*h**2)
    for i in range(1,node_number-1):
        A[i,i-1] = A[i,i-1] - (grid[i]**2-grid[i-1]**2)/(2*h**2)
        A[i,i+1] = A[i,i+1] - (grid[i+1]**2-grid[i]**2)/(2*h**2)
    " Part A"
    "fie_1 = a1*r + b1"
    "fie_2 = a2*r + b2"
    def intergral(a1,b1,a2,b2,r1,r2):
        return (a1*a2/4*(r2**4-r1**4)+(a1*b2+a2*b1)/3*(r2**3-r1**3)+b1*b2/2*(r2**2-r1**2))
    "i = 0"
    a1 = -1/h; b1 = h/h;
    a2 = -1/h; b2 = h/h;
    r1 = 0; r2 = h;
    A[0,0] = A[0,0] + intergral(a1,b1,a2,b2,r1,r2);
    a1 = -1/h; b1 = h/h;
    a2 = 1/h; b2 = -0/h;
    r1 = 0; r2 = h;
    A[0,1] = A[0,1] + intergral(a1,b1,a2,b2,r1,r2);
    "i = N"
    a1 = 1/h; b1 = -(N-1)*h/h;
    a2 = 1/h; b2 = -(N-1)*h/h;
    r1 = (N-1)*h; r2 = N*h;
    A[N,N] = A[N,N] + intergral(a1,b1,a2,b2,r1,r2);
    a1 = 1/h; b1 = -(N-1)*h/h;
    a2 = -1/h; b2 = (N)*h/h;
    r1 = (N-1)*h; r2 = N*h;
    A[N,N-1] = A[N,N-1] + intergral(a1,b1,a2,b2,r1,r2);
    for i in range(1,node_number-1):
        "A[i,i]"
        a1 = 1/h; b1 = -((i-1)*h)/h;
        a2 = 1/h; b2 = -((i-1)*h)/h;
        r1 = (i-1)*h; r2 = i*h;
        A[i,i] = A[i,i] + intergral(a1,b1,a2,b2,r1,r2);
        a1 = -1/h; b1 = ((i+1)*h)/h;
        a2 = -1/h; b2 = ((i+1)*h)/h;
        r1 = i*h; r2 = (i+1)*h;
        A[i,i] = A[i,i] + intergral(a1,b1,a2,b2,r1,r2);
        "A[i,i-1]"
        a1 = 1/h; b1 = -(i-1)*h/h;
        a2 = -1/h; b2 = i*h/h;
        r1 = (i-1)*h; r2 = i*h;
        A[i,i-1] = A[i,i-1] + intergral(a1,b1,a2,b2,r1,r2);
        "A[i,i+1]"
        a1 = -1/h; b1 = (i+1)*h/h;
        a2 = 1/h; b2 = -i*h/h;
        r1 = i*h; r2 = (i+1)*h;
        A[i,i+1] = A[i,i+1] + intergral(a1,b1,a2,b2,r1,r2);
    "%%%"
    " b vector"
    b = np.array([(np.zeros(node_number))])
    for i in range(0,node_number):
        b[0,i] = 0
    b = np.transpose(b)
    b[node_number-1,0] = 1
    for i in range(0,node_number):
        A[N,i] = 0
    A[N,N] = 1
    "%%%"
    "Solve A u = b"
    u_num = np.linalg.solve(A,b)
    " =============================Exact solution"
    """
    fig, ax = pyplot.subplots()
    r_vector = np.linspace(0,1,100)
    u_exact_vector = special.i0(r_vector) / special.i0(1)
    ax.plot(r_vector,u_exact_vector,'r',label = 'Exact Solution')
    ax.plot(grid,u_num,'-o',label = 'Numerical Solution')
    ax.legend()
    ax.set_xlabel('r')
    ax.set_ylabel('u')
    fig.savefig('pde.pdf')
    """
    " =============================Calculate the Error"
    u_exact = np.array([special.i0(grid) / special.i0(1)]);
    u_exact = np.transpose(u_exact)
    ur_exact = special.i1(grid) / special.i0(1);
    " L2 norm"
    u_error_L2 = np.sqrt(np.sum((u_num - u_exact)**2*h))
    " Calculate ur"
    ur_num = np.zeros(node_number)
    ur_num[0] = (u_num[1]-u_num[0])/h
    ur_num[N] = (u_num[N]-u_num[N-1])/h
    for i in range(1,N):
        ur_num[i] = (u_num[i+1]-u_num[i-1])/2/h
    ur_error_L2 = np.sqrt(np.sum((ur_num - ur_exact)**2*h))
    " H1 norm"
    u_error_H1 = np.sqrt(u_error_L2**2+ur_error_L2**2)
    return [u_error_L2,u_error_H1]
element_vector = []
L2 = []
H1 = []
for i in range(0,25):
    element_number = i * 4 + 4
    element_vector.append(element_number)
    error = FEM(element_number)
    L2.append(error[0])
    H1.append(error[1])
fig, ax = pyplot.subplots()
ax.plot(np.log(element_vector),np.log(L2),'-',label = 'L2 norm Error')
ax.plot(np.log(element_vector),np.log(H1),'--',label = 'H1 norm Error')
ax.legend()
ax.set_xlabel('ln(N)')
ax.set_ylabel('ln(Error)')
fig.savefig('pde.pdf')