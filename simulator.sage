# reference: Improved Provable Reduction of NTRU and Hypercubic Lattices by Henry Bambury and Phong Q. Nguyen
# author of the code: Henry Bambury
# year: 2024
# tested with sage 9.5 and python 3.11.4


from sage.all import RealBallField, n, exp, pi

# Precision. 1000 is overkill
RBF = RealBallField(100)


def arb_inc_beta(x,a,b):
    """
    Arbitrary precision for the regularised incomplete beta function with parameters a,b
    """
    return RBF(a).beta(b,RBF(x))/RBF(a).beta(b)

def dichotomy_arb(a,b,dim,proj_dim,Nb_targets=1,eps=1e-10):
    """
    Uses dichotomy to find the solution to (1-inc_beta(x,dim/2,(dim-proj_dimension)/2)**Nb_targets = 0.5 with precision eps, between a and b. This is a good approximation of the integral of (1-inc_beta(x,dim/2,(dim-proj_dimension)/2)**Nb_targets.
    """
    delta = 1
    while delta > eps:
        m = (a + b) / 2
        delta = abs(b - a)
        if abs((1-arb_inc_beta(m,proj_dim/2,(dim-proj_dim)/2))**Nb_targets - RBF(0.5)) < RBF(eps):
            return m
        elif ((1-arb_inc_beta(a,proj_dim/2,(dim-proj_dim)/2))**Nb_targets - RBF(0.5)) * ((1-arb_inc_beta(m,proj_dim/2,(dim-proj_dim)/2))**Nb_targets - RBF(0.5))  > RBF(0):
            a = m
        else:
            b = m
    return m

def exp_min_sq_proj(dim,proj_dim,Nb_targets=1,eps=1e-10):
    """
    Computes the expectation for the minimal squared norm among Nb_targets projections of independant unit vectors in dimension dim onto dimension proj_dim
    """
    if(dim==proj_dim):
        return 1
    return dichotomy_arb(0, 1, dim,proj_dim,Nb_targets=Nb_targets,eps=eps)    
    
def f_primal(dim,beta,sq_norm_short_vector,vol,Nb_targets=1,eps=1e-10):
    """
    Derives the primal attack equation: the expected norm of the shortest projection of a shortest vector versus the GSA estimates for BKZ-reduced bases.
    """
    beta = RBF(beta)
    vol = RBF(vol)
    sq_min_proj = exp_min_sq_proj(dim,beta,Nb_targets=Nb_targets,eps=eps)
    return ((RBF(sq_min_proj)*RBF(sq_norm_short_vector))/(((beta/(RBF(2)*RBF(exp(1))*RBF(pi)))*(beta*RBF(pi))**(1/beta))**((RBF(2)*beta-RBF(dim)-1)/(beta-1)))/(vol**(RBF(2)/RBF(dim)))-RBF(1))    

def dichotomy_primal(a,b,dim,sq_norm_short_vector,vol,Nb_targets=1,eps=1e-10):
    """
    Finds a solution to the primal attack equation between a and b. The approximation used for BKZ vectors with small blocksize (<50) is very bad and should be replaced by hardcoded values, so don't expect anything if the answer is <50. 
    """
    delta = 1
    while delta > eps:
        m = (a + b) / 2
        delta = abs(b - a)
        if abs(f_primal(dim,m,sq_norm_short_vector,vol,Nb_targets=Nb_targets,eps=eps)) < RBF(eps):
            return m
        elif f_primal(dim,a,sq_norm_short_vector,vol,Nb_targets=Nb_targets,eps=eps)*f_primal(dim,m,sq_norm_short_vector,vol,Nb_targets=Nb_targets,eps=eps)  > RBF(0):
            a = m
        else:
            b = m
    return m

def gen_hypercubic():
    """
    Generates intersection beta for the hypercubic lattice
    """
    n_list = list(range(200,1000,25))
    for n in n_list:
        print("Dimension",n,"| Optimal beta with 1 vector:",dichotomy_primal(0.2*n,0.8*n,n,1,1,Nb_targets=1,eps=1e-10),"| Optimal beta with "+str(n)+" vectors:",dichotomy_primal(0.2*n,0.8*n,n,1,1,Nb_targets=n,eps=1e-10))

def main():
    print("This code helps simulate the impact of the number of short vectors in a lattice, following the primal attack heuristic.")
    s = input("Select mode of operation: hypercubic_auto (a) or manual (m)"+"\n")
    if(s == 'a'):
        gen_hypercubic()
    elif(s == 'm'):
        n = int(input("Dimension of L:"))
        vol = float(input("Volume of L:"))
        lambda1 = float(input("Norm of the short vector:"))
        N = int(input("Number of targets:"))
        try:
            assert(n>=200)
        except:
            print("Small dimensions not implemented, estimates will fail.")
            return 0
        s = ''
        if(N>1):
            s = 's'
        print("Dimension",n,"| Optimal beta with",N,"vector"+s+":",dichotomy_primal(0.2*n,0.8*n,n,lambda1**2,vol,Nb_targets=N,eps=1e-10)) 
    else:
        print("Invalid option")
        
# usage: sage simulator.sage
if __name__ == "__main__": 
    main() 
