from numpy import *
from scipy.special import factorial
import matplotlib.pyplot as plt

def get_coefficients(n):
    # This function computes the coefficients a_0,a_1,...,a_n
    # of a power series
    return array([a(k) for k in range(n+1)])

def get_power_series(n):
    # This function constructs a series approximation of degree n
    coefficients = get_coefficients(n)
    return lambda x: array([ak*x**k for k,ak in enumerate(coefficients)]).sum(axis=0)

def plot_series_approximations(x,y_exact=None,k0=0,L=-5,U=5):
    
    # This function plots 9 different series approximations
    
    # x is an array of x-values at which the approximations will be evaluated.
    
    # y_exact is an exact function that is being approximated (optional).
    # If provided, the exact solution will be plotted as a dashed black line.
    
    # k0 is the starting approximation degree (optional). 
    # The function will plot the approximations of degree k0, k0+1, k0+2, ... k0 + 9.
    
    # L and U are the lower and upper y bounds for the plots (optional).
    
    fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(10,10))

    for k,ax in enumerate(axes.flatten()):
        power_series = get_power_series(k+k0)
        y = power_series(x)
        ax.plot(x,y,'C{:}'.format(k),linewidth=3)
        if y_exact is not None: ax.plot(x,y_exact,'k--',linewidth=2)
        ax.grid()
        ax.set_xlim((x[0],x[-1]))
        ax.set_ylim((L,U))
        ax.set_title('Degree-{:} approximation'.format(k+k0))

    for ax in axes[:,1:].flatten(): ax.set_yticklabels('')
    for ax in axes[:-1,:].flatten(): ax.set_xticklabels('')
    for ax in axes[:,0]: ax.set_ylabel('$y$')
    for ax in axes[-1,:]: ax.set_xlabel('$x$')

    fig.tight_layout()