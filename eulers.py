from resources306 import *
from numpy import linspace, exp, sin, cos
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Slope field and exact solutions

def f(x,y):                   # Defines the right-hand side of the differential
    return -6*x*y             # equation y' = f(x,y)

def y_exact(x,C):             # Defines a 1-parameter family of solutions
    return C*exp(-3*x**2)     # parametrized by C

def get_C(x0,y0):             # Function to determine the correct C for a given
    return y0/exp(-3*x0**2)   # initial condition y(x0) = y0

def draw_slopefield_with_solutions():
    a = -2                                              # Left edge of plotting window (also lower edge)
    b =  2                                              # Right edge of plotting window (also upper edge)
    c = -5
    d = 5

    x = linspace(a,b,400)                               # Create an x-grid from a to b with 400 equally spaced points

    plt.figure(figsize=(8,8))                           # Initialize a figure for plotting
    slopefieldplot( f, a,b, c,d, .3, lw=1)              # Generate a slope field for y' = f(x,y)

    x0s = [0,0,-.5]                                     # These are the initial values of x
    y0s = [-3,1,2]                                      # These are the initial values of y.

    # The for-loop below plots particular solutions with initial condition y(x0) = y0,
    # where x0 and y0 come from the lists x0s and y0s respectively.
    for i,(x0,y0) in enumerate(zip(x0s,y0s)):
        C = get_C(x0,y0)                                # Determines the particular choice of C that satisfies the IC
        y = y_exact(x,C)                                # Gives the exact solution to the IVP
        plt.plot(x0,y0,'ko',ms = 10,zorder=2)
        plt.plot(x0,y0,'C{:}o'.format(i),markersize=8,zorder=3) # Plots the initial condition as a point
        plt.plot(x,y,'C{:}'.format(i),linewidth=4,label='Solution with $y({0}) = {1}$'.format(x0,y0),zorder=1)      # Plots the solution of the IVP
        
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    
    
def draw_linear_approximation():
    def y(x):                                     # Defines an example function y(x) = e^(x/2)
        return exp(.5*x)

    def La(x,a):                                  # Defines the linearization of y(x) = e^(x/2) at x = a
        return exp(.5*a) + .5*exp(.5*a)*(x - a)   # L_a(x) = e^(a/2) + 1/2 * e^(a/2) * (x - a)

    a = 2                                         # Set the a-value at which to take the linear approximation
    h = 1                                         # Set the step size for how far away we will approximate
    x = linspace(0,4,400)                         # Create an x-grid from 0 to 4 with 400 equally spaced points

    plt.figure(figsize=(8,8))                     # Initialize a figure for plotting

    plt.plot(x,y(x),'C0',linewidth=4,label='$y(x)$')             # Plot the function y(x)
    plt.plot(x,La(x,a),'k--',label='Linear approximation')             # Plot the linearization of y(x) at x=a
    plt.plot(a,y(a),'ko',ms=10)
    plt.plot(a,y(a),'C0o',markersize=8)          # Plot the coordinate (a, y(a))

    plt.plot(a+h,y(a+h),'ko',ms=10)
    plt.plot(a+h,y(a+h),'C1o',markersize=8)      # Plot the exact value y(a+h)
    plt.plot(a+h,La(a+h,a),'ko',ms=10)
    plt.plot(a+h,La(a+h,a),'C1o',markersize=8)   # Plot the approximation of y(a+h) using the linearization

    props = {'arrowstyle':'->'}
    plt.annotate('$(a,y(a))$',(a,y(a)),xytext=(1,3),fontsize=15,arrowprops=props)
    plt.annotate('$(a+h,y(a+h))$',(a+h,y(a+h)),xytext=(2,5),fontsize=15,arrowprops=props)
    plt.annotate('$(a+h,L_a(a+h))$',(a+h,La(a+h,a)),xytext=(2.75,3),fontsize=15,arrowprops=props)

    plt.xlim(0,4)
    plt.ylim(1,6)

    plt.xlabel('$x$')
    plt.ylabel('$y$');
    
    plt.legend()
    

def draw_Euler_animation():
    #######################################################################
    # Select initial condition
    #######################################################################

    x0s = [0,0,-.5]                # These are the initial values of x
    y0s = [-3,1,2]                 # These are the initial values of y.

    a = -2                                              # Left edge of plotting window
    b =  2                                              # Right edge of plotting window
    c = -5                                              # Lower edge of plotting window
    d =  5                                              # Upper edge of plotting window

    x = linspace(a,b,400)                               # Create an x-grid from a to b with 400 equally spaced points

    ys = []
    xlists = []
    ylists = []

    #######################################################################
    # Use Euler's method for a selected initial condition and step size
    #######################################################################

    h = .25                        # Choose the step size for Euler's method

    for (x0,y0) in zip(x0s, y0s):                                # Determines the particular choice of C that satisfies the IC
        ys.append(y_exact(x,get_C(x0,y0)))                                   # Gives the exact solution to the IVP
        
        xlist = [x0]                                     # Initialize a list of x-values at which we will  
                                                         # approximate the solution 
        ylist_Eul = [y0]                                 # Initialize a list of y-value approximations using Euler's method
        while xlist[-1] < b:                                 # Build up our list until we've stepped past the plotting window
            xx = xlist[-1]                                   # Grab the most recent x-value that was approximated at
            yy_Eul = ylist_Eul[-1]                           # Grab the most recent Euler y-value approximation
            yy_Eul = yy_Eul + f(xx,yy_Eul)*h                 # Use Euler's method to approximate the solution at the next x-value

            xx = xx + h                                      # Update the x-value to the next position
            xlist.append(xx)                                 # Add the new x-value to the list of x-values
            ylist_Eul.append(yy_Eul)                         # Add the new Euler approximation to the list of approximations
            
        xlists.append(xlist)
        ylists.append(ylist_Eul)
        
    #######################################################################
    # Generate a slopefield plot with exact solutions
    #######################################################################

    fig = plt.figure(figsize=(8,8))                     # Initialize a figure for plotting
    slopefieldplot( f, a,b, c,d, .3, lw=1)              # Generate a slope field for y' = f(x,y)

    alpha_low = .2

    y0_plots = []
    ex_plots = []
    for i in range(len(y0s)):
        y0_plot, = plt.plot(x0s[i],y0s[i],f'C{i}o',markersize=10,alpha=1)     # Plots the initial condition as a point
        ex_plot, = plt.plot(x,ys[i],f'C{i}',linewidth=4,label='Exact solution',zorder=1,alpha=1) # Plots the solution of the IVP
        y0_plots.append(y0_plot)
        ex_plots.append(ex_plot)

    # Plot the Euler approximation in red
    Eul_plot, = plt.plot([],[],'C3o--',linewidth=2,markersize=8,zorder=3,label='Euler approximation')
    lin_plot, = plt.plot([],[], 'k--', linewidth=2,alpha=.5,zorder=2,label='Linearization')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Euler\'s method animation')
    leg = plt.legend(loc='lower right')

    for y0_plot, ex_plot in zip(y0_plots, ex_plots):
        y0_plot.set_alpha(alpha_low)
        ex_plot.set_alpha(alpha_low)

    plt.grid()

    num_steps = 8
    num_frames_per_case = 3 * num_steps + 2
    num_frames = len(y0s) * num_frames_per_case

    def update(frame):
        i = frame // num_frames_per_case
        j3 = frame % num_frames_per_case
        j = j3//3
        
        xlist = xlists[i]
        ylist_Eul = ylists[i]
        y0_plot = y0_plots[i]
        ex_plot = ex_plots[i]
        
        if j3 == 0:
            for ii in range(len(ex_plots)):
                if ii == i:
                    ex_plots[ii].set_alpha(1)
                    y0_plots[ii].set_alpha(1)
                else:
                    ex_plots[ii].set_alpha(alpha_low)
                    y0_plots[ii].set_alpha(alpha_low)
        
        if j3%3 == 0:
            Eul_plot.set_xdata(xlist[:min(len(xlist),j+1)])
            Eul_plot.set_ydata(ylist_Eul[:min(len(xlist),j+1)])
        elif j3%3 == 1:
            lin_plot.set_xdata([])
            lin_plot.set_ydata([])
        else:
            mj = f(xlist[j], ylist_Eul[j])
            Lj = ylist_Eul[j] + mj*(x - xlist[j])
            lin_plot.set_xdata(x)
            lin_plot.set_ydata(Lj)
        if frame == num_frames-1:
            plt.close()
            
    anim = FuncAnimation(fig=fig, func=update, frames=num_frames, interval=1000, repeat=False)
    
    return anim