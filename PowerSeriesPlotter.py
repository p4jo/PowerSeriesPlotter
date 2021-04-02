
#region imports
from numba import jit, njit, typeof, typed, types
from numba.typed import Dict, List
import numpy as np
import matplotlib.figure as mplfig
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
# import timeit
import time
#endregion

#region definitions
class sliderConf(object):
    def __init__(self, name:str, init:float = 0, min:float = 0, max:float = 10, step:float = 1, redrawImmediately:bool = True, color:str = 'lightgoldenrodyellow'):
        self.name = name
        self.init = init
        self.step = step
        self.min = min
        self.max = max
        self.redrawImmediately = redrawImmediately
        self.color = color
        


@njit
def calcCoeff(a0, a1, terms, par):
    n = int(terms.real)
    if n < 2:
        n = 2
    res = np.empty(n,dtype=np.complex128)
    res[0]=a0
    res[1]=a1
    return res

@njit
def psi(x, a, params):
    return polynom(a,x)

@njit   
def updateExtraParams(par):
    pass

windowTitle = "Power series plotter"
x_list = np.linspace(0, 40, 261)
a = np.empty(1,dtype=np.complex)
fig:mplfig.Figure
# ax:plt.axes.Axes
plots:list #unneccesary
resetButton:Button
redrawButton:Button
sliders:list
parameters:Dict = Dict.empty(types.string,types.complex128)
# calculatedParameters:Dict = Dict.empty(str,complex)

#endregion

# These parameters can be changed in run-time using sliders
sliderSetup = [
    # sliderConf(name:str, init:float = 0, min:float = 0, max:float = 10, step:float = 1, redrawImmediately:bool = True, color:str = 'lightgoldenrodyellow')

    sliderConf('a0',0.1, -0.7, 0.7, 0.01, color = 'orange'), # initial value ψ(0)
    sliderConf('Re(a1)',0,-1, 1,0.01, color = 'orange'), # second initial value (if needed for the recursion formula)
    sliderConf('Im(a1)',0,-1, 1,0.01, color = 'orange'),
    sliderConf('terms',1000,1,1000000,1, False, color = "lightgreen"),
] 

def setup(calculateCoefficients,ψ, updateParameters = None,paramSliderSetup = None):
    """
    Opens a interactive window for plotting a power series where only the recursion formula for the coefficients is known. 


     calcCoeff: C x C x N x params --> C^terms, (a0, a1, alpha, terms, parameters) -> coefficients
     psi: R x C^terms x params --> C, (x, a, parameters) -> ψ(x), the desired wave function in terms of polynom(a, x), the power series truncated at terms terms.
        Example: ψ(x) = polynom(a,np.tanh(x))
     calcCoeff and psi
     updateExtraParams: a function of (parameters) that updates (some) calculatedParameters
     paramSliderSetup: A list of sliderConf objects
     x_list: calculate ψ at these x values. Example: numpy.linspace(0,10,200)
    """
    global psi,calcCoeff,sliderSetup,updateExtraParams
    if not calculateCoefficients is None:
        calcCoeff = njit(calculateCoefficients)
    if not ψ is None:
        # @njit
        # def PSI(x, a, params):
        #     return ψ(x)
        psi = njit(ψ)
    if not paramSliderSetup is None:
        sliderSetup.extend(paramSliderSetup)
    if not updateParameters is None:
        updateExtraParams = jit(updateParameters)

        

#endregion

#region Calculations
@njit
def absSquared(x:complex):
    return x.real*x.real + x.imag*x.imag

@njit
def polynom(a, x): 
    """compute value of polynom with coefficients a at point x. # Can be replaced by np.polyval"""
    res = 0
    for i in range(len(a) -1, -1, -1):
        res = res * x + a[i]
    return res

output_values:List
output_names = ["Re ψ","Im ψ","|ψ|²","Re p(z)"]
output_colors = ["green","orange","blue","purple"]

@njit
def CalcYValues(a, parameters):
    psi_values = [psi(x, a, parameters) for x in x_list]
    

    return np.array([
        [psi.real for psi in psi_values],
        [psi.imag for psi in psi_values],
        [absSquared(psi) for psi in psi_values],
        [polynom(a,x/x_list[-1]).real for x in x_list]
    ])


#endregion

#region Update from Slider

def calculateExtraParameters():
    global a, parameters
    updateExtraParams(parameters)
    print("Calculating coefficients", end = " ... ", flush=True)

    a = calcCoeff(complex(parameters["a0"]), parameters["Re(a1)"] + 1j * parameters["Im(a1)"], int(parameters["terms"].real), parameters)

    print(" Done. The last coefficient (a_{}) is {}.".format(len(a)-1,a[-1]))

def updateParams(newValue, key: str, redrawImmediately: bool = True):
    global parameters
    parameters[key] = newValue
    calculateExtraParameters()
    if redrawImmediately:
        redraw()

def createUpdateFunction(na,rI):
    return lambda v: updateParams(v,na,rI)

def redraw():
    global output_values, plots, fig, ax, parameters
    # 
    print("Calculating ψ values", end = " ... ", flush=True)
    before = time.time()
    output_values = CalcYValues(a, parameters)
    print("Took "+str(time.time()-before)+"seconds")
    print("Redrawing image", end = " ... ", flush=True)
    ax.clear()
    plots = [
        ax.plot(x_list, output_values[i], label=output_names[i], color=output_colors[i], linewidth=1)[0]
        for i in range(len(output_values))]
        #plots[i].set_ydata(output_values[i])
    #plt.ylim(top=max(values[0],values[-1]))
    #plt.draw()
    ax.legend()
    ax.grid()
    fig.canvas.draw_idle()
    print("Done.")

#endregion

#region Plots

def configure_plot():
    global fig, ax
    #plt.figure(dpi=1200)
    # plt.xkcd()

    fig, ax = plt.subplots(num=windowTitle)
    #plt.subplots_adjust(left=0.25, bottom=0.3)
    
    configure_sliders()
    configure_buttons()

    
    redraw()
    plt.show()

def configure_sliders():
    global sliders, sliderSetup,createUpdateFunction
    #ax.margins(x=1)
    
    sliders = []
    number = len(sliderSetup)
    lastposition = [0.25, 0.01 + (number + 1) * 0.04, 0.65, 0.03]

    def next_pos():
        lastposition[1] -= 0.04
        return lastposition
    
    for s in sliderSetup:
        sliders.append(Slider(plt.axes(next_pos(), facecolor=s.color), s.name, s.min, s.max, valinit=s.init, valstep=s.step))
        sliders[-1].on_changed(createUpdateFunction(na = s.name, rI = s.redrawImmediately))
    

    plt.subplots_adjust(bottom= 0.02 + number * 0.06)

def configure_buttons():
    global redrawButton, resetButton

    resetax = plt.axes([0.05, 0.025, 0.1, 0.04])
    resetButton = Button(resetax, 'Reset', color='red', hovercolor='0.5')
    resetButton.on_clicked(reset)
    
    redrawax = plt.axes([0.05, 0.08, 0.1, 0.04])
    redrawButton = Button(redrawax, 'redraw', color='green', hovercolor='0.5')
    redrawButton.on_clicked(lambda e: redraw())

def reset(event):
    global sliders
    for slider in sliders:
        slider.reset()

#endregion

#region Main 
def start():
    for s in sliderSetup:
        parameters[s.name] = s.init

    calculateExtraParameters()

    # print(a)
    configure_plot()

#endregion