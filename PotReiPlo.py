#region Math/Import/Declarations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


alpha:float
a:list

#endregion

#region Parameters

x_list:list = np.linspace(0, 40, 261) # calculate at these x values

# These parameters can be changed in run-time using sliders

a0:complex = 0.4+0j # initial value ψ(0) 
a1:complex = 0+0j # second initial value (if needed for the recursion formula)
n:float = 4
E:float = 1
µ:float = 2
terms:int = 1000 # calculate this many coefficients

def declareFunctions():
    global calcCoeff, psi, alphaCalc 
    # calcCoeff : C x C x R x N --> C^terms, (a0, a1, alpha, terms) -> coefficients 
    # p : C --> C, z -> sum_k a_k z^k
    # psi : R --> C, the desired wave function in terms of p.
    # alphaCalc : R x N --> R, (E,n) -> alpha, how the parameter alpha for calcCoeff is computed from n and E
    # all of these can also depend on µ (use global µ)
    calcCoeff = calcCoeffB8A2
    psi = psiB8A2
    alphaCalc = alphaCalcB8

#endregion

#region B1A3
def calcCoeffB1A3(a0: complex, a1: complex, alpha: float, terms: int):
    a2 = -alpha * a0
    a3 = -(alpha-1) * a1 / 3
    a = [a0, a1, a2, a3]
    for i in range (2,terms+1):
        a.append( (i*(3-i)*a[i-2] + 2*(i*i- alpha)*a[i]) / ( (i+2) * (i+1)))
    return a

def psiB1A3(x):
    global µ # µ ist b (Breite)
    return p(np.tanh(x/µ))
    
def alphaCalcJanB1(E, n):
    return E + 1 - 1/2 * n**2

def alphaCalcJohannesB1(E, n):
    return E + 1
#endregion 

#region B7A2
def calcCoeffB7A2(a0:complex,a1:complex,alpha:float,terms: int):
    global n,µ
    a = [a0]
    for i in range (1,terms+1):
        a.append((µ+i-n) / (µ+i/2) * a[i-1]/i)
        if a[i] == 0:
            break
    return a

def psiB7A2(x):
    global µ,n
    return p(alpha/(2*n+1) * x) * x**µ * np.exp(-alpha/(2*n+1)*x)
#endregion

#region B8A2
def calcCoeffB8A2(a0: complex, a1: complex, alpha: float, terms: int):
    global µ # µ ist γ in der Aufgabenbearbeitung (x-Richtung-Wellenzahl) # alpha ist a0 in der Aufgabenbearbeitung ~ B (Magnetfeld) (E ist B hier)
    A1 = 2j/3 * np.sqrt(alpha*µ)
    A2 = np.sqrt(µ*µ-alpha/2)
    a = [a0, 0, 0, -A1/3 * a0, 0]
    for i in range (5,terms+1):
        a.append( (-4*(i-2)*a[i-2] + (6*i-15) * A1 * a[i-3] + 12 * A1 * A2 * a[i-5]) / i / (3-2*i) )
    return a

def psiB8A2(x):
    global µ,alpha
    A1 = 2j/3 * np.sqrt(alpha*µ)
    A2 = np.sqrt(µ*µ-alpha/2)

    return np.exp(- A2 * x + A1 * x**(3/2)) * p(np.sqrt(x))

def alphaCalcB8(E,n):
    return E
#endregion

#region Calculations

def absSquared(x:complex):
    return x.real*x.real + x.imag*x.imag

def polynom(a, x): # compute value of polynom with coefficients a at point x.
    res = 0
    for i in range(len(a) -1, -1, -1):
        res = res * x + a[i]
    return res

def p(z):
    global a
    return polynom(a, z)

output_values:list
output_names = ["Re ψ","Im ψ","|ψ|²","p(z)"]
output_colors = ["green","orange","blue","purple"]

def CalcYValues():
    global x_list, output_values
    psi_values = [psi(x) for x in x_list]
    output_values = [
        [psi.real for psi in psi_values],
        [psi.imag for psi in psi_values],
        [absSquared(psi) for psi in psi_values],
        [p(x/x_list[-1]) for x in x_list]
    ]

#endregion

#region Update from Slider

def updateKoeff():
    updateAlpha()
    global a, a0, a1, alpha, terms
    a = calcCoeff(a0, a1, alpha, terms)
    print(a[-1])

def updaten(val):
    global n
    n = val
    
    updateKoeff()
    redraw()


def updateE(val):
    global E
    E = val
    
    updateKoeff()
    redraw()

def updateµ(val):
    global µ
    µ = val
    redraw()

    
def updateA0(val):
    global a0
    a0 = val + 0*1j
    updateKoeff()
    redraw()
        
def updateReA1(val):
    global a1
    a1 = val + a1.imag
    updateKoeff()
    redraw()
    
def updateImA1(val):
    global a1
    a1 = a1.real + val * 1j
    updateKoeff()
    redraw()

def updateTerms(val):
    global terms
    terms = int(val)
    updateKoeff()
    redraw()

def redraw():
    global x_list, output_values, plots, fig, ax
    
    CalcYValues()

    ax.clear()
    plots = [ax.plot(x_list, output_values[i], label=output_names[i], color=output_colors[i], linewidth=1)[0]
        for i in range(len(output_values))]
        #plots[i].set_ydata(output_values[i])
    #plt.ylim(top=max(values[0],values[-1]))
    #plt.draw()
    ax.legend()
    ax.grid()
    fig.canvas.draw_idle()

def updateAlpha():
    global alpha,E,n
    alpha = alphaCalc(E,n)
#endregion

#region Plots

def configure_plot():
    global fig, ax
    #plt.figure(dpi=1200)
    
    fig, ax = plt.subplots()
    #plt.subplots_adjust(left=0.25, bottom=0.3)
    
    configure_sliders()
    configure_reset()

    
    redraw()
    plt.show()

def configure_sliders():
    global sliders
    #ax.margins(x=1)
    
    sliders = []
    number = 7
    lastposition= [0.25, 0.01 + (number + 1) * 0.04, 0.65, 0.03]


    def next_pos():
        lastposition[1] -= 0.04
        return lastposition

    color = 'lightgoldenrodyellow'

    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'n', 0, 20, valinit=n, valstep=1))
    sliders[-1].on_changed(updaten)
    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'E', -0.1, 1, valinit=E, valstep=0.001))
    sliders[-1].on_changed(updateE)
    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'µ', 0.1, 4, valinit=µ, valstep=0.1))
    sliders[-1].on_changed(updateµ)

    color = 'orange'

    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'a0', -0.7, 0.7, valinit=a0.real, valstep=0.01))
    sliders[-1].on_changed(updateA0)
    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'Re(a1)', -1, 1, valinit=a1.real, valstep=0.01))
    sliders[-1].on_changed(updateReA1)
    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'Im(a1)', -1, 1, valinit=a1.imag, valstep=0.01))
    sliders[-1].on_changed(updateImA1)

    color = "lightgreen"

    sliders.append(Slider(plt.axes(next_pos(), facecolor=color), 'terms', 4, 2000000, valinit=terms, valstep=1))
    sliders[-1].on_changed(updateTerms)

    
    plt.subplots_adjust(bottom= 0.02 + number* 0.06)

def configure_reset():
    global button

    resetax = plt.axes([0.05, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color='red', hovercolor='0.975')
    
    button.on_clicked(reset)

def reset(event):
    global sliders
    for slider in sliders:
        slider.reset()
#endregion

#region Main 

if __name__ != "__main__":
    print ("main called, but not as __main__")

declareFunctions()

updateKoeff()
print(a)
configure_plot()

#endregion