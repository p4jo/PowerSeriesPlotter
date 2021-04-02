import PowerSeriesPlotter as psp
from PowerSeriesPlotter import polynom
import numpy as np
import warnings
warnings.filterwarnings(action="error", category=np.ComplexWarning)

paramSliderSetup = [
    psp.sliderConf('E',0.1,-0.1, 1,0.001),
    psp.sliderConf('µ',2,0.1,4,0.1),
    psp.sliderConf('n'),
]



#region B1A3
def calcCoeff(a0: complex, a1: complex, terms: int, par):
    a2 = -par['alpha']  * a0
    a3 = -(par['alpha'] - 1) * a1 / 3
    a = [a0, a1, a2, a3]
    for i in range (2,terms+1):
        a.append( (i*(3-i)*a[i-2] + 2*(i*i- par['alpha'] )*a[i]) / ( (i+2) * (i+1)))
    return a

def psi(x, a, par):
    # par['µ']  ist b (Breite)
    return polynom(a, np.tanh(x/par['µ'] ))
    
#Jan
# def alphaCalc(E, n):
#     return E + 1 - 1/2 * n**2

#Johannes
def alphaCalc(E, n):
    return E + 1
#endregion 

#region B7A2
# def calcCoeff(a0:complex,a1:complex,terms: int, par):
#     a = [a0]
#     for i in range (1,terms+1):
#         a.append((par['µ'] +i-par['n'] ) / (par['µ'] +i/2) * a[i-1]/i)
#         if a[i] == 0:
#             break
#     return a

# def psi(x, a, par):
#     return polynom(a, par['alpha'] /(2*par['n'] +1) * x) * x**par['µ']  * np.exp(-par['alpha'] /(2*par['n'] +1)*x)
#endregion

#region B8A2
# def calcCoeff(a0: complex, a1: complex, terms: int, par):
#     # par['µ']  ist γ in der Aufgabenbearbeitung (x-Richtung-Wellenzahl) # par['alpha']  ist a0 in der Aufgabenbearbeitung ~ B (Magnetfeld) (E ist B hier)
#     A1 = 2j/3 * np.sqrt(par['alpha'] *par['µ'] )
#     A2 = np.sqrt(par['µ'] *par['µ'] -par['alpha'] /2)
#     a = [a0, 0, 0, -A1/3 * a0, 0]
#     for i in range (5,terms+1):
#         a.append( (-4*(i-2)*a[i-2] + (6*i-15) * A1 * a[i-3] + 12 * A1 * A2 * a[i-5]) / i / (3-2*i) )
#     return a

# def psi(x, a, par):
#     A1 = 2j/3 * np.sqrt(par['alpha'] *par['µ'] )
#     A2 = np.sqrt(par['µ'] *par['µ'] -par['alpha'] /2)

#     return np.exp(- A2 * x + A1 * x**(3/2)) * polynom(a, np.sqrt(x))

# def alphaCalc(E,n):
#     return E
# #endregion

def updateExtraParams(par):
    par['alpha'] = alphaCalc(par['E'],par['n'])


psp.setup(
    calculateCoefficients = calcCoeff,
    ψ = psi,
    updateParameters = updateExtraParams,
    paramSliderSetup = paramSliderSetup)

psp.start()