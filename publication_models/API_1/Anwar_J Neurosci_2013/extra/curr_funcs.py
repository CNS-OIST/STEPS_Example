
import math

F = 96485
R = 8.314


def getGHKI(P, V, z, T, Si, So):
    vs = (z*V*F)/(R*T)
    den = P*z*F*vs*(Si-(So*math.exp(-vs)))
    num = 1-math.exp(-vs)
    return (den/num)

def getOhmI(V, Erev, g):
    return g*(V-Erev)

