import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import pyvista as pv
import os

def binomial_coeff(n,k):
    ans = math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
    return ans

def derivative(coefficients):
    length = len(coefficients)
    new_coefficients = coefficients[0:length-1]
    for i in range(1, len(coefficients)):
        new_coefficients[i-1] = coefficients[i] * i

    return new_coefficients

def legendre_polynomial(l):
    # generate the legendre_polynomials 
    coeff = []
    for i in range(0,l+1):
        coeff.append(binomial_coeff(l,i) * (-1)**(i+l) * 1 / (2**l*math.factorial(l)))
        if i<l:
            coeff.append(0)
    coeff = derivative_degree(coeff,l)

    return coeff

def derivative_degree(c, deg):
    if deg == 0:
        return c
    else:
        return derivative_degree(derivative(c), deg - 1)
    
def associated_legendre(l, m, x):
    # x is cos(theta) from a meshgrid of theta,phi
    # associated legedre is a function of m and legendra polynomial of x
    v = np.power(1 - np.power(x, 2), m / 2)
    coefficients = legendre_polynomial(l)
    coefficients = derivative_degree(coefficients, m)
    LEGX = np.zeros_like(x)
    
    for i in range(len(coefficients)):
        LEGX += coefficients[i] * np.power(x, i)

    result = v * LEGX * (-1)**m

    return result

def sphere_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r *np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

  
def Spherical_harmonics(theta, phi, l, m, radius):
    polar = associated_legendre(l, abs(m), np.cos(theta))
    angular = (-1)**m * np.cos(abs(m) * phi)
    if m < 0:
        angular = np.sin(abs(m) * phi)
        polar *= (-1)**m * math.factorial(l-m) / math.factorial(l+m)
        
    spherical_harmonic =  angular * polar
    amplitude = np.max(spherical_harmonic)
    spherical_harmonic /= amplitude

    if radius == 0:
        x,y,z = sphere_to_cartesian(np.abs(spherical_harmonic + radius),theta,phi)
    else:
        x,y,z = sphere_to_cartesian(spherical_harmonic + radius,theta,phi)

    magnitude = pv.PolyData(np.column_stack([x.flatten(),y.flatten(), z.flatten()]))
    magnitude['Color'] = (spherical_harmonic.flatten()+1 - radius)/2
    return magnitude


if __name__ == "__main__":

    resolution = 500
    max_l = 7
    mesh_container = []
    theta = np.linspace(0, math.pi, resolution)
    phi = np.linspace(0, 2*math.pi, resolution)  
    Theta, Phi = np.meshgrid(theta,phi)
    for l in range(0, max_l+1):
        for m in range(-l, l+1):
            smesh = Spherical_harmonics(Theta, Phi, l, m, 2)
            mesh_container.append(smesh)
            name = f'mesh_{l},{m}.vtk'
            smesh.save(name)
    
    plotter = pv.Plotter(shape=(max_l + 1, 2 * max_l + 1), window_size=[1000,1200])
    index = 0
    for l in range(0, max_l+1):
        for m in range(-l, l+1):
            plotter.subplot(l, m + max_l)
            plotter.add_mesh(mesh_container[index], cmap = 'seismic' ,scalars = 'Color' , opacity = 1.0, show_scalar_bar=False)
            plotter.camera.azimuth = 30
            plotter.camera.zoom(2.0)
            index+=1
            

    plotter.show()

