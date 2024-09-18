import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import pyvista as pv

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
    coeff = []
    for i in range(0,l+1):
        coeff.append(binomial_coeff(l,i) * (-1)**(i+l) * 1 / (2**l*math.factorial(l)))
        if i<l:
            coeff.append(0)
    coeff = derivative_degree(coeff,l)

    return coeff

def derivative_degree(c, deg):
    for i in range(0,deg):
        c = derivative(c)

    return c

def associated_legendre(l, m, x):
    v = np.power(1 - np.power(x, 2), m / 2)
    coefficients = legendre_polynomial(l)
    coefficients = derivative_degree(coefficients, m)
    LEGX = np.zeros_like(x)
    
    for i in range(len(coefficients)):
        LEGX += coefficients[i] * np.power(x, i)

    result = v * LEGX * (-1)**m
    print(l,m)
    print(LEGX)
    return result

def sphere_to_cartesian(r, theta, phi):
    x = r * math.sin(phi) * math.cos(theta)
    y = r * math.sin(phi) * math.sin(theta)
    z = r * math.cos(phi)
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
    spherical_harmonic += radius
    if radius == 0:
        x = np.abs(spherical_harmonic) * np.sin(theta) * np.cos(phi)
        y = np.abs(spherical_harmonic) * np.sin(theta) * np.sin(phi)
        z = np.abs(spherical_harmonic) * np.cos(theta)
    else:
        x = spherical_harmonic * np.sin(theta) * np.cos(phi)
        y = spherical_harmonic * np.sin(theta) * np.sin(phi)
        z = spherical_harmonic * np.cos(theta)

    points = np.column_stack([x.flatten(),y.flatten(), z.flatten()])
    magnitude = pv.PolyData(points)
    spherical_harmonic -= radius
    magnitude['Color'] = (spherical_harmonic.flatten()+1)/2
    return magnitude


if __name__ == "__main__":

    resolution = 500
    max_l = 2
    plotter = pv.Plotter(shape=(max_l + 1, 2 * max_l + 1), window_size=[1000,1200])
    mesh_container = []
    for l in range(0, max_l+1):
        for m in range(-l, l+1):
            '''
            theta = np.linspace(0, math.pi, resolution)
            phi = np.linspace(0, 2*math.pi, resolusotion)
            y = associated_legendre(l,abs(m),np.cos(theta), resolution)
            x = (-1)**m*np.cos(abs(m) * phi)
            if m < 0:
                x = np.sin(abs(m) * phi)
                y *= (-1)**m * math.factorial(l-m)/math.factorial(l+m)
            
            amplitudes = np.outer(y,x)
            amplitude = np.max(amplitudes)
            amplitudes = amplitudes / amplitude
            Theta, Phi = np.meshgrid(theta,phi)
            

            #print(np.min(amplitudes))
            #print(np.max(amplitudes))
            colors = amplitudes.flatten()
            colmap = [(0, 'blue'), (0.5, 'white'), (1, 'red')]
            cmap1 = LinearSegmentedColormap.from_list('custom', colmap, N = 256)
            amplitudes = (amplitudes + 2) 

            x = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.cos(phi)).flatten()
            y = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.sin(phi)).flatten()
            z = np.array(np.abs(amplitudes) * np.cos(theta[:, np.newaxis])).flatten()
 
            pos = plotter.camera.position
            

            smesh = pv.StructuredGrid(x,y,z)
            '''
            theta = np.linspace(0, math.pi, resolution)
            phi = np.linspace(0, 2*math.pi, resolution)  
            Theta, Phi = np.meshgrid(theta,phi)
            smesh = Spherical_harmonics(Theta, Phi, l, m, 2)


            
            colmap = [(0, 'blue'), (0.5, 'white'), (1, 'red')]
            cmap1 = LinearSegmentedColormap.from_list('custom', colmap, N = 256)

            
            
            #plotter.add_points(np.c_[x,y,z], scalars = colors , cmap = 'seismic', point_size = 0.5, show_scalar_bar = False)
            plotter.subplot(l, m + max_l)
            plotter.add_mesh(smesh, cmap = 'seismic' ,scalars = 'Color' , opacity = 1.0, show_scalar_bar=False)
            
            #plotter.camera_position = (pos[0]*1/4, pos[1]*1/4, pos[2] * 1/4)



    #def update_time():
        
    plotter.render()
    plotter.show()



    '''
    plotter = pv.Plotter(shape=(max_l + 1, 2 * max_l + 2))

    for l in range(0, max_l+1):
        for m in range(0, l+1):
            theta = np.linspace(0, math.pi, resolution)
            phi = np.linspace(0, 2*math.pi, resolution)
            y = associated_legendre(l,abs(m),np.cos(theta), resolution)
            x = (-1)**m*np.cos(abs(m) * phi)
            if m < 0:
                x = np.sin(abs(m) * phi)
                y *= (-1)**m * math.factorial(l-m)/math.factorial(l+m)
            
            amplitudes = np.outer(y,x)
            amplitude = np.max(amplitudes)
            amplitudes = amplitudes / amplitude
            amplitudes1 = (amplitudes + 2) 

            x = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.cos(phi)).flatten()
            y = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.sin(phi)).flatten()
            z = np.array(np.abs(amplitudes) * np.cos(theta[:, np.newaxis])).flatten()

            xh = np.array(np.abs(amplitudes1) * np.sin(theta[:, np.newaxis]) * np.cos(phi)).flatten()
            yh = np.array(np.abs(amplitudes1) * np.sin(theta[:, np.newaxis]) * np.sin(phi)).flatten()
            zh = np.array(np.abs(amplitudes1) * np.cos(theta[:, np.newaxis])).flatten()

            colors = amplitudes.flatten()
            smesh = pv.StructuredGrid(x,y,z)
            hmesh = pv.StructuredGrid(xh,yh,zh)
            plotter.subplot(l, m + max_l+1 )
            plotter.add_mesh(hmesh, cmap = 'seismic',scalars = colors, opacity = 1.0, show_scalar_bar=False)


            #plotter.add_points(np.c_[x,y,z], scalars = colors , cmap = 'seismic', point_size = 0.5, show_scalar_bar = False)
            plotter.subplot(l,max_l - m)
            plotter.add_mesh(smesh, cmap = 'seismic',scalars = colors, opacity = 1.0, show_scalar_bar=False)


    plotter.show()



    '''

    '''
    l  = 3
    m = 2
    plotter = pv.Plotter(window_size=[1000,1200])
    theta = np.linspace(0, math.pi, resolution)
    phi = np.linspace(0, 2*math.pi, resolution)
    y = associated_legendre(l,abs(m),np.cos(theta), resolution)
    x = (-1)**m*np.cos(abs(m) * phi)
    if m < 0:
        x = np.sin(abs(m) * phi)
        y *= (-1)**m * math.factorial(l-m)/math.factorial(l+m)
    
    amplitudes = np.outer(y,x)
    amplitude = np.max(amplitudes)
    amplitudes = amplitudes / amplitude
    amplitudes = (amplitudes * math.sin(math.pi*5/6) + 2) 

    x = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.cos(phi)).flatten()
    y = np.array(np.abs(amplitudes) * np.sin(theta[:, np.newaxis]) * np.sin(phi)).flatten()
    z = np.array(np.abs(amplitudes) * np.cos(theta[:, np.newaxis])).flatten()
    
    colors = amplitudes.flatten()
    smesh = pv.StructuredGrid(x,y,z)
    #plotter.add_points(np.c_[x,y,z], scalars = colors , cmap = 'seismic', point_size = 0.5, show_scalar_bar = False)
    plotter.add_mesh(smesh, cmap = 'seismic',scalars = colors, opacity = 1.0, show_scalar_bar=False)
    
    camr = 8
    camphi = math.pi/3
    camtheta = 0
    camx, camy, camz = sphere_to_cartesian(camr, camtheta, camphi)
    plotter.camera_position = (camx,camy, camz)
    print(plotter.camera.position)
    plotter.render()
    plotter.show()
    '''
    