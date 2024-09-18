import pyvista as pv


if __name__ == '__main__':
    max_l = 4
    mesh_container = []
    for l in range(0, max_l+1):
        for m in range(-l, l+1):
            mesh_container.append(pv.read(f'mesh/mesh_{l},{m}.vtk'))
            
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