n_timesteps = 5550
dt = 0.1

with open("solution.pvd", "w") as f:
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('  <Collection>\n')
    for i in range(n_timesteps):
        time = i * dt
        filename = f"Solution/solution-{i}.vtu"
        f.write(f'    <DataSet timestep="{time}" group="" part="0" file="{filename}"/>\n')
    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')
