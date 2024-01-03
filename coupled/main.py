import driver as dr 
import numpy as np 

c_to_k = float(273.15)
rs = 287.0

def initialize_temperature(T: np.ndarray, outside_temperature: float, Nx: int, Ny: int, Nz: int, nGhosts: int, from_save=False) -> np.ndarray: 
    if not from_save: 
        T = np.full((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), outside_temperature, dtype='f4') # initialize temperature grid 
        T[nGhosts:Nx+nGhosts, nGhosts:Ny+nGhosts, nGhosts:Nz+nGhosts] = outside_temperature + 1.0 # Set the plastic shell to 1.0
        T[nGhosts+1:Nx+nGhosts-1, nGhosts+1:Ny+nGhosts-1, nGhosts+1:Nz+nGhosts-1] = outside_temperature + 2.0 # Set internal 
    else: 
        T[:nGhosts, :, :] = outside_temperature
        T[-nGhosts:, :, :] = outside_temperature
        T[:, :nGhosts, :] = outside_temperature
        T[:, -nGhosts:, :] = outside_temperature
        T[:, :, :nGhosts] = outside_temperature
        T[:, :, -nGhosts:] = outside_temperature
    return T
def initialize_sources(sources: np.ndarray, computer_power: float, lamp_power: float, 
                       lamp_len_x: int, lamp_len_y: int, lamp_len_z: 
                       int, lamp_x: int, lamp_y: int, lamp_z: int, 
                       Nx: int, Ny: int, Nz: int, nGhosts: int, ds: float) -> np.ndarray: 
    lamp_power_density = lamp_power / ((lamp_len_x*lamp_len_y*lamp_len_z)*ds**3)
    computer_power_density_per_grid = float(computer_power / ((Ny/2)*ds**3))

    sources[nGhosts:Nx//2, :nGhosts+1,  nGhosts:] = computer_power_density_per_grid

    x_range = np.arange(Nx + 2*nGhosts)
    y_range = np.arange(Ny + 2*nGhosts)
    z_range = np.arange(Nz + 2*nGhosts)

    # Using broadcasting to create a 3D boolean mask
    lamp_mask = ((x_range[:, None, None] >= lamp_x) & (x_range[:, None, None] < lamp_x + lamp_len_x) &
                (y_range[None, :, None] >= lamp_y) & (y_range[None, :, None] < lamp_y + lamp_len_y) &
                (z_range[None, None, :] >= lamp_z) & (z_range[None, None, :] < lamp_z + lamp_len_z))

    # Apply the lamp power density to the lamp area
    sources[lamp_mask] += lamp_power_density
    return sources

def set_all_arrays_as_fortran(T, sources, rho, p, vx, vy, vz) -> list: 
    vx = np.asfortranarray(vx, dtype='f4')
    vy = np.asfortranarray(vy, dtype='f4')
    vz = np.asfortranarray(vz, dtype='f4')
    rho = np.asfortranarray(rho, dtype='f4')
    p = np.asfortranarray(p, dtype='f4')
    T = np.asfortranarray(T, dtype='f4')
    sources = np.asfortranarray(sources, dtype='f4')
    return T, sources, rho, p, vx, vy, vz

def save_arrays(T, rho, p, vx, vy, vz, i, results_dir='./results/'): 
    for name, arr in zip(['T', 'rho', 'p', 'vx', 'vy', 'vz'], [T, rho, p, vx, vy, vz]): 
        fname = results_dir + f'{name}_{i}'
        np.save(fname, arr)
import yaml 
import sys 
def read_yaml_input_file(fname: str) -> dict: 
    with open(fname, 'r') as f: 
        return yaml.safe_load(f)

def initialize_everything(config: dict, nGhosts_cond: int, nGhosts_euler: int): 

    Nx = int(config['Nx'])
    Ny = int(config['Ny'])
    Nz = int(config['Nz'])
    lamp_power = float(config['lamp_power'])
    lamp_len_x = config['lamp_dimensions']['len_x']
    lamp_len_y = config['lamp_dimensions']['len_y']
    lamp_len_z = config['lamp_dimensions']['len_z']
    lamp_x = config['lamp_position']['x']
    lamp_y = config['lamp_position']['y']
    lamp_z = config['lamp_position']['z']
    outside_temperature = float(config['outside_temperature'])
    computer_power = float(config['computer_power'])
    ds = float(config['ds'])

    T_cond = np.zeros((Nx+2*nGhosts_cond, Ny+2*nGhosts_cond, Nz+2*nGhosts_cond), dtype='f4', order='F') # initialize temperature grid  
    T_euler = np.zeros((Nx+2*nGhosts_euler, Ny+2*nGhosts_euler, Nz+2*nGhosts_euler), dtype='f4', order='F') # initialize temperature grid  
    T_euler[nGhosts_euler:Nx+nGhosts_euler, nGhosts_euler:Ny+nGhosts_euler, nGhosts_euler:Nz+nGhosts_euler] = T_cond[1:-1, 1:-1, 1:-1]
    sources = np.zeros((Nx+2*nGhosts_cond, Ny+2*nGhosts_cond, Nz+2*nGhosts_cond), dtype='f4', order='F')
    T_cond = initialize_temperature(T_cond, outside_temperature, Nx, Ny, Nz, nGhosts_cond)
    sources = initialize_sources(sources, computer_power, lamp_power, lamp_len_x, lamp_len_y, lamp_len_z, lamp_x, lamp_y, lamp_z, Nx, Ny, Nz, nGhosts_cond, ds)

    rho = 1.225 # kg/m^3
    rho = np.full((Nx+2*nGhosts_euler, Ny+2*nGhosts_euler, Nz+2*nGhosts_euler), rho, dtype='f4', order='F')
    p = (T_euler+c_to_k)*rs*rho # Pa

    vx = np.zeros((Nx+2*nGhosts_euler, Ny+2*nGhosts_euler, Nz+2*nGhosts_euler), dtype='f4', order='F')
    vy = np.zeros((Nx+2*nGhosts_euler, Ny+2*nGhosts_euler, Nz+2*nGhosts_euler), dtype='f4', order='F')
    vz = np.zeros((Nx+2*nGhosts_euler, Ny+2*nGhosts_euler, Nz+2*nGhosts_euler), dtype='f4', order='F')
    return T_cond, T_euler, sources, rho, p, vx, vy, vz
def drive():
    config = read_yaml_input_file(sys.argv[1])
    Nx = int(config['Nx'])
    Ny = int(config['Ny'])
    Nz = int(config['Nz'])
    wall_thickness = float(config['wall_thickness'])
    nGhosts_cond, nGhosts_euler = int(1), int(2)
    dt = float(0.1)
    ds = float(0.05)
    
    time = 0.0
    time_max = 1.0

    # setup arrays 
    T_cond, T_euler, sources, rho, p, vx, vy, vz = initialize_everything(config, nGhosts_cond, nGhosts_euler)

    save_arrays(T_cond, rho, p, vx, vy, vz, 0)
    for i in range(1, 15): 
        print(f'Iteration {i}')

        T_cond, sources, rho, p, vx, vy, vz = set_all_arrays_as_fortran(T_cond, sources, rho, p, vx, vy, vz)  
        print('Condctive step')
        dr.run_conductive(T_cond, Nx, Ny, Nz, nGhosts_cond, dt, ds, sources, wall_thickness, time_max)
        
        # copy the T cond temperature except for the border 
        T_euler[nGhosts_euler:Nx+nGhosts_euler, nGhosts_euler:Ny+nGhosts_euler, nGhosts_euler:Nz+nGhosts_euler] = T_cond[1:-1, 1:-1, 1:-1]
        rho = p / ((T_euler+c_to_k)*rs)
        dr.euler_driver(rho,vx,vy,vz,p,Nx,Ny,Nz,nGhosts_euler,ds,time,time_max)

        # update temperature or pressure before passing to the next iteration 
        T_euler = p/(rs*rho) - c_to_k
        T_cond[1:-1, 1:-1, 1:-1] = T_euler[nGhosts_euler:Nx+nGhosts_euler, nGhosts_euler:Ny+nGhosts_euler, nGhosts_euler:Nz+nGhosts_euler]

        # save arrays 
        save_arrays(T_cond, rho, p, vx, vy, vz, i)
if __name__ == '__main__': 
    drive()