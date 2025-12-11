#Code based on equation (13) in the following review: https://www.mdpi.com/2073-8994/15/11/2025

import matplotlib.pyplot as plt
import numpy as np

beta=0.4
gamma=0.1
D=0.20 #diffusion coefficient

nx, ny=100, 100
susceptible_grid=np.ones((nx, ny))
#modeling portions of Mission, Potrero Hill, and South of Market in San Francisco (https://statisticalatlas.com/place/California/San-Francisco/Population)
mission_scaled=56872/(100*50) #assume population spread evenly throughout neighborhoods
susceptible_grid[:, 0:50]*=mission_scaled
south_scaled=24608/(25*50)
susceptible_grid[0:25, 50:100]*=south_scaled
potrero_scaled=10825/(75*50)
susceptible_grid[25:100, 50:100]*=potrero_scaled


infected_grid=np.zeros((nx, ny)) #represents concentration of people infected
infected_grid[nx//2-6:nx//2+6, ny//2-6:ny//2+6]=10 #people at the center are infected

recovered_grid=np.zeros((nx, ny))
N=np.sum(susceptible_grid+infected_grid+recovered_grid)

#simulation with time steps
visualize_times=[0, 150, 300, 450] #time steps to visualize diffusion
visualize_idx=0 #points to next visualize time
fig, ax=plt.subplots(2, 2)
snapshots=dict()
snapshots_recovered=dict()

for step in range(451):
    new_susceptible_grid, new_infected_grid, new_recovered_grid=susceptible_grid.copy(), infected_grid.copy(), recovered_grid.copy()

    if visualize_idx<len(visualize_times) and step==visualize_times[visualize_idx]: #snapshot the susceptible people at this time
        snapshots[visualize_idx]=new_susceptible_grid.copy()
        snapshots_recovered[visualize_idx]=recovered_grid.copy()
        visualize_idx+=1

    #loop over old grids and update according to diffusion equation w/ laplacian operator
    #assume boundaries are absolute + cannot have diffusion there (eliminates edge cases)
    #laplacian[x][y] is about (grid[x-1][y]+grid[x+1][y]+grid[x][y-1]+grid[x][y+1]-4*grid[x][y])/1.0 (step size=1.0)
    laplacians=(susceptible_grid[2:, 1:-1]+susceptible_grid[0:-2, 1:-1]+susceptible_grid[1:-1, 2:]+susceptible_grid[1:-1, 0:-2]-4*susceptible_grid[1:-1, 1:-1])/1.0 #computes all laplacians at inside points much faster
    new_susceptible_grid[1:-1, 1:-1]=new_susceptible_grid[1:-1, 1:-1]+D*laplacians
    new_susceptible_grid=new_susceptible_grid-beta*susceptible_grid*infected_grid/N #SIR equation (-beta*S*I/N)


    laplacians=(infected_grid[2:, 1:-1]+infected_grid[0:-2, 1:-1]+infected_grid[1:-1, 2:]+infected_grid[1:-1, 0:-2]-4*infected_grid[1:-1, 1:-1])/1.0
    new_infected_grid[1:-1, 1:-1]+=D*laplacians
    new_infected_grid=new_infected_grid+beta*susceptible_grid*infected_grid/N-gamma*infected_grid #SIR equation (beta*S*I/N - gamma*I)

    laplacians=(recovered_grid[2:, 1:-1]+recovered_grid[0:-2, 1:-1]+recovered_grid[1:-1, 2:]+recovered_grid[1:-1, 0:-2]-4*recovered_grid[1:-1, 1:-1])/1.0
    new_recovered_grid[1:-1, 1:-1]+=D*laplacians
    new_recovered_grid=new_recovered_grid+gamma*infected_grid #SIR eqaution (gamma*I)

    #update old grids
    susceptible_grid=new_susceptible_grid.copy()
    infected_grid=new_infected_grid.copy()
    recovered_grid=new_recovered_grid.copy()


#visualization code! - plot three heatmaps for susceptible people
vmin=min(np.min(s) for s in snapshots.values())
vmax=max(np.max(s) for s in snapshots.values())

for i, step in enumerate(visualize_times):
    row, col=i//2, i%2
    ref=ax[row, col].imshow(snapshots[i], cmap='viridis', vmin=vmin, vmax=vmax)
    ax[row, col].set_title(f"Susceptible People: Day {step}", fontsize=10)
    ax[row, col].set_xlabel("X-Coordinate", fontsize=9, labelpad=-1)
    ax[row, col].set_ylabel("Y-Coordinate", fontsize=9)

fig.colorbar(ref, ax=ax.ravel().tolist(), orientation='vertical',
             fraction=0.046, pad=0.046, label="Number of People", location="right")
fig.suptitle("Number of Susceptible Individuals over Time")
plt.show()


#visualization code! - plot three heatmaps for recovered people
vmin=min(np.min(s) for s in snapshots_recovered.values())
vmax=max(np.max(s) for s in snapshots_recovered.values())

for i, step in enumerate(visualize_times):
    row, col=i//2, i%2
    ref=ax[row, col].imshow(snapshots_recovered[i], cmap='viridis', vmin=vmin, vmax=vmax)
    ax[row, col].set_title(f"Recovered People: Day {step}", fontsize=10)
    ax[row, col].set_xlabel("X-Coordinate", fontsize=9, labelpad=-1)
    ax[row, col].set_ylabel("Y-Coordinate", fontsize=9)

fig.colorbar(ref, ax=ax.ravel().tolist(), orientation='vertical',
             fraction=0.046, pad=0.046, label="Number of People", location="right")

fig.suptitle("Number of Recovered People over Time")
plt.show()

