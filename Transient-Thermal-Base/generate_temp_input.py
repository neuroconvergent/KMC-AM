import matplotlib.pyplot as plt
import numpy as np

# Parameters
mass = 0.1                # kg
cp = 500                  # J/kg·K
area = 0.025                # m²
h = 10                    # W/m²·K
epsilon = 0.3             # Emissivity
sigma = 5.670e-8          # W/m²·K⁴
T_amb = 20 + 273.15       # Ambient temp in K

# Time settings
dt = 0.1
t_heating = 5
t_cooling = 550
t_total = t_heating + t_cooling

times = []
temperatures = []

# Heating: linear ramp from 25 to 900 °C over t_heating
for t in np.arange(0, t_heating + dt, dt):
    T_C = 25 + (900 - 25) * (t / t_heating)
    times.append(t)
    temperatures.append(T_C)

# Cooling phase
T = 900 + 273.15  # Kelvin initial temp after heating

for step in np.arange(dt, t_cooling + dt, dt):
    t = t_heating + step
    # Calculate heat losses in Kelvin
    q_conv = h * area * (T - T_amb)
    q_rad = epsilon * sigma * area * (T**4 - T_amb**4)
    q_total = q_conv + q_rad
    dTdt = -q_total / (mass * cp)
    T += dTdt * dt
    times.append(t)
    temperatures.append(T - 273.15)  # Store Celsius for plotting

# Save data
with open("temperature_input.dat", "w") as f:
    f.write("# Time(s)\tTemperature(C)\n")
    for t, T_C in zip(times, temperatures):
        f.write(f"{t:.2f}\t{T_C:.2f}\n")

# Plot
plt.figure(figsize=(10, 5))
plt.plot(times, temperatures, label="Heating + Cooling")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.title("Thermal Profile")
plt.grid(True)
plt.axvline(x=t_heating, color='gray', linestyle='--', label="End of Heating")
plt.legend()
plt.tight_layout()
plt.savefig("temperature_profile.png", dpi=300)
plt.show()
