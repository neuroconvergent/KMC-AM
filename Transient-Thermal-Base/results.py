import os
import glob
import re
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.offline as pyo


def read_gpl_file(filepath):
    """
    Reads a .gpl file and extracts node positions and temperature.
    Assumes format: x y z temperature (space-separated).
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    data = []
    for line in lines:
        if line.strip().startswith('#') or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) >= 4:
            x, y, z, temp = map(float, parts[:4])
            data.append(temp)
    return data

def collect_temperature_time_series(folder):
    """
    Parses all .gpl files in a folder and builds a list of temperature arrays per time step.
    Returns:
        node_time_series: list of temperature lists per node
        timesteps: list of timestep indices
    """
    file_list = sorted(glob.glob(os.path.join(folder, "solution-*.gpl")),
                       key=lambda x: int(re.search(r"solution-(\d+)\.gpl", x).group(1)))

    timesteps = []
    all_data = []

    for f in file_list:
        timestep = int(re.search(r"solution-(\d+)\.gpl", f).group(1))
        temps = read_gpl_file(f)
        timesteps.append(timestep*0.1)
        all_data.append(temps)

    # Transpose: per-node temperature vs time
    node_time_series = list(map(list, zip(*all_data)))
    return node_time_series, timesteps

def plot_node_temperatures(node_time_series, timesteps, max_plots=12):
    """
    Plot temperature vs. time for each node using Plotly.
    """
    num_nodes = len(node_time_series)
    rows = min(max_plots, num_nodes)
    fig = make_subplots(rows=rows, cols=1,
                        shared_xaxes=True,
                        subplot_titles=[f"Node {i}" for i in range(rows)],
                        vertical_spacing=0.05)

    for i in range(rows):
        trace = go.Scatter(x=timesteps,
                           y=node_time_series[i],
                           mode='lines+markers',
                           name=f'Node {i}')
        fig.add_trace(trace, row=i + 1, col=1)

    fig.update_layout(
        height=500 * rows,
        title_text="Temperature vs. Time for Each Node",
        showlegend=False,
    )
    for i in range(rows):
        fig.update_xaxes(
                showticklabels = True,
                dtick = 55,
                title_text = "Time (s)",
                row=i,
                col=1
                )
        fig.update_yaxes(
                showticklabels = True,
                dtick = 100,
                title_text = "Temp (°C)",
                row=i,
                col=1
                )
    pyo.plot(fig, filename="node_temperatures.html", auto_open=True)
# ---------- Run ----------

folder = "Solution"  # Change if your .gpl files are in another directory
node_time_series, timesteps = collect_temperature_time_series(folder)
plot_node_temperatures(node_time_series, timesteps, max_plots=10)  # Adjust max_plots as needed

