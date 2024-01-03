import concurrent.futures
import subprocess
import os
import sys 
def run_simulation(config_file):
    """Run a single simulation with the given configuration file."""
    subprocess.run(['python3', 'main.py', config_file])

if __name__ == "__main__":
    config_files = sys.argv[1:]

    # Run simulations in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(run_simulation, config_files)