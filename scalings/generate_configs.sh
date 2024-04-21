#!/bin/bash

# Define the ranges and step sizes for each parameter
declare -a temp_range=(-10 -5 0 10)
declare -a wall_thickness_range=(0.01 0.04 0.07 0.1)
declare -a comp_power_range=(100 400 700 1000)

# Create a directory for the configuration files
mkdir -p run

# Function to generate a YAML configuration file
generate_yaml() {
    local temp=$1
    local thickness=$2
    local power=$3
    local filename="run/scaling_temp${temp}_thickness${thickness}_power${power}.yaml"

    echo "Nx: 20" > $filename
    echo "Ny: 36" >> $filename
    echo "Nz: 40" >> $filename
    echo "nGhosts_cond: 1" >> $filename
    echo "nGhosts_euler: 2" >> $filename
    echo "wall_thickness: $thickness" >> $filename
    echo "computer_power: $power" >> $filename
    echo "outside_temperature: $temp" >> $filename
    echo "lamp_power: 500.0" >> $filename
    echo "lamp_dimensions:" >> $filename
    echo "  len_x: 8" >> $filename
    echo "  len_y: 12" >> $filename
    echo "  len_z: 5" >> $filename
    echo "lamp_position:" >> $filename
    echo "  x: 3" >> $filename
    echo "  y: 12" >> $filename
    echo "  z: 30" >> $filename
}

# Generate the configuration files
for temp in "${temp_range[@]}"; do
    for thickness in "${wall_thickness_range[@]}"; do
        for power in "${comp_power_range[@]}"; do
            generate_yaml $temp $thickness $power
        done
    done
done

echo "Configuration files generated in config_files directory."
