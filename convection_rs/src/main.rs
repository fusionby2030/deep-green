mod boundaries;
mod globals;
mod grid;
mod init;
mod io;
mod physics;
use crate::grid::grid::*;
use crate::init::init::*;
use crate::io::io::*;
use crate::physics::physics::*;
use clap::Parser;
use std::process::ExitCode;
use std::time::Instant;
use strum::IntoEnumIterator;

#[derive(Parser, Debug)]
#[clap(about = "A compressible CFD code to model greenhouses for cheap!")]
#[clap(author = "Kostis Papadakis and Adam Kit 2024")]
#[command(version, about,author, long_about = None)]
struct Args {
    /// Sets the demo project to run
    #[clap(short, long)]
    demo: String,

    /// Sets the output interval
    #[clap(short, long)]
    tout: f32,

    /// Sets the total time to simulate
    #[clap(short, long)]
    total_time: f32,

    /// Sets the total steps to simulate
    #[clap(short, long)]
    max_steps: usize,
}

fn compute(
    info: Simulation,
    grid: &mut SimulationGrid,
    total_time: f32,
    max_steps: usize,
    tout: f32,
) {
    let mut tnow: f32 = 0.0;
    let mut step: usize = 0;
    let mut dt: f32 = 0.5;
    let mut wtime = 0.0;
    let mut wstep = 0;
    calc_conserved(grid);
    loop {
        if tnow > total_time || step >= max_steps {
            break;
        }
        println!("Computing step {step} and time {tnow} with dt {dt}. ");
        compute_step(&info, grid, &mut dt);
        if wtime > tout || wstep == 0 {
            let padded = format!("{:05}", wstep);
            let fname = "state.".to_owned() + &padded.to_string() + ".h5";
            save_to_disk(&grid, &fname);
            wtime = 0.0;
            wstep += 1;
        }
        tnow += dt;
        wtime += dt;
        step += 1;
    }
}

fn main() -> ExitCode {
    let args = Args::parse();

    let total_time = args.total_time;
    let tout = args.tout;
    let max_steps = args.max_steps;
    let demo = args.demo;

    let t0 = Instant::now();
    if demo == "KELVIN_HELMHOLTZ" {
        let (simulation, mut grid) = set_up_khi_demo();
        compute(simulation, &mut grid, total_time, max_steps, tout);
    } else if demo == "THERMAL_RISING_BUBBLE" {
        let (simulation, mut grid) = set_up_trb_demo();
        compute(simulation, &mut grid, total_time, max_steps, tout);
    } else if demo == "THERMAL_RISING_BUBBLES" {
        let (simulation, mut grid) = set_up_trbs_demo();
        compute(simulation, &mut grid, total_time, max_steps, tout);
    } else if demo == "EXPLOSION" {
        let (simulation, mut grid) = set_up_2d_explosion_demo();
        compute(simulation, &mut grid, total_time, max_steps, tout);
    } else if demo == "RELAXATION" {
        let (simulation, mut grid) = set_up_3d_relaxation_demo();
        compute(simulation, &mut grid, total_time, max_steps, tout);
    } else {
        eprintln!("Demo not found! Available demos are:");
        let available_demos: Vec<RunType> = RunType::iter().collect();
        eprintln!("{:?}", available_demos);
    }

    let t1 = t0.elapsed();
    println!("Simulation done in {:2?}!", t1);
    return ExitCode::from(0);
}
