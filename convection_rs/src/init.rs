#[allow(dead_code)]
pub mod init {
    use ndarray::s;

    use crate::{
        calc_conserved, globals::globals, BoundaryType, RunType, Simulation, SimulationGrid,
    };

    pub fn set_up_khi_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
        ];
        let nx = 256;
        let ny = 2;
        let nz = 128;
        let nghosts = 2;
        let ds = 1.0 / nx as f32;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::KELVIN_HELMHOLTZ,
            bcs,
            [true, true, true],
            nghosts,
            false,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_khi(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn set_up_explosion_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
            BoundaryType::PERIODIC,
        ];
        let nx = 64;
        let ny = 64;
        let nz = 64;
        let nghosts = 2;
        let ds = 1.0 / nx as f32;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::EXPLOSION,
            bcs,
            [true, true, true],
            nghosts,
            false,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_explosion(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn set_up_2d_explosion_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::WALL,     //x-
            BoundaryType::WALL,     //x+
            BoundaryType::PERIODIC, //y-
            BoundaryType::PERIODIC, //y+
            BoundaryType::WALL,     //z-
            BoundaryType::WALL,     //z+
        ];
        let nx = 128;
        let ny = 2;
        let nz = 128;
        let nghosts = 2;
        let ds = 1.0 / nx as f32;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::EXPLOSION,
            bcs,
            [false, true, false],
            nghosts,
            false,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_explosion_2d(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn set_up_3d_relaxation_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::WALL,    //x-
            BoundaryType::WALL,    //x+
            BoundaryType::WALL,    //y-
            BoundaryType::WALL,    //y+
            BoundaryType::WALL,    //z-
            BoundaryType::OUTFLOW, //z+
        ];
        let nx = 30;
        let ny = 30;
        let nz = 30;
        let nghosts = 2;
        let ds = 0.1;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::RELAXATION,
            bcs,
            [false, false, false],
            nghosts,
            true,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_relaxation_3d(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn set_up_trb_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::WALL,     //x-
            BoundaryType::WALL,     //x+
            BoundaryType::PERIODIC, //y-
            BoundaryType::PERIODIC, //y+
            BoundaryType::WALL,     //z-
            BoundaryType::OUTFLOW,  //z+
        ];
        let nx = 150;
        let ny = 2;
        let nz = 512;
        let nghosts = 2;
        let ds = 1.0;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::THERMAL_RISING_BUBBLE,
            bcs,
            [false, true, false],
            nghosts,
            true,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_trb_demo(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn set_up_trbs_demo() -> (Simulation, SimulationGrid) {
        let bcs = [
            BoundaryType::WALL,     //x-
            BoundaryType::WALL,     //x+
            BoundaryType::PERIODIC, //y-
            BoundaryType::PERIODIC, //y+
            BoundaryType::WALL,     //z-
            BoundaryType::OUTFLOW,  //z+
        ];
        let nx = 128;
        let ny = 2;
        let nz = 500;
        let nghosts = 2;
        let ds = 2.0;

        let simulation: Simulation = Simulation::new(
            nx,
            ny,
            nz,
            ds,
            RunType::THERMAL_RISING_BUBBLES,
            bcs,
            [false, true, false],
            nghosts,
            true,
        )
        .unwrap();
        let mut grid: SimulationGrid = SimulationGrid::new(simulation);
        println!("Running simulation {:?}.", simulation);
        init_trbs_demo(&mut grid);
        calc_conserved(&mut grid);
        grid.update_primitive_ghost_cells();
        (simulation, grid)
    }

    pub fn init_uniform(grid: &mut SimulationGrid) {
        grid.primitives.rho.map_inplace(|x| *x = globals::DENSITY);
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);
        grid.primitives.p.map_inplace(|x| *x = globals::PRESSURE);
    }

    pub fn init_uniform_split(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let xpos = i as f32 * ds;
                    let length = nx as f32 * ds;
                    grid.primitives.rho[[i, j, k]] = if xpos < 0.3 * length || xpos > 0.7 * length {
                        1.
                    } else {
                        2.
                    };
                }
            }
        }
        grid.primitives.vx.map_inplace(|x| *x = 1000.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);
        grid.primitives.p.map_inplace(|x| *x = 2.5);
    }

    pub fn init_khi(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;

        let kappa = 0.1;
        let s = 0.05 / f32::sqrt(2.0);
        let logic2dbl = |condition: bool| -> f32 {
            if condition {
                1.0
            } else {
                0.0
            }
        };

        grid.primitives.p.map_inplace(|x| *x = 2.5);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let xs = i as f32 * ds;
                    let ys = 0.0 as f32 * ds;
                    let zs = k as f32 * ds;

                    grid.primitives.rho[[i, j, k]] = 1.0
                        + logic2dbl(f32::abs(zs - ds * nz as f32 / 2.0) < ds * nz as f32 / 4.0)
                        + 0.0 * ys;

                    grid.primitives.vx[[i, j, k]] =
                        -0.5 + logic2dbl(f32::abs(zs - ds * nz as f32 / 2.) < ds * nz as f32 / 4.0);

                    grid.primitives.vz[[i, j, k]] = kappa
                        * f32::sin(4.0 * std::f32::consts::PI * xs)
                        * (f32::exp(
                            -f32::powf(zs - ds * nz as f32 / 4.0, 2.0) / (2.0 * s.powi(2)),
                        ) + f32::exp(
                            -f32::powf(zs - 3.0 * ds * nz as f32 / 4.0, 2.0) / (2.0 * s.powi(2)),
                        ));
                }
            }
        }
    }

    pub fn init_explosion(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;

        grid.primitives.rho.map_inplace(|x| *x = 1.0);
        grid.primitives.p.map_inplace(|x| *x = 2.5);
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);

        let gaussian_3d = |i: f32, j: f32, k: f32, sigma: f32| -> f32 {
            let a: f32 = 1.0;
            let exponent = -(i.powi(2) + j.powi(2) + k.powi(2)) / (2.0 * sigma.powi(2));
            let result = a * (exponent).exp();
            result
        };

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let xs = ((i as f32 - nx as f32 / 2.) as f32 * ds) / nx as f32;
                    let ys = ((j as f32 - ny as f32 / 2.) as f32 * ds) / nx as f32;
                    let zs = ((k as f32 - nz as f32 / 2.) as f32 * ds) / nx as f32;
                    grid.primitives.p[[i, j, k]] += 0.1 * gaussian_3d(xs, ys, zs, 0.001);
                }
            }
        }
    }

    pub fn init_explosion_2d(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;

        grid.primitives.rho.map_inplace(|x| *x = 1.0);
        grid.primitives.p.map_inplace(|x| *x = 2.5);
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);

        let gaussian_3d = |i: f32, j: f32, k: f32, sigma: f32| -> f32 {
            let a: f32 = 1.0;
            let exponent = -(i.powi(2) + j.powi(2) + k.powi(2)) / (2.0 * sigma.powi(2));
            let result = a * (exponent).exp();
            result
        };

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let xs = ((i as f32 - nx as f32 / 2.) as f32 * ds) / nx as f32;
                    let ys = ((j as f32 - ny as f32 / 2.) as f32 * ds) / nx as f32;
                    let zs = ((k as f32 - nz as f32 / 2.) as f32 * ds) / nx as f32;
                    grid.primitives.p[[i, j, k]] += 0.1 * gaussian_3d(xs, ys, zs, 0.0001);
                }
            }
        }
    }

    pub fn init_relaxation_3d(grid: &mut SimulationGrid) {
        grid.primitives.rho.map_inplace(|x| *x = globals::DENSITY);
        grid.primitives.p.map_inplace(|x| *x = globals::PRESSURE);
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);
    }

    pub fn init_trb_demo(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;
        let ghosts = grid.info.nghosts;
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);
        let p0 = 101000.0;
        let t0 = 293.0;

        //First we build an artificial equilibrium.
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let zs: f32 = (nz as f32 - 2.0 * ghosts as f32) * ds as f32 + (0.0 * ds)
                        - (k as f32 - 2.0) * ds;
                    grid.primitives.p[[i, j, k]] = p0 + globals::G.abs() * zs;
                    grid.primitives.rho[[i, j, k]] =
                        grid.primitives.p[[i, j, k]] / (globals::RS * t0);
                }
            }
        }
        //Slices
        let mut ss = grid.primitives.p.slice(s![.., .., 2]).to_owned();
        let mut ss2 = grid.primitives.p.slice(s![.., .., 2]).to_owned();
        ss.mapv_inplace(|x| x + globals::G.abs() * ds);
        ss2.mapv_inplace(|x| x + 2.0 * globals::G.abs() * ds);
        grid.primitives.p.slice_mut(s![.., .., 1]).assign(&ss);
        grid.primitives.p.slice_mut(s![.., .., 0]).assign(&ss2);

        //Then we add the bubble;
        let rad_outter = 25.0;
        let rad_inner = 20.0;
        let center = [75.0, 1.0, 128.0];

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mag = f32::sqrt(
                        (i as f32 - center[0]).powi(2)
                            + (j as f32 - center[1]).powi(2)
                            + (k as f32 - center[2]).powi(2),
                    );
                    if mag < rad_outter {
                        grid.primitives.rho[[i, j, k]] = 1.0;
                    }
                    if mag < rad_inner {
                        grid.primitives.rho[[i, j, k]] = 0.8;
                    }
                }
            }
        }
    }

    pub fn init_trbs_demo(grid: &mut SimulationGrid) {
        let nx = grid.info.nx;
        let ny = grid.info.ny;
        let nz = grid.info.nz;
        let ds = grid.info.ds;
        let ghosts = grid.info.nghosts;
        grid.primitives.vx.map_inplace(|x| *x = 0.0);
        grid.primitives.vy.map_inplace(|x| *x = 0.0);
        grid.primitives.vz.map_inplace(|x| *x = 0.0);
        let p0 = 101000.0;
        let t0 = 293.0;

        //First we build an artificial equilibrium.
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let zs: f32 = (nz as f32 - 2.0 * ghosts as f32) * ds as f32 + (0.0 * ds)
                        - (k as f32 - 2.0) * ds;
                    grid.primitives.p[[i, j, k]] = p0 + globals::G.abs() * zs;
                    grid.primitives.rho[[i, j, k]] =
                        grid.primitives.p[[i, j, k]] / (globals::RS * t0);
                }
            }
        }
        //Slices
        let mut ss = grid.primitives.p.slice(s![.., .., 2]).to_owned();
        let mut ss2 = grid.primitives.p.slice(s![.., .., 2]).to_owned();
        ss.mapv_inplace(|x| x + globals::G.abs() * ds);
        ss2.mapv_inplace(|x| x + 2.0 * globals::G.abs() * ds);
        grid.primitives.p.slice_mut(s![.., .., 1]).assign(&ss);
        grid.primitives.p.slice_mut(s![.., .., 0]).assign(&ss2);

        //Then we add the bubble;
        let rad_outter = 15.0;
        let rad_inner = 10.0;
        let center1 = [32.0 + 8.0, 1.0, 32.0];
        let center2 = [88.0 - 8.0, 1.0, 32.0];
        let center3 = [64.0, 1.0, 100.0];

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mag = f32::sqrt(
                        (i as f32 - center1[0]).powi(2)
                            + (j as f32 - center1[1]).powi(2)
                            + (k as f32 - center1[2]).powi(2),
                    );
                    if mag < rad_outter {
                        grid.primitives.rho[[i, j, k]] = 1.0;
                    }
                    if mag < rad_inner {
                        grid.primitives.rho[[i, j, k]] = 0.8;
                    }
                }
            }
        }

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mag = f32::sqrt(
                        (i as f32 - center2[0]).powi(2)
                            + (j as f32 - center2[1]).powi(2)
                            + (k as f32 - center2[2]).powi(2),
                    );
                    if mag < rad_outter {
                        grid.primitives.rho[[i, j, k]] = 1.0;
                    }
                    if mag < rad_inner {
                        grid.primitives.rho[[i, j, k]] = 0.8;
                    }
                }
            }
        }

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mag = f32::sqrt(
                        (i as f32 - center3[0]).powi(2)
                            + (j as f32 - center3[1]).powi(2)
                            + (k as f32 - center3[2]).powi(2),
                    );
                    if mag < 1.5 * rad_outter {
                        grid.primitives.rho[[i, j, k]] = 1.0;
                    }
                    if mag < 1.5 * rad_inner {
                        grid.primitives.rho[[i, j, k]] = 0.8;
                    }
                }
            }
        }
    }
}
