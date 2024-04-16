#[allow(dead_code)]
pub mod physics {

    use crate::globals::globals;
    use crate::{Gradients, Primitives, Simulation, SimulationGrid};
    use ndarray::{s, Array3, Array4, Zip};
    use ndarray_stats::QuantileExt;

    pub fn update_ghost_cells(mesh: &mut Array3<f32>, ghosts: usize) {
        let shape = mesh.shape();
        let nx = shape[0];
        let ny = shape[1];
        let nz = shape[2];
        assert!(ghosts == 2);

        {
            let x0 = mesh.slice(s![nx - ghosts - 2, .., ..]).to_owned();
            let x1 = mesh.slice(s![nx - ghosts - 1, .., ..]).to_owned();
            mesh.slice_mut(s![0, .., ..]).assign(&x0);
            mesh.slice_mut(s![1, .., ..]).assign(&x1);
        }

        {
            let x0 = mesh.slice(s![2, .., ..]).to_owned();
            let x1 = mesh.slice(s![3, .., ..]).to_owned();
            mesh.slice_mut(s![nx - 2, .., ..]).assign(&x0);
            mesh.slice_mut(s![nx - 1, .., ..]).assign(&x1);
        }

        {
            let y0 = mesh.slice(s![.., ny - ghosts - 2, ..]).to_owned();
            let y1 = mesh.slice(s![.., ny - ghosts - 1, ..]).to_owned();
            mesh.slice_mut(s![.., 0, ..]).assign(&y0);
            mesh.slice_mut(s![.., 1, ..]).assign(&y1);
        }

        {
            let y0 = mesh.slice(s![.., 2, ..]).to_owned();
            let y1 = mesh.slice(s![.., 3, ..]).to_owned();
            mesh.slice_mut(s![.., ny - 2, ..]).assign(&y0);
            mesh.slice_mut(s![.., ny - 1, ..]).assign(&y1);
        }

        {
            let z0 = mesh.slice(s![.., .., nz - ghosts - 2]).to_owned();
            let z1 = mesh.slice(s![.., .., nz - ghosts - 1]).to_owned();
            mesh.slice_mut(s![.., .., 0]).assign(&z0);
            mesh.slice_mut(s![.., .., 1]).assign(&z1);
        }

        {
            let z0 = mesh.slice(s![.., .., 2]).to_owned();
            let z1 = mesh.slice(s![.., .., 3]).to_owned();
            mesh.slice_mut(s![.., .., nz - 2]).assign(&z0);
            mesh.slice_mut(s![.., .., nz - 1]).assign(&z1);
        }
    }

    pub fn update_ghost_cells4(mesh: &mut Array4<f32>, ghosts: usize) {
        let shape = mesh.shape();
        let nx = shape[0];
        let ny = shape[1];
        let nz = shape[2];
        assert!(ghosts == 2);

        {
            let x0 = mesh.slice(s![nx - ghosts - 2, .., .., ..]).to_owned();
            let x1 = mesh.slice(s![nx - ghosts - 1, .., .., ..]).to_owned();
            mesh.slice_mut(s![0, .., .., ..]).assign(&x0);
            mesh.slice_mut(s![1, .., .., ..]).assign(&x1);
        }

        {
            let x0 = mesh.slice(s![2, .., .., ..]).to_owned();
            let x1 = mesh.slice(s![3, .., .., ..]).to_owned();
            mesh.slice_mut(s![nx - 2, .., .., ..]).assign(&x0);
            mesh.slice_mut(s![nx - 1, .., .., ..]).assign(&x1);
        }

        {
            let y0 = mesh.slice(s![.., ny - ghosts - 2, .., ..]).to_owned();
            let y1 = mesh.slice(s![.., ny - ghosts - 1, .., ..]).to_owned();
            mesh.slice_mut(s![.., 0, .., ..]).assign(&y0);
            mesh.slice_mut(s![.., 1, .., ..]).assign(&y1);
        }

        {
            let y0 = mesh.slice(s![.., 2, .., ..]).to_owned();
            let y1 = mesh.slice(s![.., 3, .., ..]).to_owned();
            mesh.slice_mut(s![.., ny - 2, .., ..]).assign(&y0);
            mesh.slice_mut(s![.., ny - 1, .., ..]).assign(&y1);
        }

        {
            let z0 = mesh.slice(s![.., .., nz - ghosts - 2, ..]).to_owned();
            let z1 = mesh.slice(s![.., .., nz - ghosts - 1, ..]).to_owned();
            mesh.slice_mut(s![.., .., 0, ..]).assign(&z0);
            mesh.slice_mut(s![.., .., 1, ..]).assign(&z1);
        }

        {
            let z0 = mesh.slice(s![.., .., 2, ..]).to_owned();
            let z1 = mesh.slice(s![.., .., 3, ..]).to_owned();
            mesh.slice_mut(s![.., .., nz - 2, ..]).assign(&z0);
            mesh.slice_mut(s![.., .., nz - 1, ..]).assign(&z1);
        }
    }

    pub fn calc_primitives(grid: &mut SimulationGrid) {
        let cell_volume = grid.info.ds.powi(3);

        Zip::from(&mut grid.primitives.rho)
            .and(&grid.conserved.mass)
            .par_for_each(|a, b| *a = b / cell_volume);

        Zip::from(&mut grid.primitives.vx)
            .and(&grid.conserved.momentum_x)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b / (c * cell_volume));

        Zip::from(&mut grid.primitives.vy)
            .and(&grid.conserved.momentum_y)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b / (c * cell_volume));

        Zip::from(&mut grid.primitives.vz)
            .and(&grid.conserved.momentum_z)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b / (c * cell_volume));

        Zip::from(&mut grid.primitives.p)
            .and(&grid.conserved.energy)
            .and(&grid.primitives.rho)
            .and(&grid.primitives.vx)
            .and(&grid.primitives.vy)
            .and(&grid.primitives.vz)
            .par_for_each(|p, e, r, vx, vy, vz| {
                *p = (e / cell_volume - 0.5 * r * (vx * vx + vy * vy + vz * vz))
                    * (globals::GAMMA - 1.0)
            });
    }

    pub fn calc_conserved(grid: &mut SimulationGrid) {
        let cell_volume = grid.info.ds.powi(3);

        Zip::from(&mut grid.conserved.mass)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b| *a = b * cell_volume);

        Zip::from(&mut grid.conserved.momentum_x)
            .and(&grid.primitives.vx)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b * c * cell_volume);

        Zip::from(&mut grid.conserved.momentum_y)
            .and(&grid.primitives.vy)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b * c * cell_volume);

        Zip::from(&mut grid.conserved.momentum_z)
            .and(&grid.primitives.vz)
            .and(&grid.primitives.rho)
            .par_for_each(|a, b, c| *a = b * c * cell_volume);

        Zip::from(&mut grid.conserved.energy)
            .and(&grid.primitives.rho)
            .and(&grid.primitives.vx)
            .and(&grid.primitives.vy)
            .and(&grid.primitives.vz)
            .and(&grid.primitives.p)
            .par_for_each(|e, r, vx, vy, vz, p| {
                *e = cell_volume
                    * (r * (vx * vx + vy * vy + vz * vz) / 2.0 + p / (globals::GAMMA - 1.0))
            });
    }

    pub fn calc_timestep(grid: &SimulationGrid) -> f32 {
        let ds = grid.info.ds;
        let mut temp = grid.primitives.p.clone();

        Zip::from(&mut temp)
            .and(&grid.primitives.rho)
            .and(&grid.primitives.vx)
            .and(&grid.primitives.vy)
            .and(&grid.primitives.vz)
            .par_for_each(|v, r, vx, vy, vz| {
                *v = ds
                    / (f32::sqrt(globals::GAMMA * (*v) / r)
                        + f32::sqrt(vx * vx + vy * vy + vz * vz))
            });
        let min_val = match temp.min() {
            Ok(v) => v.to_owned(),
            Err(_err) => 1.0,
        };
        globals::CFL * min_val
    }

    pub fn calc_gradients(mesh: &Array3<f32>, gradient: &mut Array4<f32>, info: &Simulation) {
        let ds = info.ds;
        let (nx, ny, nz) = (mesh.shape()[0], mesh.shape()[1], mesh.shape()[2]);
        Zip::indexed(gradient).par_for_each(|(i, j, k, dim), val| {
            if i >= 1 && i < nx - 1 {
                if j >= 1 && j < ny - 1 {
                    if k >= 1 && k < nz - 1 {
                        *val = match dim {
                            0 => (mesh[[i + 1, j, k]] - mesh[[i - 1, j, k]]) / (2.0 * ds),
                            1 => (mesh[[i, j + 1, k]] - mesh[[i, j - 1, k]]) / (2.0 * ds),
                            2 => (mesh[[i, j, k + 1]] - mesh[[i, j, k - 1]]) / (2.0 * ds),
                            _ => panic!(),
                        };
                    }
                }
            }
        });
    }

    pub fn extrapolate_primitives_t(
        p_xtr: &mut Primitives,
        p: &Primitives,
        g: &Gradients,
        dt: f32,
    ) {
        let (nx, ny, nz) = (p.rho.shape()[0], p.rho.shape()[1], p.rho.shape()[2]);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let rho = p.rho[[i, j, k]];
                    let vx = p.vx[[i, j, k]];
                    let vy = p.vy[[i, j, k]];
                    let vz = p.vz[[i, j, k]];

                    let drho_dx = g.drho[[i, j, k, 0]];
                    let dvx_dx = g.dvx[[i, j, k, 0]];
                    let dvy_dx = g.dvy[[i, j, k, 0]];
                    let dvz_dx = g.dvz[[i, j, k, 0]];

                    let drho_dy = g.drho[[i, j, k, 1]];
                    let dvx_dy = g.dvx[[i, j, k, 1]];
                    let dvy_dy = g.dvy[[i, j, k, 1]];
                    let dvz_dy = g.dvz[[i, j, k, 1]];

                    let drho_dz = g.drho[[i, j, k, 2]];
                    let dvx_dz = g.dvx[[i, j, k, 2]];
                    let dvy_dz = g.dvy[[i, j, k, 2]];
                    let dvz_dz = g.dvz[[i, j, k, 2]];

                    let dp_dx = g.dp[[i, j, k, 0]];
                    let dp_dy = g.dp[[i, j, k, 1]];
                    let dp_dz = g.dp[[i, j, k, 2]];

                    p_xtr.rho[[i, j, k]] = rho
                        - 0.5
                            * dt
                            * (vx * drho_dx
                                + rho * dvx_dx
                                + vy * drho_dy
                                + rho * dvy_dy
                                + vz * drho_dz
                                + rho * dvz_dz);
                    p_xtr.vx[[i, j, k]] = vx
                        - 0.5
                            * dt
                            * (vx * dvx_dx + vy * dvx_dy + vz * dvx_dz + (1.0 / rho) * dp_dx);
                    p_xtr.vy[[i, j, k]] = vy
                        - 0.5
                            * dt
                            * (vx * dvy_dx + vy * dvy_dy + vz * dvy_dz + (1.0 / rho) * dp_dy);
                    p_xtr.vz[[i, j, k]] = vz
                        - 0.5
                            * dt
                            * (vx * dvz_dx + vy * dvz_dy + vz * dvz_dz + (1.0 / rho) * dp_dz);
                    p_xtr.p[[i, j, k]] = p.p[[i, j, k]]
                        - 0.5
                            * dt
                            * (globals::GAMMA * p.p[[i, j, k]] * (dvx_dx + dvy_dy + dvz_dz)
                                + vx * dp_dx
                                + vy * dp_dy
                                + vz * dp_dz);
                }
            }
        }
    }

    pub fn calc_flux_reconstruction(
        info: &Simulation,
        mass_flux_x: &mut Array4<f32>,
        momentum_x_flux_x: &mut Array4<f32>,
        momentum_y_flux_x: &mut Array4<f32>,
        momentum_z_flux_x: &mut Array4<f32>,
        energy_flux_x: &mut Array4<f32>,
        drho: &Array4<f32>,
        dvx: &Array4<f32>,
        dvy: &Array4<f32>,
        dvz: &Array4<f32>,
        dp: &Array4<f32>,
        rho: &Array3<f32>,
        vx: &Array3<f32>,
        vy: &Array3<f32>,
        vz: &Array3<f32>,
        p: &Array3<f32>,
        offsets: [usize; 3],
        dim: usize,
    ) {
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;
        let ds = info.ds;
        let ghosts = info.nghosts;

        for i in 1..(nx - 1) {
            for j in 1..(ny - 1) {
                for k in 1..(nz - 1) {
                    let rl = rho[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                        - (drho[[i + offsets[0], j + offsets[1], k + offsets[2], dim]]) * (ds / 2.);

                    let rr = rho[[i, j, k]] + (drho[[i, j, k, dim]]) * (ds / 2.);

                    let vxl = vx[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                        - (dvx[[i + offsets[0], j + offsets[1], k + offsets[2], dim]]) * (ds / 2.);

                    let vxr = vx[[i, j, k]] + (dvx[[i, j, k, dim]]) * (ds / 2.);

                    let vyl = vy[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                        - (dvy[[i + offsets[0], j + offsets[1], k + offsets[2], dim]]) * (ds / 2.);

                    let vyr = vy[[i, j, k]] + (dvy[[i, j, k, dim]]) * (ds / 2.);

                    let vzl = vz[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                        - (dvz[[i + offsets[0], j + offsets[1], k + offsets[2], dim]]) * (ds / 2.);

                    let vzr = vz[[i, j, k]] + (dvz[[i, j, k, dim]]) * (ds / 2.);

                    let pl = p[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                        - (dp[[i + offsets[0], j + offsets[1], k + offsets[2], dim]]) * (ds / 2.);

                    let pr = p[[i, j, k]] + (dp[[i, j, k, dim]]) * (ds / 2.);

                    let en_left = 0.5 * (rl * (vxl.powi(2) + vyl.powi(2) + vzl.powi(2)))
                        + (pl / (globals::GAMMA - 1.0));
                    let en_right = 0.5 * (rr * (vxr.powi(2) + vyr.powi(2) + vzr.powi(2)))
                        + (pr / (globals::GAMMA - 1.0));

                    let rho_star = (rl + rr) / 2.0;
                    let momentum_x_star = (rl * vxl + rr * vxr) / 2.0;
                    let momentum_y_star = (rl * vyl + rr * vyr) / 2.0;
                    let momentum_z_star = (rl * vzl + rr * vzr) / 2.0;
                    let en_star = 0.5 * (en_left + en_right);
                    let p_star = (globals::GAMMA - 1.0)
                        * (en_star
                            - 0.5
                                * (momentum_x_star.powi(2)
                                    + momentum_y_star.powi(2)
                                    + momentum_z_star.powi(2))
                                / rho_star);

                    mass_flux_x[[i, j, k, dim]] = momentum_x_star;
                    momentum_x_flux_x[[i, j, k, dim]] =
                        (momentum_x_star.powi(2) / rho_star) + p_star;
                    momentum_y_flux_x[[i, j, k, dim]] =
                        momentum_x_star * momentum_y_star / rho_star;
                    momentum_z_flux_x[[i, j, k, dim]] =
                        momentum_x_star * momentum_z_star / rho_star;
                    energy_flux_x[[i, j, k, dim]] =
                        (en_star + p_star) * (momentum_x_star / rho_star);

                    //Gravity term if enabled
                    if offsets[2] == 1 && info.gravity == true {
                        let h = (nz as f32 - 2.0 * ghosts as f32) * ds - (k as f32 - 2.0) * ds;
                        momentum_x_flux_x[[i, j, k, dim]] += 1.0 * rho_star * globals::G * h;
                        energy_flux_x[[i, j, k, dim]] +=
                            0.5 * rho_star * (vxl + vxr) * globals::G * h;
                    }

                    let gamma = globals::GAMMA;
                    let c_l = (gamma * pl / rl).sqrt() + vxl.abs();
                    let c_r = (gamma * pr / rr).sqrt() + vxr.abs();
                    let c_star = c_l.max(c_r);

                    mass_flux_x[[i, j, k, dim]] = momentum_x_star - 0.5 * c_star * (rl - rr);
                    momentum_x_flux_x[[i, j, k, dim]] -= c_star * (rl * vxl - rr * vxr) / 2.0;
                    momentum_y_flux_x[[i, j, k, dim]] -= c_star * (rl * vyl - rr * vyr) / 2.0;
                    momentum_z_flux_x[[i, j, k, dim]] -= c_star * (rl * vzl - rr * vzr) / 2.0;
                    energy_flux_x[[i, j, k, dim]] -= c_star * (en_left - en_right) / 2.0;
                    mass_flux_x[[i, j, k, dim]] *= ds;
                    momentum_x_flux_x[[i, j, k, dim]] *= ds;
                    momentum_y_flux_x[[i, j, k, dim]] *= ds;
                    momentum_z_flux_x[[i, j, k, dim]] *= ds;
                    energy_flux_x[[i, j, k, dim]] *= ds;
                }
            }
        }
    }

    pub fn _calc_flux_reconstruction(
        info: &Simulation,
        mass_flux_x: &mut Array4<f32>,
        momentum_x_flux_x: &mut Array4<f32>,
        momentum_y_flux_x: &mut Array4<f32>,
        momentum_z_flux_x: &mut Array4<f32>,
        energy_flux_x: &mut Array4<f32>,
        drho: &Array4<f32>,
        dvx: &Array4<f32>,
        dvy: &Array4<f32>,
        dvz: &Array4<f32>,
        dp: &Array4<f32>,
        rho: &Array3<f32>,
        vx: &Array3<f32>,
        vy: &Array3<f32>,
        vz: &Array3<f32>,
        p: &Array3<f32>,
        offsets: [usize; 3],
        dim: usize,
    ) {
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;
        let ds = info.ds;
        let ghosts = info.nghosts;
        let mut mass_flux_view = mass_flux_x.slice_mut(s![.., .., .., dim]);
        let mut momx_flux_view = momentum_x_flux_x.slice_mut(s![.., .., .., dim]);
        let mut momy_flux_view = momentum_y_flux_x.slice_mut(s![.., .., .., dim]);
        let mut momz_flux_view = momentum_z_flux_x.slice_mut(s![.., .., .., dim]);
        let mut e_flux_view = energy_flux_x.slice_mut(s![.., .., .., dim]);
        Zip::indexed(&mut mass_flux_view)
            .and(&mut momx_flux_view)
            .and(&mut momy_flux_view)
            .and(&mut momz_flux_view)
            .and(&mut e_flux_view)
            .par_for_each(|(i, j, k), q1, q2, q3, q4, q5| {
                if i >= 1 && i < nx - 1 {
                    if j >= 1 && j < ny - 1 {
                        if k >= 1 && k < nz - 1 {
                            let rl = rho[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                                - (drho[[i + offsets[0], j + offsets[1], k + offsets[2], dim]])
                                    * (ds / 2.);

                            let rr = rho[[i, j, k]] + (drho[[i, j, k, dim]]) * (ds / 2.);

                            let vxl = vx[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                                - (dvx[[i + offsets[0], j + offsets[1], k + offsets[2], dim]])
                                    * (ds / 2.);

                            let vxr = vx[[i, j, k]] + (dvx[[i, j, k, dim]]) * (ds / 2.);

                            let vyl = vy[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                                - (dvy[[i + offsets[0], j + offsets[1], k + offsets[2], dim]])
                                    * (ds / 2.);

                            let vyr = vy[[i, j, k]] + (dvy[[i, j, k, dim]]) * (ds / 2.);

                            let vzl = vz[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                                - (dvz[[i + offsets[0], j + offsets[1], k + offsets[2], dim]])
                                    * (ds / 2.);

                            let vzr = vz[[i, j, k]] + (dvz[[i, j, k, dim]]) * (ds / 2.);

                            let pl = p[[i + offsets[0], j + offsets[1], k + offsets[2]]]
                                - (dp[[i + offsets[0], j + offsets[1], k + offsets[2], dim]])
                                    * (ds / 2.);

                            let pr = p[[i, j, k]] + (dp[[i, j, k, dim]]) * (ds / 2.);

                            let en_left = 0.5 * (rl * (vxl.powi(2) + vyl.powi(2) + vzl.powi(2)))
                                + (pl / (globals::GAMMA - 1.0));
                            let en_right = 0.5 * (rr * (vxr.powi(2) + vyr.powi(2) + vzr.powi(2)))
                                + (pr / (globals::GAMMA - 1.0));

                            let rho_star = (rl + rr) / 2.0;
                            let momentum_x_star = (rl * vxl + rr * vxr) / 2.0;
                            let momentum_y_star = (rl * vyl + rr * vyr) / 2.0;
                            let momentum_z_star = (rl * vzl + rr * vzr) / 2.0;
                            let en_star = 0.5 * (en_left + en_right);
                            let p_star = (globals::GAMMA - 1.0)
                                * (en_star
                                    - 0.5
                                        * (momentum_x_star.powi(2)
                                            + momentum_y_star.powi(2)
                                            + momentum_z_star.powi(2))
                                        / rho_star);

                            let mut q1_new = momentum_x_star;
                            let mut q2_new = (momentum_x_star.powi(2) / rho_star) + p_star;
                            let mut q3_new = momentum_x_star * momentum_y_star / rho_star;
                            let mut q4_new = momentum_x_star * momentum_z_star / rho_star;
                            let mut q5_new = (en_star + p_star) * (momentum_x_star / rho_star);

                            //Gravity term if enabled
                            if offsets[2] == 1 && info.gravity == true {
                                let h =
                                    (nz as f32 - 2.0 * ghosts as f32) * ds - (k as f32 - 2.0) * ds;
                                q1_new += 1.0 * rho_star * globals::G * h;
                                q5_new += 0.5 * rho_star * (vxl + vxr) * globals::G * h;
                            }

                            let gamma = globals::GAMMA;
                            let c_l = (gamma * pl / rl).sqrt() + vxl.abs();
                            let c_r = (gamma * pr / rr).sqrt() + vxr.abs();
                            let c_star = c_l.max(c_r);

                            q1_new = momentum_x_star - 0.5 * c_star * (rl - rr);
                            q2_new -= c_star * (rl * vxl - rr * vxr) / 2.0;
                            q3_new -= c_star * (rl * vyl - rr * vyr) / 2.0;
                            q4_new -= c_star * (rl * vzl - rr * vzr) / 2.0;
                            q5_new -= c_star * (en_left - en_right) / 2.0;
                            q1_new *= ds;
                            q2_new *= ds;
                            q3_new *= ds;
                            q4_new *= ds;
                            q5_new *= ds;

                            *q1 = q1_new;
                            *q2 = q2_new;
                            *q3 = q3_new;
                            *q4 = q4_new;
                            *q5 = q5_new;
                        }
                    }
                }
            });
    }

    pub fn add_fluxes(info: &Simulation, grid: &mut SimulationGrid, dt: f32) {
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;
        let ds = info.ds;

        for i in 1..(nx - 1) {
            for j in 1..(ny - 1) {
                for k in 1..(nz - 1) {
                    grid.conserved.mass[[i, j, k]] -= dt
                        * ds
                        * (grid.fluxes.mass[[i, j, k, 0]] - grid.fluxes.mass[[i - 1, j, k, 0]]
                            + grid.fluxes.mass[[i, j, k, 1]]
                            - grid.fluxes.mass[[i, j - 1, k, 1]]
                            + grid.fluxes.mass[[i, j, k, 2]]
                            - grid.fluxes.mass[[i, j, k - 1, 2]]);

                    grid.conserved.momentum_x[[i, j, k]] -= dt
                        * ds
                        * (grid.fluxes.momentum_x[[i, j, k, 0]]
                            - grid.fluxes.momentum_x[[i - 1, j, k, 0]]
                            + grid.fluxes.momentum_x[[i, j, k, 1]]
                            - grid.fluxes.momentum_x[[i, j - 1, k, 1]]
                            + grid.fluxes.momentum_x[[i, j, k, 2]]
                            - grid.fluxes.momentum_x[[i, j, k - 1, 2]]);

                    grid.conserved.momentum_y[[i, j, k]] -= dt
                        * ds
                        * (grid.fluxes.momentum_y[[i, j, k, 0]]
                            - grid.fluxes.momentum_y[[i - 1, j, k, 0]]
                            + grid.fluxes.momentum_y[[i, j, k, 1]]
                            - grid.fluxes.momentum_y[[i, j - 1, k, 1]]
                            + grid.fluxes.momentum_y[[i, j, k, 2]]
                            - grid.fluxes.momentum_y[[i, j, k - 1, 2]]);

                    grid.conserved.momentum_z[[i, j, k]] -= dt
                        * ds
                        * (grid.fluxes.momentum_z[[i, j, k, 0]]
                            - grid.fluxes.momentum_z[[i - 1, j, k, 0]]
                            + grid.fluxes.momentum_z[[i, j, k, 1]]
                            - grid.fluxes.momentum_z[[i, j - 1, k, 1]]
                            + grid.fluxes.momentum_z[[i, j, k, 2]]
                            - grid.fluxes.momentum_z[[i, j, k - 1, 2]]);

                    grid.conserved.energy[[i, j, k]] -= dt
                        * ds
                        * (grid.fluxes.energy[[i, j, k, 0]] - grid.fluxes.energy[[i - 1, j, k, 0]]
                            + grid.fluxes.energy[[i, j, k, 1]]
                            - grid.fluxes.energy[[i, j - 1, k, 1]]
                            + grid.fluxes.energy[[i, j, k, 2]]
                            - grid.fluxes.energy[[i, j, k - 1, 2]]);
                }
            }
        }
    }

    pub fn compute_step(info: &Simulation, grid: &mut SimulationGrid, dt: &mut f32) {
        // grid.update_primitive_ghost_cells();
        calc_primitives(grid);
        grid.update_primitive_ghost_cells();
        *dt = calc_timestep(&grid);

        calc_gradients(&grid.primitives.rho, &mut grid.gradients.drho, &info);
        calc_gradients(&grid.primitives.vx, &mut grid.gradients.dvx, &info);
        calc_gradients(&grid.primitives.vy, &mut grid.gradients.dvy, &info);
        calc_gradients(&grid.primitives.vz, &mut grid.gradients.dvz, &info);
        calc_gradients(&grid.primitives.p, &mut grid.gradients.dp, &info);
        grid.update_gradients_ghost_cells();

        extrapolate_primitives_t(
            &mut grid.primitives_xtr,
            &grid.primitives,
            &grid.gradients,
            *dt,
        );
        grid.update_primitive_xtr_ghost_cells();

        calc_flux_reconstruction(
            &info,
            &mut grid.fluxes.mass,
            &mut grid.fluxes.momentum_x,
            &mut grid.fluxes.momentum_y,
            &mut grid.fluxes.momentum_z,
            &mut grid.fluxes.energy,
            &grid.gradients.drho,
            &grid.gradients.dvx,
            &grid.gradients.dvy,
            &grid.gradients.dvz,
            &grid.gradients.dp,
            &grid.primitives_xtr.rho,
            &grid.primitives_xtr.vx,
            &grid.primitives_xtr.vy,
            &grid.primitives_xtr.vz,
            &grid.primitives_xtr.p,
            [1, 0, 0],
            0,
        );

        calc_flux_reconstruction(
            &info,
            &mut grid.fluxes.mass,
            &mut grid.fluxes.momentum_y,
            &mut grid.fluxes.momentum_x,
            &mut grid.fluxes.momentum_z,
            &mut grid.fluxes.energy,
            &grid.gradients.drho,
            &grid.gradients.dvy,
            &grid.gradients.dvx,
            &grid.gradients.dvz,
            &grid.gradients.dp,
            &grid.primitives_xtr.rho,
            &grid.primitives_xtr.vy,
            &grid.primitives_xtr.vx,
            &grid.primitives_xtr.vz,
            &grid.primitives_xtr.p,
            [0, 1, 0],
            1,
        );

        calc_flux_reconstruction(
            &info,
            &mut grid.fluxes.mass,
            &mut grid.fluxes.momentum_z,
            &mut grid.fluxes.momentum_x,
            &mut grid.fluxes.momentum_y,
            &mut grid.fluxes.energy,
            &grid.gradients.drho,
            &grid.gradients.dvz,
            &grid.gradients.dvx,
            &grid.gradients.dvy,
            &grid.gradients.dp,
            &grid.primitives_xtr.rho,
            &grid.primitives_xtr.vz,
            &grid.primitives_xtr.vx,
            &grid.primitives_xtr.vy,
            &grid.primitives_xtr.p,
            [0, 0, 1],
            2,
        );
        grid.update_fluxes_ghost_cells();

        add_fluxes(info, grid, *dt);
    }
}
