pub mod boundaries {
    use crate::{BoundaryType, Simulation};
    use ndarray::{s, Array3, Array4};

    pub fn apply_boundary_conditions(info: &Simulation, mesh: &mut Array3<f32>) {
        let ghosts = info.nghosts;
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;
        match info.bcs[0] {
            BoundaryType::OUTFLOW => {
                let x0 = mesh.slice(s![2, .., ..]).to_owned();
                mesh.slice_mut(s![0, .., ..]).assign(&x0);
                mesh.slice_mut(s![1, .., ..]).assign(&x0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let x0 = mesh.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                let x1 = mesh.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                mesh.slice_mut(s![0, .., ..]).assign(&x0);
                mesh.slice_mut(s![1, .., ..]).assign(&x1);
            }
        };

        match info.bcs[1] {
            BoundaryType::OUTFLOW => {
                let x0 = mesh.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                mesh.slice_mut(s![nx - 2, .., ..]).assign(&x0);
                mesh.slice_mut(s![nx - 1, .., ..]).assign(&x0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let x0 = mesh.slice(s![2, .., ..]).to_owned();
                let x1 = mesh.slice(s![3, .., ..]).to_owned();
                mesh.slice_mut(s![nx - 2, .., ..]).assign(&x0);
                mesh.slice_mut(s![nx - 1, .., ..]).assign(&x1);
            }
        };

        match info.bcs[2] {
            BoundaryType::OUTFLOW => {
                let y0 = mesh.slice(s![.., 2, ..]).to_owned();
                mesh.slice_mut(s![.., 0, ..]).assign(&y0);
                mesh.slice_mut(s![.., 1, ..]).assign(&y0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let y0 = mesh.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                let y1 = mesh.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                mesh.slice_mut(s![.., 0, ..]).assign(&y0);
                mesh.slice_mut(s![.., 1, ..]).assign(&y1);
            }
        };

        match info.bcs[3] {
            BoundaryType::OUTFLOW => {
                let y0 = mesh.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                mesh.slice_mut(s![.., ny - 2, ..]).assign(&y0);
                mesh.slice_mut(s![.., ny - 1, ..]).assign(&y0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let y0 = mesh.slice(s![.., 2, ..]).to_owned();
                let y1 = mesh.slice(s![.., 3, ..]).to_owned();
                mesh.slice_mut(s![.., ny - 2, ..]).assign(&y0);
                mesh.slice_mut(s![.., ny - 1, ..]).assign(&y1);
            }
        };

        match info.bcs[4] {
            BoundaryType::OUTFLOW => {
                let z0 = mesh.slice(s![.., .., 2]).to_owned();
                mesh.slice_mut(s![.., .., 0]).assign(&z0);
                mesh.slice_mut(s![.., .., 1]).assign(&z0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let z0 = mesh.slice(s![.., .., nz - ghosts - 2]).to_owned();
                let z1 = mesh.slice(s![.., .., nz - ghosts - 1]).to_owned();
                mesh.slice_mut(s![.., .., 0]).assign(&z0);
                mesh.slice_mut(s![.., .., 1]).assign(&z1);
            }
        };

        match info.bcs[5] {
            BoundaryType::OUTFLOW => {
                let z0 = mesh.slice(s![.., .., nz - ghosts - 1]).to_owned();
                mesh.slice_mut(s![.., .., nz - 2]).assign(&z0);
                mesh.slice_mut(s![.., .., nz - 1]).assign(&z0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let z0 = mesh.slice(s![.., .., 2]).to_owned();
                let z1 = mesh.slice(s![.., .., 3]).to_owned();
                mesh.slice_mut(s![.., .., nz - 2]).assign(&z0);
                mesh.slice_mut(s![.., .., nz - 1]).assign(&z1);
            }
        };
    }

    pub fn apply_boundary_conditions4(info: &Simulation, mesh: &mut Array4<f32>) {
        let ghosts = info.nghosts;
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;
        match info.bcs[0] {
            BoundaryType::OUTFLOW => {
                let x0 = mesh.slice(s![2, .., .., ..]).to_owned();
                mesh.slice_mut(s![0, .., .., ..]).assign(&x0);
                mesh.slice_mut(s![1, .., .., ..]).assign(&x0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let x0 = mesh.slice(s![nx - ghosts - 2, .., .., ..]).to_owned();
                let x1 = mesh.slice(s![nx - ghosts - 1, .., .., ..]).to_owned();
                mesh.slice_mut(s![0, .., .., ..]).assign(&x0);
                mesh.slice_mut(s![1, .., .., ..]).assign(&x1);
            }
        };

        match info.bcs[1] {
            BoundaryType::OUTFLOW => {
                let x0 = mesh.slice(s![nx - ghosts - 1, .., .., ..]).to_owned();
                mesh.slice_mut(s![nx - 2, .., .., ..]).assign(&x0);
                mesh.slice_mut(s![nx - 1, .., .., ..]).assign(&x0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let x0 = mesh.slice(s![2, .., .., ..]).to_owned();
                let x1 = mesh.slice(s![3, .., .., ..]).to_owned();
                mesh.slice_mut(s![nx - 2, .., .., ..]).assign(&x0);
                mesh.slice_mut(s![nx - 1, .., .., ..]).assign(&x1);
            }
        };

        match info.bcs[2] {
            BoundaryType::OUTFLOW => {
                let y0 = mesh.slice(s![.., 2, .., ..]).to_owned();
                mesh.slice_mut(s![.., 0, .., ..]).assign(&y0);
                mesh.slice_mut(s![.., 1, .., ..]).assign(&y0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let y0 = mesh.slice(s![.., ny - ghosts - 2, .., ..]).to_owned();
                let y1 = mesh.slice(s![.., ny - ghosts - 1, .., ..]).to_owned();
                mesh.slice_mut(s![.., 0, .., ..]).assign(&y0);
                mesh.slice_mut(s![.., 1, .., ..]).assign(&y1);
            }
        };

        match info.bcs[3] {
            BoundaryType::OUTFLOW => {
                let y0 = mesh.slice(s![.., ny - ghosts - 1, .., ..]).to_owned();
                mesh.slice_mut(s![.., ny - 2, .., ..]).assign(&y0);
                mesh.slice_mut(s![.., ny - 1, .., ..]).assign(&y0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let y0 = mesh.slice(s![.., 2, .., ..]).to_owned();
                let y1 = mesh.slice(s![.., 3, .., ..]).to_owned();
                mesh.slice_mut(s![.., ny - 2, .., ..]).assign(&y0);
                mesh.slice_mut(s![.., ny - 1, .., ..]).assign(&y1);
            }
        };

        match info.bcs[4] {
            BoundaryType::OUTFLOW => {
                let z0 = mesh.slice(s![.., .., 2, ..]).to_owned();
                mesh.slice_mut(s![.., .., 0, ..]).assign(&z0);
                mesh.slice_mut(s![.., .., 1, ..]).assign(&z0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let z0 = mesh.slice(s![.., .., nz - ghosts - 2, ..]).to_owned();
                let z1 = mesh.slice(s![.., .., nz - ghosts - 1, ..]).to_owned();
                mesh.slice_mut(s![.., .., 0, ..]).assign(&z0);
                mesh.slice_mut(s![.., .., 1, ..]).assign(&z1);
            }
        };

        match info.bcs[5] {
            BoundaryType::OUTFLOW => {
                let z0 = mesh.slice(s![.., .., nz - ghosts - 1, ..]).to_owned();
                mesh.slice_mut(s![.., .., nz - 2, ..]).assign(&z0);
                mesh.slice_mut(s![.., .., nz - 1, ..]).assign(&z0);
            }
            BoundaryType::WALL => {}
            BoundaryType::PERIODIC => {
                let z0 = mesh.slice(s![.., .., 2, ..]).to_owned();
                let z1 = mesh.slice(s![.., .., 3, ..]).to_owned();
                mesh.slice_mut(s![.., .., nz - 2, ..]).assign(&z0);
                mesh.slice_mut(s![.., .., nz - 1, ..]).assign(&z1);
            }
        };
    }

    pub fn handle_walls(
        info: &Simulation,
        rho: &mut Array3<f32>,
        vx: &mut Array3<f32>,
        vy: &mut Array3<f32>,
        vz: &mut Array3<f32>,
        p: &mut Array3<f32>,
    ) {
        let ghosts = info.nghosts;
        let nx = info.nx;
        let ny = info.ny;
        let nz = info.nz;

        match info.bcs[0] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![2, .., ..]).to_owned();
                    let mut vxt = vx.slice(s![2, .., ..]).to_owned();
                    let vyt = vy.slice(s![2, .., ..]).to_owned();
                    let vzt = vz.slice(s![2, .., ..]).to_owned();
                    let pt = p.slice(s![2, .., ..]).to_owned();
                    vxt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![1, .., ..]).assign(&rt);
                    vx.slice_mut(s![1, .., ..]).assign(&vxt);
                    vy.slice_mut(s![1, .., ..]).assign(&vyt);
                    vz.slice_mut(s![1, .., ..]).assign(&vzt);
                    p.slice_mut(s![1, .., ..]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![3, .., ..]).to_owned();
                    let mut vxt = vx.slice(s![3, .., ..]).to_owned();
                    let vyt = vy.slice(s![3, .., ..]).to_owned();
                    let vzt = vz.slice(s![3, .., ..]).to_owned();
                    let pt = p.slice(s![3, .., ..]).to_owned();
                    vxt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![0, .., ..]).assign(&rt);
                    vx.slice_mut(s![0, .., ..]).assign(&vxt);
                    vy.slice_mut(s![0, .., ..]).assign(&vyt);
                    vz.slice_mut(s![0, .., ..]).assign(&vzt);
                    p.slice_mut(s![0, .., ..]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };

        match info.bcs[1] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                    let mut vxt = vx.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                    let vyt = vy.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                    let vzt = vz.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                    let pt = p.slice(s![nx - ghosts - 1, .., ..]).to_owned();
                    vxt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![nx - ghosts, .., ..]).assign(&rt);
                    vx.slice_mut(s![nx - ghosts, .., ..]).assign(&vxt);
                    vy.slice_mut(s![nx - ghosts, .., ..]).assign(&vyt);
                    vz.slice_mut(s![nx - ghosts, .., ..]).assign(&vzt);
                    p.slice_mut(s![nx - ghosts, .., ..]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                    let mut vxt = vx.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                    let vyt = vy.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                    let vzt = vz.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                    let pt = p.slice(s![nx - ghosts - 2, .., ..]).to_owned();
                    vxt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![nx - 1, .., ..]).assign(&rt);
                    vx.slice_mut(s![nx - 1, .., ..]).assign(&vxt);
                    vy.slice_mut(s![nx - 1, .., ..]).assign(&vyt);
                    vz.slice_mut(s![nx - 1, .., ..]).assign(&vzt);
                    p.slice_mut(s![nx - 1, .., ..]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };

        match info.bcs[2] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![.., 2, ..]).to_owned();
                    let vxt = vx.slice(s![.., 2, ..]).to_owned();
                    let mut vyt = vy.slice(s![.., 2, ..]).to_owned();
                    let vzt = vz.slice(s![.., 2, ..]).to_owned();
                    let pt = p.slice(s![.., 2, ..]).to_owned();
                    vyt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., 1, ..]).assign(&rt);
                    vx.slice_mut(s![.., 1, ..]).assign(&vxt);
                    vy.slice_mut(s![.., 1, ..]).assign(&vyt);
                    vz.slice_mut(s![.., 1, ..]).assign(&vzt);
                    p.slice_mut(s![.., 1, ..]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![.., 3, ..]).to_owned();
                    let vxt = vx.slice(s![.., 3, ..]).to_owned();
                    let mut vyt = vy.slice(s![.., 3, ..]).to_owned();
                    let vzt = vz.slice(s![.., 3, ..]).to_owned();
                    let pt = p.slice(s![.., 3, ..]).to_owned();
                    vyt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., 0, ..]).assign(&rt);
                    vx.slice_mut(s![.., 0, ..]).assign(&vxt);
                    vy.slice_mut(s![.., 0, ..]).assign(&vyt);
                    vz.slice_mut(s![.., 0, ..]).assign(&vzt);
                    p.slice_mut(s![.., 0, ..]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };

        match info.bcs[3] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                    let vxt = vx.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                    let mut vyt = vy.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                    let vzt = vz.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                    let pt = p.slice(s![.., ny - ghosts - 1, ..]).to_owned();
                    vyt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., ny - ghosts, ..]).assign(&rt);
                    vx.slice_mut(s![.., ny - ghosts, ..]).assign(&vxt);
                    vy.slice_mut(s![.., ny - ghosts, ..]).assign(&vyt);
                    vz.slice_mut(s![.., ny - ghosts, ..]).assign(&vzt);
                    p.slice_mut(s![.., ny - ghosts, ..]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                    let vxt = vx.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                    let mut vyt = vy.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                    let vzt = vz.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                    let pt = p.slice(s![.., ny - ghosts - 2, ..]).to_owned();
                    vyt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., ny - 1, ..]).assign(&rt);
                    vx.slice_mut(s![.., ny - 1, ..]).assign(&vxt);
                    vy.slice_mut(s![.., ny - 1, ..]).assign(&vyt);
                    vz.slice_mut(s![.., ny - 1, ..]).assign(&vzt);
                    p.slice_mut(s![.., ny - 1, ..]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };

        match info.bcs[4] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![.., .., 2]).to_owned();
                    let vxt = vx.slice(s![.., .., 2]).to_owned();
                    let vyt = vy.slice(s![.., .., 2]).to_owned();
                    let mut vzt = vz.slice(s![.., .., 2]).to_owned();
                    let pt = p.slice(s![.., .., 2]).to_owned();
                    vzt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., .., 1]).assign(&rt);
                    vx.slice_mut(s![.., .., 1]).assign(&vxt);
                    vy.slice_mut(s![.., .., 1]).assign(&vyt);
                    vz.slice_mut(s![.., .., 1]).assign(&vzt);
                    p.slice_mut(s![.., .., 1]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![.., .., 3]).to_owned();
                    let vxt = vx.slice(s![.., .., 3]).to_owned();
                    let vyt = vy.slice(s![.., .., 3]).to_owned();
                    let mut vzt = vz.slice(s![.., .., 3]).to_owned();
                    let pt = p.slice(s![.., .., 3]).to_owned();
                    vzt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., .., 0]).assign(&rt);
                    vx.slice_mut(s![.., .., 0]).assign(&vxt);
                    vy.slice_mut(s![.., .., 0]).assign(&vyt);
                    vz.slice_mut(s![.., .., 0]).assign(&vzt);
                    p.slice_mut(s![.., .., 0]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };

        match info.bcs[5] {
            BoundaryType::WALL => {
                {
                    let rt = rho.slice(s![.., .., nz - ghosts - 1]).to_owned();
                    let vxt = vx.slice(s![.., .., nz - ghosts - 1]).to_owned();
                    let vyt = vy.slice(s![.., .., nz - ghosts - 1]).to_owned();
                    let mut vzt = vz.slice(s![.., .., nz - ghosts - 1]).to_owned();
                    let pt = p.slice(s![.., .., nz - ghosts - 1]).to_owned();
                    vzt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., .., nz - ghosts]).assign(&rt);
                    vx.slice_mut(s![.., .., nz - ghosts]).assign(&vxt);
                    vy.slice_mut(s![.., .., nz - ghosts]).assign(&vyt);
                    vz.slice_mut(s![.., .., nz - ghosts]).assign(&vzt);
                    p.slice_mut(s![.., .., nz - ghosts]).assign(&pt);
                }
                {
                    let rt = rho.slice(s![.., .., nz - ghosts - 2]).to_owned();
                    let vxt = vx.slice(s![.., .., nz - ghosts - 2]).to_owned();
                    let vyt = vy.slice(s![.., .., nz - ghosts - 2]).to_owned();
                    let mut vzt = vz.slice(s![.., .., nz - ghosts - 2]).to_owned();
                    let pt = p.slice(s![.., .., nz - ghosts - 2]).to_owned();
                    vzt.mapv_inplace(|x| -1.0 * x);
                    rho.slice_mut(s![.., .., nz - 1]).assign(&rt);
                    vx.slice_mut(s![.., .., nz - 1]).assign(&vxt);
                    vy.slice_mut(s![.., .., nz - 1]).assign(&vyt);
                    vz.slice_mut(s![.., .., nz - 1]).assign(&vzt);
                    p.slice_mut(s![.., .., nz - 1]).assign(&pt);
                }
            }
            BoundaryType::PERIODIC => {}
            BoundaryType::OUTFLOW => {}
        };
    }
}
