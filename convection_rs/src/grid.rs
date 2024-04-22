#[allow(dead_code)]
#[allow(non_camel_case_types)]
pub mod grid {

    use ndarray::{Array3, Array4};

    use crate::boundaries::boundaries::apply_boundary_conditions;
    use crate::boundaries::boundaries::apply_boundary_conditions4;
    use crate::boundaries::boundaries::handle_walls;
    use strum_macros::EnumIter;

    #[derive(Debug, Copy, Clone, EnumIter)]
    pub enum RunType {
        RELAXATION,
        EXPLOSION,
        KELVIN_HELMHOLTZ,
        THERMAL_RISING_BUBBLE,
        THERMAL_RISING_BUBBLES,
        SHOCK_TUBE,
        RAYLEIGH_TAYLOR,
    }

    #[derive(Debug, Copy, Clone, PartialEq, Eq)]
    pub enum BoundaryType {
        PERIODIC,
        OUTFLOW,
        WALL,
    }

    #[derive(Debug, Copy, Clone)]
    pub struct Simulation {
        pub nx: usize,
        pub ny: usize,
        pub nz: usize,
        pub ds: f32,
        pub run_type: RunType,
        pub bcs: [BoundaryType; 6],
        periodicity: [bool; 3],
        pub nghosts: usize,
        pub gravity: bool,
    }

    impl Simulation {
        pub fn new(
            nx: usize,
            ny: usize,
            nz: usize,
            ds: f32,
            run_type: RunType,
            bcs: [BoundaryType; 6],
            periodicity: [bool; 3],
            nghosts: usize,
            gravity: bool,
        ) -> Option<Self> {
            if periodicity[0] {
                assert!(bcs[0] == BoundaryType::PERIODIC && bcs[1] == BoundaryType::PERIODIC);
            }
            if periodicity[1] {
                assert!(bcs[2] == BoundaryType::PERIODIC && bcs[3] == BoundaryType::PERIODIC);
            }
            if periodicity[2] {
                assert!(bcs[4] == BoundaryType::PERIODIC && bcs[5] == BoundaryType::PERIODIC);
            }
            let sim = Simulation {
                nx: nx + 2 * nghosts,
                ny: ny + 2 * nghosts,
                nz: nz + 2 * nghosts,
                ds,
                run_type,
                bcs,
                periodicity,
                nghosts,
                gravity,
            };
            Some(sim)
        }
    }

    pub struct Fluxes {
        pub mass: Array4<f32>,
        pub momentum_x: Array4<f32>,
        pub momentum_y: Array4<f32>,
        pub momentum_z: Array4<f32>,
        pub energy: Array4<f32>,
    }

    pub struct Conserved {
        pub mass: Array3<f32>,
        pub momentum_x: Array3<f32>,
        pub momentum_y: Array3<f32>,
        pub momentum_z: Array3<f32>,
        pub energy: Array3<f32>,
    }

    pub struct Primitives {
        pub rho: Array3<f32>,
        pub vx: Array3<f32>,
        pub vy: Array3<f32>,
        pub vz: Array3<f32>,
        pub p: Array3<f32>,
    }

    pub struct Gradients {
        pub drho: Array4<f32>,
        pub dvx: Array4<f32>,
        pub dvy: Array4<f32>,
        pub dvz: Array4<f32>,
        pub dp: Array4<f32>,
    }

    pub struct SimulationGrid {
        pub info: Simulation,
        pub primitives: Primitives,
        pub primitives_xtr: Primitives,
        pub conserved: Conserved,
        pub gradients: Gradients,
        pub fluxes: Fluxes,
    }

    impl SimulationGrid {
        pub fn update_primitive_ghost_cells(&mut self) {
            handle_walls(
                &self.info,
                &mut self.primitives.rho,
                &mut self.primitives.vx,
                &mut self.primitives.vy,
                &mut self.primitives.vz,
                &mut self.primitives.p,
            );
            apply_boundary_conditions(&self.info, &mut self.primitives.rho);
            apply_boundary_conditions(&self.info, &mut self.primitives.vx);
            apply_boundary_conditions(&self.info, &mut self.primitives.vy);
            apply_boundary_conditions(&self.info, &mut self.primitives.vz);
            apply_boundary_conditions(&self.info, &mut self.primitives.p);
        }

        pub fn update_primitive_xtr_ghost_cells(&mut self) {
            handle_walls(
                &self.info,
                &mut self.primitives_xtr.rho,
                &mut self.primitives_xtr.vx,
                &mut self.primitives_xtr.vy,
                &mut self.primitives_xtr.vz,
                &mut self.primitives_xtr.p,
            );
            apply_boundary_conditions(&self.info, &mut self.primitives_xtr.rho);
            apply_boundary_conditions(&self.info, &mut self.primitives_xtr.vx);
            apply_boundary_conditions(&self.info, &mut self.primitives_xtr.vy);
            apply_boundary_conditions(&self.info, &mut self.primitives_xtr.vz);
            apply_boundary_conditions(&self.info, &mut self.primitives_xtr.p);
        }

        pub fn update_gradients_ghost_cells(&mut self) {
            apply_boundary_conditions4(&self.info, &mut self.gradients.drho);
            apply_boundary_conditions4(&self.info, &mut self.gradients.dvx);
            apply_boundary_conditions4(&self.info, &mut self.gradients.dvy);
            apply_boundary_conditions4(&self.info, &mut self.gradients.dvz);
            apply_boundary_conditions4(&self.info, &mut self.gradients.dp);
        }

        pub fn update_fluxes_ghost_cells(&mut self) {
            apply_boundary_conditions4(&self.info, &mut self.fluxes.mass);
            apply_boundary_conditions4(&self.info, &mut self.fluxes.momentum_x);
            apply_boundary_conditions4(&self.info, &mut self.fluxes.momentum_y);
            apply_boundary_conditions4(&self.info, &mut self.fluxes.momentum_z);
            apply_boundary_conditions4(&self.info, &mut self.fluxes.energy);
        }

        pub fn new(info: Simulation) -> Self {
            let nx = info.nx;
            let ny = info.ny;
            let nz = info.nz;

            let p = Primitives {
                rho: Array3::<f32>::zeros((nx, ny, nz)),
                vx: Array3::<f32>::zeros((nx, ny, nz)),
                vy: Array3::<f32>::zeros((nx, ny, nz)),
                vz: Array3::<f32>::zeros((nx, ny, nz)),
                p: Array3::<f32>::zeros((nx, ny, nz)),
            };

            let p_xtr = Primitives {
                rho: Array3::<f32>::zeros((nx, ny, nz)),
                vx: Array3::<f32>::zeros((nx, ny, nz)),
                vy: Array3::<f32>::zeros((nx, ny, nz)),
                vz: Array3::<f32>::zeros((nx, ny, nz)),
                p: Array3::<f32>::zeros((nx, ny, nz)),
            };

            let c = Conserved {
                mass: Array3::<f32>::zeros((nx, ny, nz)),
                momentum_x: Array3::<f32>::zeros((nx, ny, nz)),
                momentum_y: Array3::<f32>::zeros((nx, ny, nz)),
                momentum_z: Array3::<f32>::zeros((nx, ny, nz)),
                energy: Array3::<f32>::zeros((nx, ny, nz)),
            };

            let g = Gradients {
                drho: Array4::<f32>::zeros((nx, ny, nz, 3)),
                dvx: Array4::<f32>::zeros((nx, ny, nz, 3)),
                dvy: Array4::<f32>::zeros((nx, ny, nz, 3)),
                dvz: Array4::<f32>::zeros((nx, ny, nz, 3)),
                dp: Array4::<f32>::zeros((nx, ny, nz, 3)),
            };

            let f = Fluxes {
                mass: Array4::<f32>::zeros((nx, ny, nz, 3)),
                momentum_x: Array4::<f32>::zeros((nx, ny, nz, 3)),
                momentum_y: Array4::<f32>::zeros((nx, ny, nz, 3)),
                momentum_z: Array4::<f32>::zeros((nx, ny, nz, 3)),
                energy: Array4::<f32>::zeros((nx, ny, nz, 3)),
            };
            SimulationGrid {
                info: info,
                primitives: p,
                primitives_xtr: p_xtr,
                conserved: c,
                gradients: g,
                fluxes: f,
            }
        }
    }
}
