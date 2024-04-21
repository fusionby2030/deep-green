pub mod io {
    use crate::SimulationGrid;
    use hdf5::File;

    pub fn save_to_disk(grid: &SimulationGrid, fname: &str) {
        let shape = grid.primitives.rho.shape();
        let gradient_shape = grid.gradients.drho.shape();
        let fluxes_shape = grid.fluxes.mass.shape();
        let file = File::create(fname).unwrap();
        let group_primitives = file.create_group("primitives").unwrap();
        let group_conserved = file.create_group("conserved").unwrap();
        let group_gradients = file.create_group("gradients").unwrap();
        let group_fluxes = file.create_group("fluxes").unwrap();
        {
            let dataset = group_primitives
                .new_dataset::<f64>()
                .shape(shape)
                .create("rho")
                .unwrap();
            dataset.write(&grid.primitives.rho).unwrap();
        }

        {
            let dataset = group_primitives
                .new_dataset::<f64>()
                .shape(shape)
                .create("vx")
                .unwrap();
            dataset.write(&grid.primitives.vx).unwrap();
        }

        {
            let dataset = group_primitives
                .new_dataset::<f64>()
                .shape(shape)
                .create("vy")
                .unwrap();
            dataset.write(&grid.primitives.vy).unwrap();
        }

        {
            let dataset = group_primitives
                .new_dataset::<f64>()
                .shape(shape)
                .create("vz")
                .unwrap();
            dataset.write(&grid.primitives.vz).unwrap();
        }

        {
            let dataset = group_primitives
                .new_dataset::<f64>()
                .shape(shape)
                .create("pressure")
                .unwrap();
            dataset.write(&grid.primitives.p).unwrap();
        }

        //Conserved
        {
            let dataset = group_conserved
                .new_dataset::<f64>()
                .shape(shape)
                .create("mass")
                .unwrap();
            dataset.write(&grid.conserved.mass).unwrap();
        }

        {
            let dataset = group_conserved
                .new_dataset::<f64>()
                .shape(shape)
                .create("momentum_x")
                .unwrap();
            dataset.write(&grid.conserved.momentum_x).unwrap();
        }

        {
            let dataset = group_conserved
                .new_dataset::<f64>()
                .shape(shape)
                .create("momentum_y")
                .unwrap();
            dataset.write(&grid.conserved.momentum_y).unwrap();
        }

        {
            let dataset = group_conserved
                .new_dataset::<f64>()
                .shape(shape)
                .create("momentum_z")
                .unwrap();
            dataset.write(&grid.conserved.momentum_z).unwrap();
        }

        {
            let dataset = group_conserved
                .new_dataset::<f64>()
                .shape(shape)
                .create("energy")
                .unwrap();
            dataset.write(&grid.conserved.energy).unwrap();
        }

        //Gradients
        {
            let dataset = group_gradients
                .new_dataset::<f64>()
                .shape(gradient_shape)
                .create("rho")
                .unwrap();
            dataset.write(&grid.gradients.drho).unwrap();
        }

        {
            let dataset = group_gradients
                .new_dataset::<f64>()
                .shape(gradient_shape)
                .create("vx")
                .unwrap();
            dataset.write(&grid.gradients.dvx).unwrap();
        }

        {
            let dataset = group_gradients
                .new_dataset::<f64>()
                .shape(gradient_shape)
                .create("vy")
                .unwrap();
            dataset.write(&grid.gradients.dvy).unwrap();
        }

        {
            let dataset = group_gradients
                .new_dataset::<f64>()
                .shape(gradient_shape)
                .create("vz")
                .unwrap();
            dataset.write(&grid.gradients.dvz).unwrap();
        }

        {
            let dataset = group_gradients
                .new_dataset::<f64>()
                .shape(gradient_shape)
                .create("pressure")
                .unwrap();
            dataset.write(&grid.gradients.dp).unwrap();
        }

        //Fluxes
        {
            let dataset = group_fluxes
                .new_dataset::<f64>()
                .shape(fluxes_shape)
                .create("mass")
                .unwrap();
            dataset.write(&grid.fluxes.mass).unwrap();
        }

        {
            let dataset = group_fluxes
                .new_dataset::<f64>()
                .shape(fluxes_shape)
                .create("momentum_x")
                .unwrap();
            dataset.write(&grid.fluxes.momentum_x).unwrap();
        }

        {
            let dataset = group_fluxes
                .new_dataset::<f64>()
                .shape(fluxes_shape)
                .create("momentum_y")
                .unwrap();
            dataset.write(&grid.fluxes.momentum_y).unwrap();
        }

        {
            let dataset = group_fluxes
                .new_dataset::<f64>()
                .shape(fluxes_shape)
                .create("momentum_z")
                .unwrap();
            dataset.write(&grid.fluxes.momentum_z).unwrap();
        }

        {
            let dataset = group_fluxes
                .new_dataset::<f64>()
                .shape(fluxes_shape)
                .create("energy")
                .unwrap();
            dataset.write(&grid.fluxes.energy).unwrap();
        }
    }
}
