mod body;
mod quad_tree;
mod simulation;

use quad_tree::Boundary;
use simulation::{generate_random_bodies, Simulation};

fn main() {
    // Define the simulation boundary
    let boundary = Boundary {
        x_min: -1.0e6, // Minimum x-coordinate
        x_max: 1.0e6,  // Maximum x-coordinate
        y_min: -1.0e6, // Minimum y-coordinate
        y_max: 1.0e6,  // Maximum y-coordinate
    };

    // Generate random bodies for the simulation
    let num_bodies = 100; // Number of bodies
    let bodies = generate_random_bodies(num_bodies, &boundary);

    // Simulation parameters
    let theta = 0.5; // Barnes-Hut approximation threshold
    let time_step = 1.0; // Time step in seconds
    let num_steps = 10; // Number of simulation steps

    // Initialize the simulation
    let mut simulation = Simulation::new(bodies, boundary, theta, time_step);

    // Run the simulation
    println!("Starting simulation with {} bodies...", num_bodies);
    simulation.run(num_steps);
}
