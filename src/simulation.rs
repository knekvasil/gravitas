use crate::body::Body;
use crate::quad_tree::{Boundary, QuadTree};
use std::f64::consts::PI;

pub struct Simulation {
    pub bodies: Vec<Body>,
    pub quad_tree: QuadTree,
    pub theta: f64,     // Threshold for Barnes-Hut approximation
    pub time_step: f64, // Time step for the simulation
}

impl Simulation {
    pub fn new(bodies: Vec<Body>, boundary: Boundary, theta: f64, time_step: f64) -> Self {
        let quad_tree = QuadTree::new(boundary);
        Self {
            bodies,
            quad_tree,
            theta,
            time_step,
        }
    }

    pub fn update(&mut self) {
        // Clear and rebuild the quadtree
        self.quad_tree = QuadTree::new(self.quad_tree.get_boundary().clone());
        for body in &self.bodies {
            self.quad_tree.insert(body.clone());
        }

        // Reset accelerations
        for body in &mut self.bodies {
            body.acceleration = (0.0, 0.0);
        }

        // Calculate forces and update bodies
        for body in &mut self.bodies {
            let force = self.quad_tree.calculate_force(body, self.theta);
            body.apply_force(force);
        }

        // Update positions and velocities
        for body in &mut self.bodies {
            body.update_velocity(self.time_step);
            body.update_position(self.time_step);
        }
    }

    pub fn run(&mut self, steps: usize) {
        for step in 0..steps {
            println!("Step: {}", step + 1);
            self.update();
            self.print_summary();
        }
    }

    fn print_summary(&self) {
        println!("Simulation Summary:");
        for (i, body) in self.bodies.iter().enumerate() {
            println!(
                "Body {}: Position ({:.2}, {:.2}), Velocity ({:.2}, {:.2})",
                i, body.position.0, body.position.1, body.velocity.0, body.velocity.1
            );
        }
    }
}

// Example of creating a random set of bodies (can be replaced with actual data)
pub fn generate_random_bodies(count: usize, boundary: &Boundary) -> Vec<Body> {
    let mut bodies = Vec::new();
    for _ in 0..count {
        let position = (
            rand::random::<f64>() * (boundary.x_max - boundary.x_min) + boundary.x_min,
            rand::random::<f64>() * (boundary.y_max - boundary.y_min) + boundary.y_min,
        );
        let velocity = (0.0, 0.0);
        let mass = rand::random::<f64>() * 1e5 + 1e3; // Random mass between 1e3 and 1e5
        let _radius = (mass / (4.0 / 3.0 * PI)).cbrt(); // Approximate radius based on mass
        bodies.push(Body::new(position, velocity, mass, _radius));
    }
    bodies
}
