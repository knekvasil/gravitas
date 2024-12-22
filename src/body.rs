// body.rs
#[derive(Clone)]
pub struct Body {
    pub position: (f64, f64),
    pub velocity: (f64, f64),
    pub acceleration: (f64, f64),
    pub mass: f64,
    pub _radius: f64,
}

impl Body {
    pub fn new(position: (f64, f64), velocity: (f64, f64), mass: f64, _radius: f64) -> Self {
        Self {
            position,
            velocity,
            acceleration: (0.0, 0.0),
            mass,
            _radius: 1.0,
        }
    }

    pub fn update_position(&mut self, dt: f64) {
        self.position.0 += self.velocity.0 * dt + 0.5 * self.acceleration.0 * dt * dt;
        self.position.1 += self.velocity.1 * dt + 0.5 * self.acceleration.1 * dt * dt;
    }

    pub fn update_velocity(&mut self, dt: f64) {
        self.velocity.0 += self.acceleration.0 * dt;
        self.velocity.1 += self.acceleration.1 * dt;
    }

    pub fn apply_force(&mut self, force: (f64, f64)) {
        self.acceleration.0 = force.0 / self.mass;
        self.acceleration.1 = force.1 / self.mass;
    }
}
