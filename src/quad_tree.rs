use crate::body::Body;

#[derive(Clone, Copy)]
pub struct Boundary {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

impl Boundary {
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.x_min && x <= self.x_max && y >= self.y_min && y <= self.y_max
    }

    pub fn center(&self) -> (f64, f64) {
        (
            (self.x_min + self.x_max) / 2.0,
            (self.y_min + self.y_max) / 2.0,
        )
    }

    pub fn subdivide(&self) -> [Boundary; 4] {
        let (x_center, y_center) = self.center();
        [
            Boundary {
                x_min: self.x_min,
                x_max: x_center,
                y_min: self.y_min,
                y_max: y_center,
            },
            Boundary {
                x_min: x_center,
                x_max: self.x_max,
                y_min: self.y_min,
                y_max: y_center,
            },
            Boundary {
                x_min: self.x_min,
                x_max: x_center,
                y_min: y_center,
                y_max: self.y_max,
            },
            Boundary {
                x_min: x_center,
                x_max: self.x_max,
                y_min: y_center,
                y_max: self.y_max,
            },
        ]
    }
}

pub enum QuadTreeNode {
    Empty,
    Leaf {
        body: Body,
    },
    Internal {
        center_of_mass: (f64, f64),
        total_mass: f64,
        children: [Option<Box<QuadTreeNode>>; 4],
    },
}

impl QuadTreeNode {
    pub fn insert(&mut self, body: Body, boundary: Boundary) {
        match self {
            QuadTreeNode::Empty => {
                *self = QuadTreeNode::Leaf { body };
            }
            QuadTreeNode::Leaf {
                body: existing_body,
            } => {
                if existing_body.position == body.position {
                    // Avoid infinite recursion by ignoring identical positions
                    return;
                }

                let mut children = [None, None, None, None];
                let quadrants = boundary.subdivide();

                // Insert existing body
                for (i, quadrant) in quadrants.iter().enumerate() {
                    if quadrant.contains(existing_body.position.0, existing_body.position.1) {
                        children[i] = Some(Box::new(QuadTreeNode::Leaf {
                            body: existing_body.clone(),
                        }));
                        break;
                    }
                }

                // Insert new body
                for (i, quadrant) in quadrants.iter().enumerate() {
                    if quadrant.contains(body.position.0, body.position.1) {
                        if children[i].is_none() {
                            children[i] = Some(Box::new(QuadTreeNode::Leaf { body: body.clone() }));
                        } else {
                            children[i]
                                .as_mut()
                                .unwrap()
                                .insert(body.clone(), *quadrant);
                        }
                        break;
                    }
                }

                *self = QuadTreeNode::Internal {
                    center_of_mass: Self::calculate_center_of_mass(&children),
                    total_mass: body.mass + existing_body.mass,
                    children,
                };
            }
            QuadTreeNode::Internal {
                ref mut children,
                ref mut center_of_mass,
                ref mut total_mass,
                ..
            } => {
                let quadrants = boundary.subdivide();
                for (i, quadrant) in quadrants.iter().enumerate() {
                    if quadrant.contains(body.position.0, body.position.1) {
                        if children[i].is_none() {
                            children[i] = Some(Box::new(QuadTreeNode::Empty));
                        }
                        children[i]
                            .as_mut()
                            .unwrap()
                            .insert(body.clone(), *quadrant);
                        break;
                    }
                }
                *center_of_mass = Self::calculate_center_of_mass(children);
                *total_mass += body.mass;
            }
        }
    }

    fn calculate_center_of_mass(children: &[Option<Box<QuadTreeNode>>; 4]) -> (f64, f64) {
        let mut total_mass = 0.0;
        let mut x_mass_sum = 0.0;
        let mut y_mass_sum = 0.0;

        for child in children.iter() {
            if let Some(node) = child {
                match **node {
                    QuadTreeNode::Empty => continue,
                    QuadTreeNode::Leaf { ref body } => {
                        total_mass += body.mass;
                        x_mass_sum += body.position.0 * body.mass;
                        y_mass_sum += body.position.1 * body.mass;
                    }
                    QuadTreeNode::Internal {
                        center_of_mass,
                        total_mass: mass,
                        ..
                    } => {
                        total_mass += mass;
                        x_mass_sum += center_of_mass.0 * mass;
                        y_mass_sum += center_of_mass.1 * mass;
                    }
                }
            }
        }

        if total_mass > 0.0 {
            (x_mass_sum / total_mass, y_mass_sum / total_mass)
        } else {
            (0.0, 0.0)
        }
    }

    pub fn calculate_force(&self, body: &Body, theta: f64) -> (f64, f64) {
        match self {
            QuadTreeNode::Empty => (0.0, 0.0),
            QuadTreeNode::Leaf { body: other_body } => {
                if body.position == other_body.position {
                    (0.0, 0.0)
                } else {
                    calculate_gravity(body, other_body.position, other_body.mass)
                }
            }
            QuadTreeNode::Internal {
                center_of_mass,
                total_mass,
                children,
            } => {
                let dx = center_of_mass.0 - body.position.0;
                let dy = center_of_mass.1 - body.position.1;
                let d = (dx * dx + dy * dy).sqrt();

                if d == 0.0 || (self.get_boundary_size() / d) < theta {
                    calculate_gravity(body, *center_of_mass, *total_mass)
                } else {
                    let mut force = (0.0, 0.0);
                    for child in children.iter() {
                        if let Some(child) = child {
                            let child_force = child.calculate_force(body, theta);
                            force.0 += child_force.0;
                            force.1 += child_force.1;
                        }
                    }
                    force
                }
            }
        }
    }

    fn get_boundary_size(&self) -> f64 {
        match self {
            QuadTreeNode::Internal { children, .. } => {
                if let Some(child) = &children[0] {
                    let boundary = child.get_boundary();
                    boundary.x_max - boundary.x_min
                } else {
                    0.0
                }
            }
            _ => 0.0,
        }
    }

    fn get_boundary(&self) -> Boundary {
        match self {
            QuadTreeNode::Internal { children, .. } => {
                if let Some(child) = &children[0] {
                    child.get_boundary()
                } else {
                    Boundary {
                        x_min: 0.0,
                        x_max: 0.0,
                        y_min: 0.0,
                        y_max: 0.0,
                    }
                }
            }
            _ => Boundary {
                x_min: 0.0,
                x_max: 0.0,
                y_min: 0.0,
                y_max: 0.0,
            },
        }
    }
}

pub struct QuadTree {
    root: QuadTreeNode,
    boundary: Boundary,
}

impl QuadTree {
    pub fn new(boundary: Boundary) -> Self {
        Self {
            root: QuadTreeNode::Empty,
            boundary,
        }
    }

    pub fn insert(&mut self, body: Body) {
        self.root.insert(body, self.boundary);
    }

    pub fn calculate_force(&self, body: &Body, theta: f64) -> (f64, f64) {
        self.root.calculate_force(body, theta)
    }

    pub fn get_boundary(&self) -> &Boundary {
        &self.boundary
    }
}

fn calculate_gravity(body: &Body, other_pos: (f64, f64), other_mass: f64) -> (f64, f64) {
    const G: f64 = 6.67430e-11;
    let dx = other_pos.0 - body.position.0;
    let dy = other_pos.1 - body.position.1;
    let d_squared = dx * dx + dy * dy;
    let d = d_squared.sqrt();

    if d < 1e-10 {
        return (0.0, 0.0);
    }

    let force = G * body.mass * other_mass / d_squared;
    (force * dx / d, force * dy / d)
}
