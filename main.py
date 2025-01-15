import pygame
import numpy as np
import math
import random
import datetime
import os
import imageio

# Recording variables
is_recording = False
frame_count = 0
frames_dir = "frames"
os.makedirs(frames_dir, exist_ok=True)
frames_per_stage = 10

# Constants
WIDTH, HEIGHT = 800, 800
G = 1  # Gravitational constant (larger = stronger attraction)
THETA = 0.5  # Barnes-Hut parameter (larger = faster, but larger error)
DT = 0.1  # Time step (larger = slower, but more precise position incrementation)

# Pygame setup
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("2D N-Body Simulation with Barnes-Hut")
clock = pygame.time.Clock()


# Body class
class Body:
    def __init__(self, x, y, mass, vx=0, vy=0):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.radius = max(2, int(math.log(mass)))
        self.trail = []

    def draw(self, screen):
        pygame.draw.circle(screen, (0, 0, 0), (int(self.x), int(self.y)), self.radius)

    def draw_velocity(self, screen):
        arrow_length = math.sqrt(self.vx**2 + self.vy**2) * 5
        if arrow_length > 0:
            end_x = self.x + self.vx * 5
            end_y = self.y + self.vy * 5
            pygame.draw.line(screen, (0, 255, 0), (self.x, self.y), (end_x, end_y), 1)

    def draw_trail(self, screen):
        if len(self.trail) >= 2:
            for i in range(1, len(self.trail)):
                alpha = int(255 * (i / len(self.trail)))
                color = (100, 100, 100, alpha)
                start_pos = self.trail[i - 1]
                end_pos = self.trail[i]
                pygame.draw.line(screen, color, start_pos, end_pos, 1)

    def update_position(self):
        self.x += self.vx * DT
        self.y += self.vy * DT
        self.trail.append((int(self.x), int(self.y)))
        if len(self.trail) > 100:
            self.trail.pop(0)

    def apply_force(self, fx, fy):
        self.vx += fx / self.mass * DT
        self.vy += fy / self.mass * DT


# Barnes-Hut QuadTree Node
class QuadTreeNode:
    def __init__(self, x, y, width, height):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.bodies = []
        self.children = []
        self.total_mass = 0
        self.center_of_mass = (0, 0)

    def insert(self, body):
        if not self._contains(body):
            return False

        if len(self.bodies) < 1:
            self.bodies.append(body)
            self._update_mass(body)
            return True

        if not self.children:
            self._subdivide()

        for child in self.children:
            if child.insert(body):
                return True

        return False

    def _contains(self, body):
        return (
            self.x <= body.x < self.x + self.width
            and self.y <= body.y < self.y + self.height
        )

    def _subdivide(self):
        half_width = self.width / 2
        half_height = self.height / 2
        self.children = [
            QuadTreeNode(self.x, self.y, half_width, half_height),
            QuadTreeNode(self.x + half_width, self.y, half_width, half_height),
            QuadTreeNode(self.x, self.y + half_height, half_width, half_height),
            QuadTreeNode(
                self.x + half_width, self.y + half_height, half_width, half_height
            ),
        ]

    def _update_mass(self, body):
        self.total_mass += body.mass
        self.center_of_mass = (
            (
                self.center_of_mass[0] * (self.total_mass - body.mass)
                + body.x * body.mass
            )
            / self.total_mass,
            (
                self.center_of_mass[1] * (self.total_mass - body.mass)
                + body.y * body.mass
            )
            / self.total_mass,
        )

    def compute_force(self, body):
        if not self.bodies:
            return 0, 0

        dx = self.center_of_mass[0] - body.x
        dy = self.center_of_mass[1] - body.y
        distance_sq = dx * dx + dy * dy

        if distance_sq < 1e-10:
            return 0, 0

        if self.width * self.width / distance_sq < THETA * THETA or not self.children:
            force_magnitude = G * self.total_mass * body.mass / distance_sq
            force_x = force_magnitude * dx / math.sqrt(distance_sq)
            force_y = force_magnitude * dy / math.sqrt(distance_sq)
            return force_x, force_y

        force_x, force_y = 0, 0
        for child in self.children:
            fx, fy = child.compute_force(body)
            force_x += fx
            force_y += fy
        return force_x, force_y

    def draw(self, screen):
        pygame.draw.rect(
            screen, (102, 153, 255), (self.x, self.y, self.width, self.height), 1
        )
        for child in self.children:
            child.draw(screen)


# Collision handling
def handle_collisions(bodies, e=0.8):
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            body1 = bodies[i]
            body2 = bodies[j]

            dx = body2.x - body1.x
            dy = body2.y - body1.y
            distance = math.sqrt(dx * dx + dy * dy)

            min_distance = body1.radius + body2.radius
            if distance < min_distance:
                nx = dx / distance
                ny = dy / distance

                v_rel_x = body1.vx - body2.vx
                v_rel_y = body1.vy - body2.vy

                v_rel_dot_n = v_rel_x * nx + v_rel_y * ny

                J = (1 + e) * v_rel_dot_n / (1 / body1.mass + 1 / body2.mass)

                body1.vx -= J * nx / body1.mass
                body1.vy -= J * ny / body1.mass
                body2.vx += J * nx / body2.mass
                body2.vy += J * ny / body2.mass

                # Add minor buffer (sometimes bodies overlap...)
                overlap = min_distance - distance
                body1.x -= overlap * nx / 2
                body1.y -= overlap * ny / 2
                body2.x += overlap * nx / 2
                body2.y += overlap * ny / 2


def main():
    bodies = [
        Body(random.uniform(0, WIDTH), random.uniform(0, HEIGHT), random.uniform(1, 10))
        for _ in range(500)
    ]

    show_velocity = False
    show_quadtree = False
    show_trails = False

    global is_recording, frame_count
    is_recording = False
    frame_count = 0
    simulation_steps = 0
    running = True

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_v:
                    show_velocity = not show_velocity
                if event.key == pygame.K_q:
                    show_quadtree = not show_quadtree
                if event.key == pygame.K_s:
                    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"screenshot_{timestamp}.png"
                    pygame.image.save(screen, filename)
                    print(f"Screenshot saved as {filename}")
                if event.key == pygame.K_t:
                    show_trails = not show_trails
                if event.key == pygame.K_r:
                    is_recording = not is_recording
                    if is_recording:
                        print("Recording started...")
                    else:
                        print("Recording stopped. Compiling frames...")
                        compile_frames_to_gif()
                        frame_count = 0
        screen.fill((255, 255, 255))

        root = QuadTreeNode(0, 0, WIDTH, HEIGHT)
        for body in bodies:
            root.insert(body)

        if show_quadtree:
            root.draw(screen)

        if show_trails:
            for body in bodies:
                body.draw_trail(screen)

        for body in bodies:
            fx, fy = root.compute_force(body)
            body.apply_force(fx, fy)
            body.update_position()
            body.draw(screen)
            if show_velocity:
                body.draw_velocity(screen)

        # Handle collisions after every step
        handle_collisions(bodies)

        # Capture frame if recording and at the right stage
        if is_recording and simulation_steps % frames_per_stage == 0:
            frame_filename = os.path.join(frames_dir, f"frame_{frame_count:04d}.png")
            pygame.image.save(screen, frame_filename)
            frame_count += 1

        simulation_steps += 1

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


def compile_frames_to_gif():
    frame_files = sorted(
        [
            os.path.join(frames_dir, f)
            for f in os.listdir(frames_dir)
            if f.startswith("frame_")
        ]
    )
    if not frame_files:
        print("No frames to compile.")
        return

    output_file = "simulation.gif"
    with imageio.get_writer(output_file, mode="I", fps=30) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)
            writer.append_data(image)
        print(f"GIF saved as {output_file}")

    # Wipe frames after gif is created
    for frame_file in frame_files:
        os.remove(frame_file)
    print("Frames cleaned up.")


if __name__ == "__main__":
    main()
