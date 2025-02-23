use nalgebra::Vector3;
use rand::{self, Rng};
use tammes::{lower_bound_of_distance, Grid};

pub fn generate_points_on_sphere() -> Vector3<f64> {
    let mut rng = rand::rng();
    let z: f64 = rng.random_range(-1.0..1.0);
    let theta: f64 = rng.random_range(0.0..std::f64::consts::TAU);
    let r = (1.0 - z * z).sqrt();
    let (sin_theta, cos_theta) = theta.sin_cos();
    Vector3::new(r * cos_theta, r * sin_theta, z)
}

fn main() {
    const POINTS_NUM: usize = 2000;
    const ITERATIONS: usize = 10000;

    let points = (0..POINTS_NUM)
        .map(|_| generate_points_on_sphere())
        .collect::<Vec<_>>();

    let lower_bound_of_distance = lower_bound_of_distance(POINTS_NUM);

    let mut grid = Grid::new_from_vec(&points);
    for i in 0..ITERATIONS {
        let l = (lower_bound_of_distance) * (ITERATIONS - i) as f64 * (1.0 / (ITERATIONS as f64));
        grid = grid.iterate(l);

    }

    for point in grid.to_vec() {
        println!("{{{},{},{}}},", point.x, point.y, point.z);
    }
}
