use std::f64::consts::PI;

use nalgebra::Vector3;

pub struct Grid {
    grid_size: usize,
    cell_size: f64,
    pub cells: Vec<Vec<Vector3<f64>>>,
}

#[must_use]
pub fn lower_bound_of_distance(n: usize) -> f64 {
    (4.0 - 1.0 / ((PI * n as f64) / (6 * (n - 2)) as f64).powi(2).sin()).sqrt()
}

impl Grid {
    #[must_use]
    pub fn new_from_vec(points: &[Vector3<f64>]) -> Self {
        let lower_bound_of_distance = lower_bound_of_distance(points.len());
        let grid_size = (2.0 / lower_bound_of_distance).floor() as usize;
        let cell_size = 2.0 / (grid_size as f64);
        let mut grid = Self {
            grid_size,
            cell_size,
            cells: vec![Vec::new(); grid_size * grid_size],
        };
        points.iter().for_each(|point| grid.push(*point));
        grid
    }

    #[must_use]
    pub fn to_vec(&self) -> Vec<Vector3<f64>> {
        let mut vec = Vec::new();
        self.cells.iter().for_each(|cell| {
            vec.extend(cell);
        });
        vec
    }

    fn empty_from_grid(grid: &Self) -> Self {
        Self {
            grid_size: grid.grid_size,
            cell_size: grid.cell_size,
            cells: vec![Vec::new(); grid.grid_size * grid.grid_size],
        }
    }

    fn push(&mut self, point: Vector3<f64>) {
        let index = self.index_from(self.hash(&point));
        self.cells[index].push(point);
    }

    const fn index_from(&self, (hash_x, hash_y): (usize, usize)) -> usize {
        hash_x * self.grid_size + hash_y
    }

    fn hash(&self, point: &Vector3<f64>) -> (usize, usize) {
        let hash_x = ((point.x + 1.0) / self.cell_size) as usize;
        let hash_y = ((point.y + 1.0) / self.cell_size) as usize;
        (hash_x, hash_y)
    }

    fn find_closest_non_overlap(&self, current: &Vector3<f64>) -> Option<&Vector3<f64>> {
        let (hash_x, hash_y) = self.hash(current);
        let x_from = hash_x.saturating_sub(1);
        let x_to = (hash_x + 1).min(self.grid_size - 1);
        let y_from = hash_y.saturating_sub(1);
        let y_to = (hash_y + 1).min(self.grid_size - 1);

        let (mut closest, mut max_cos_theta) = (None, -1.0);
        (x_from..=x_to).for_each(|x| {
            (y_from..=y_to).for_each(|y| {
                self.cells[self.index_from((x, y))]
                    .iter()
                    .for_each(|point| {
                        let cos_theta = current.dot(point);
                        if cos_theta > max_cos_theta && point != current {
                            (closest, max_cos_theta) = (Some(point), cos_theta);
                        }
                    });
            });
        });
        closest
    }

    #[must_use]
    pub fn iterate(&self, l: f64) -> Self {
        let mut grid = Self::empty_from_grid(self);
        self.cells.iter().for_each(|cell| {
            for current in cell {
                grid.push(
                    self.find_closest_non_overlap(current)
                        .map_or(*current, |closest| {
                            let direction = (current - closest).normalize();
                            (current + direction * l).normalize()
                        }),
                );
            }
        });
        grid
    }
}
