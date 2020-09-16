use rand::distributions::{Distribution, Standard};
use rand::seq::SliceRandom;
use rand::{rngs::ThreadRng, thread_rng, Rng};
use std::collections::HashMap;
use std::fmt;
use std::io;
use std::io::Write;
use std::{thread, time};
use termion::raw::IntoRawMode;
use termion::screen::AlternateScreen;
use termion::{cursor, terminal_size};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    N,
    S,
    E,
    W,
}

impl Distribution<Direction> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Direction {
        match rng.gen_range(0, 4) {
            0 => Direction::N,
            1 => Direction::S,
            2 => Direction::E,
            3 => Direction::W,
            _ => unreachable!(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Cell {
    x: usize,
    y: usize,
}
impl Cell {
    /// Create a neighbor cell by moving in the direction given.
    fn neighbor(&self, direction: &Direction, width: usize, height: usize) -> Option<Cell> {
        match direction {
            Direction::N => {
                if self.y > 0 {
                    Some(Cell {
                        x: self.x,
                        y: self.y - 1,
                    })
                } else {
                    None
                }
            }
            Direction::S => {
                if self.y + 1 < height {
                    Some(Cell {
                        x: self.x,
                        y: self.y + 1,
                    })
                } else {
                    None
                }
            }
            Direction::E => {
                if self.x + 1 < width {
                    Some(Cell {
                        x: self.x + 1,
                        y: self.y,
                    })
                } else {
                    None
                }
            }
            Direction::W => {
                if self.x > 0 {
                    Some(Cell {
                        x: self.x - 1,
                        y: self.y,
                    })
                } else {
                    None
                }
            }
        }
    }
}

/// Generates random cells within a grid bounds.
struct RandomCellGen {
    rng: ThreadRng,
    width: usize,
    height: usize,
}
impl RandomCellGen {
    fn new(width: usize, height: usize) -> RandomCellGen {
        RandomCellGen {
            rng: thread_rng(),
            width,
            height,
        }
    }

    /// Generate a random cell within the grid bounds.
    fn cell(&mut self) -> Cell {
        Cell {
            x: self.rng.gen_range(0, self.width),
            y: self.rng.gen_range(0, self.height),
        }
    }

    /// Generate a random cell next to the given cell, within the grid bounds,
    /// along with the direction we travelled to get there.
    fn neighbor(&mut self, cell: Cell) -> (Cell, Direction) {
        let mut all_dirs: [Direction; 4] = [Direction::N, Direction::S, Direction::E, Direction::W];
        all_dirs.shuffle(&mut self.rng);

        for d in all_dirs.iter() {
            if let Some(neighbor) = cell.neighbor(d, self.width, self.height) {
                return (neighbor, *d);
            }
        }

        unreachable!("Should always be possible to generate a neighbor cell unless it's a 1x1 grid")
    }
}

/// A maze!
///
/// The inner Vectors are rows, to make rendering easier. It looks like:
///
/// `vec[row][column]`
///
/// For, example a 2x2 grid:
///
/// ```text
/// 0: [0,1]
/// 1: [0,1]
/// ```
///
/// A cell `{x: 0, y: 1}` would access a vector like `vec[y][x]`, or `vec[1][0]`.
#[derive(Debug)]
struct Maze {
    h_walls: Vec<Vec<bool>>,
    v_walls: Vec<Vec<bool>>,
}

impl Maze {
    /// Use [Wilson's algorithm](https://en.wikipedia.org/wiki/Maze_generation_algorithm#Wilson's_algorithm) to generate a maze.
    ///
    /// # Arguments
    ///
    /// * `width` - number of cells wide
    /// * `height` - number of cells high
    ///
    /// The actual width/height in characters will be greater.
    fn new(width: usize, height: usize) -> Maze {
        let mut rand_cell = RandomCellGen::new(width, height);

        let mut h_walls = Maze::create_grid(width, height + 1, true);
        let mut v_walls = Maze::create_grid(width + 1, height, true);
        let mut cells_visited = Maze::create_grid(width, height, false);
        let mut num_cells_visited = 0;

        {
            let initial = rand_cell.cell();
            cells_visited[initial.y][initial.x] = true;
            num_cells_visited += 1;
        }

        let mut walk: HashMap<Cell, Direction> = HashMap::new();

        while num_cells_visited < width * height {
            // choose a starting point which hasn't been visited
            let mut init = rand_cell.cell();
            while cells_visited[init.y][init.x] {
                init = rand_cell.cell();
            }
            // do a random walk until we meet an already visited cell
            let mut prev = init;
            'walk: loop {
                let (cell, direction) = rand_cell.neighbor(prev);
                walk.insert(prev, direction);
                prev = cell;

                if cells_visited[cell.y][cell.x] {
                    // we've joined up with a visited cell, walk the walk and add to the maze
                    let mut c = init;
                    while let Some(&d) = walk.get(&c) {
                        cells_visited[c.y][c.x] = true;
                        num_cells_visited += 1;

                        // break down some walls
                        if d == Direction::N || d == Direction::S {
                            h_walls[c.y + if d == Direction::S { 1 } else { 0 }][c.x] = false;
                        } else {
                            v_walls[c.y][c.x + if d == Direction::E { 1 } else { 0 }] = false;
                        }
                        c = c.neighbor(&d, width, height).unwrap();
                    }

                    walk.clear();

                    break 'walk;
                }
            }
        }

        Maze { v_walls, h_walls }
    }

    fn create_grid(width: usize, height: usize, value: bool) -> Vec<Vec<bool>> {
        std::iter::repeat(std::iter::repeat(value).take(width).collect())
            .take(height)
            .collect()
    }
}

impl fmt::Display for Maze {
    /// # Example output
    ///
    /// ```text
    /// +--+--+
    /// |  |  |
    /// +--+--+
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // FIXME: the ol' switcheroo
        let mut h = self.h_walls.iter();
        let mut v = self.v_walls.iter();

        while let Some(h_row) = h.next() {
            for &wall in h_row {
                write!(f, "{}", if wall { "+--" } else { "+  " })?;
            }
            write!(f, "+\n\r")?;
            if let Some(v_row) = v.next() {
                for &wall in v_row.iter().take(v_row.len() - 1) {
                    write!(f, "{}", if wall { "|  " } else { "   " })?;
                }
                write!(f, "|")?;
            }
            write!(f, "\n\r")?;
        }

        fmt::Result::Ok(())
    }
}

fn main() {
    // Get and lock the stdios so we don't have to get the lock all the time
    let stdout = io::stdout();
    // let stdout = stdout.lock();
    // let stdin = io::stdin();
    // let stdin = stdin.lock();
    // let stderr = io::stderr();
    // let mut stderr = stderr.lock();

    let mut stdout = stdout.into_raw_mode().unwrap();

    let maze = Maze::new(2, 2);
    write!(stdout, "{}", maze).unwrap();

    // {
    //     let mut screen = AlternateScreen::from(stdout);
    //     write!(screen, "{}", cursor::Goto(1, 1)).unwrap();
    //     write!(screen, "Welcome to the alternate screen.\n\r\nPlease wait patiently until we arrive back at the main screen in a about three seconds.\n\r").unwrap();
    //     screen.flush().unwrap();

    //     let size = terminal_size().unwrap();

    //     // let maze = Maze::new(size.0.into(), size.1.into());
    //     let maze = Maze::new(16, 16);
    //     write!(screen, "{}", maze).unwrap();

    //     thread::sleep(time::Duration::from_secs(10));
    // }
}