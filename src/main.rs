use rand::distributions::{Distribution, Standard};
use rand::seq::SliceRandom;
use rand::{rngs::ThreadRng, thread_rng, Rng};
use std::collections::HashMap;
use std::convert::TryInto;
use std::io;
use std::io::Write;
use std::ops::Range;
use termion::event::Key;
use termion::input::TermRead;
use termion::raw::IntoRawMode;
use termion::screen::AlternateScreen;
use termion::{cursor, terminal_size};

#[derive(Debug, Clone)]
struct Size {
    width: usize,
    height: usize,
}
#[derive(Debug, Clone)]
struct Position {
    x: isize,
    y: isize,
}
impl Position {
    fn up(&self) -> Position {
        let mut p = self.clone();
        p.y -= 1;
        p
    }
    fn down(&self) -> Position {
        let mut p = self.clone();
        p.y += 1;
        p
    }
    fn left(&self) -> Position {
        let mut p = self.clone();
        p.x -= 1;
        p
    }
    fn right(&self) -> Position {
        let mut p = self.clone();
        p.x += 1;
        p
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Area {
    width: Range<usize>,
    height: Range<usize>,
}

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
struct CellPosition {
    x: usize,
    y: usize,
}
impl CellPosition {
    /// Create a neighbor cell by moving in the direction given.
    fn neighbor(&self, direction: &Direction, width: usize, height: usize) -> Option<CellPosition> {
        match direction {
            Direction::N => {
                if self.y > 0 {
                    Some(CellPosition {
                        x: self.x,
                        y: self.y - 1,
                    })
                } else {
                    None
                }
            }
            Direction::S => {
                if self.y + 1 < height {
                    Some(CellPosition {
                        x: self.x,
                        y: self.y + 1,
                    })
                } else {
                    None
                }
            }
            Direction::E => {
                if self.x + 1 < width {
                    Some(CellPosition {
                        x: self.x + 1,
                        y: self.y,
                    })
                } else {
                    None
                }
            }
            Direction::W => {
                if self.x > 0 {
                    Some(CellPosition {
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
    fn cell(&mut self) -> CellPosition {
        CellPosition {
            x: self.rng.gen_range(0, self.width),
            y: self.rng.gen_range(0, self.height),
        }
    }

    /// Generate a random cell next to the given cell, within the grid bounds,
    /// along with the direction we travelled to get there.
    fn neighbor(&mut self, cell: CellPosition) -> (CellPosition, Direction) {
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

enum MazeSizeUnit {
    Chars,
    #[allow(dead_code)]
    Cells,
}
impl MazeSizeUnit {
    fn to_chars(&self, (width, height): (usize, usize)) -> (usize, usize) {
        match self {
            MazeSizeUnit::Cells => MazeSizeUnit::cell_to_char((width, height)),
            MazeSizeUnit::Chars => {
                // We might ask for a maze of size 5x5 chars, but we can only fit 1x2 cells
                // which is 4x5 chars
                MazeSizeUnit::cell_to_char(MazeSizeUnit::char_to_cell((width, height)))
            }
        }
    }

    fn to_cells(&self, (width, height): (usize, usize)) -> (usize, usize) {
        match self {
            MazeSizeUnit::Cells => (width, height),
            MazeSizeUnit::Chars => MazeSizeUnit::char_to_cell((width, height)),
        }
    }

    fn char_to_cell((width, height): (usize, usize)) -> (usize, usize) {
        ((width - 1) / 3, (height - 1) / 2)
    }

    fn cell_to_char((width, height): (usize, usize)) -> (usize, usize) {
        ((width * 3) + 1, (height * 2) + 1)
    }
}

/// A maze!
///
/// # Example output
///
/// ```text
/// +--+--+
/// |  |  |
/// +--+--+
/// ```
#[derive(Debug)]
struct Maze {
    buffer: Vec<char>,
    char_width: usize,
    char_height: usize,
}

impl Maze {
    /// Use [Wilson's algorithm](https://en.wikipedia.org/wiki/Maze_generation_algorithm#Wilson's_algorithm) to generate a maze.
    fn new(size: Size, unit: MazeSizeUnit) -> Maze {
        let (width, height): (usize, usize) = (size.width.into(), size.height.into());
        let (char_width, char_height) = unit.to_chars((width, height));
        let (cells_wide, cells_high) = unit.to_cells((width, height));
        let mut rand_cell = RandomCellGen::new(cells_wide, cells_high);
        //
        // The inner Vectors are rows, to make rendering easier. It looks like:
        //
        // `vec[row][column]`
        //
        // For, example a 2x2 grid:
        //
        // ```text
        // 0: [0,1]
        // 1: [0,1]
        // ```
        //
        // A cell `{x: 0, y: 1}` would access a vector like `vec[y][x]`, or `vec[1][0]`.
        //
        let mut h_walls = Maze::create_grid(cells_wide, cells_high + 1, true);
        let mut v_walls = Maze::create_grid(cells_wide + 1, cells_high, true);
        let mut cells_visited = Maze::create_grid(cells_wide, cells_high, false);
        let mut num_cells_visited = 0;

        {
            let initial = rand_cell.cell();
            cells_visited[initial.y][initial.x] = true;
            num_cells_visited += 1;
        }

        let mut walk: HashMap<CellPosition, Direction> = HashMap::new();

        while num_cells_visited < cells_wide * cells_high {
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
                        c = c.neighbor(&d, cells_wide, cells_high).unwrap();
                    }

                    walk.clear();

                    break 'walk;
                }
            }
        }

        // entrace and exit
        let row = cells_high / 2;
        v_walls[row][0] = false;
        v_walls[row][cells_wide] = false;

        let mut buffer = Vec::with_capacity(char_width * char_height);
        Maze::draw_to_buffer(&h_walls, &v_walls, &mut buffer);

        Maze {
            buffer,
            char_width,
            char_height,
        }
    }

    fn create_grid(width: usize, height: usize, value: bool) -> Vec<Vec<bool>> {
        std::iter::repeat(std::iter::repeat(value).take(width).collect())
            .take(height)
            .collect()
    }

    /// Draw the maze into `w`.
    /// `term_size` is the size of the terminal represented by `w`.
    /// `pos` is the position of the top-left corner of the maze within the terminal.
    fn draw<W: Write>(
        &self,
        w: &mut W,
        term_size: &Size,
        pos: &Position,
    ) -> std::result::Result<(), std::io::Error> {
        // Goto is 1-based
        let term_x = (pos.x + 1).max(1).try_into().unwrap();
        let mut term_y = (pos.y + 1).max(1).try_into().unwrap();
        write!(w, "{}", cursor::Goto(term_x, term_y))?;
        if let Some(area) = self.on_screen_area(term_size, pos) {
            for y in area.height {
                for x in area.width.clone() {
                    write!(
                        w,
                        "{}", //
                        self.buffer[(y * self.char_width) + x]
                    )?;
                }
                term_y += 1;
                write!(w, "{}", cursor::Goto(term_x, term_y))?;
            }
        }

        Ok(())
    }

    /// Determines which area of the maze is visible on screen when drawing to a
    /// terminal of size `term_size`, with the top-left of the maze at `pos`.
    fn on_screen_area(&self, term_size: &Size, pos: &Position) -> Option<Area> {
        // convert maze position to terminal coordinates
        let offset_maze_x_s = pos.x;
        let offset_maze_x_e = pos.x + self.char_width as isize;
        let offset_maze_y_s = pos.y;
        let offset_maze_y_e = pos.y + self.char_height as isize;

        if offset_maze_x_s > term_size.width as isize || 0 > offset_maze_x_e {
            return None;
        }
        if offset_maze_y_s > term_size.height as isize || 0 > offset_maze_y_e {
            return None;
        }
        let x_s = offset_maze_x_s.max(0);
        let y_s = offset_maze_y_s.max(0);
        let x_e = (term_size.width as isize).min(offset_maze_x_e);
        let y_e = (term_size.height as isize).min(offset_maze_y_e);

        // convert back to maze coordinates
        Some(Area {
            width: (x_s - pos.x) as usize..(x_e - pos.x) as usize,
            height: (y_s - pos.y) as usize..(y_e - pos.y) as usize,
        })
    }

    fn draw_to_buffer(h_walls: &Vec<Vec<bool>>, v_walls: &Vec<Vec<bool>>, buf: &mut Vec<char>) {
        let mut h = h_walls.iter();
        let mut v = v_walls.iter();

        while let Some(h_row) = h.next() {
            for &wall in h_row {
                if wall {
                    buf.extend_from_slice(&['+', '-', '-']);
                } else {
                    buf.extend_from_slice(&['+', ' ', ' ']);
                }
            }
            buf.push('+');

            if let Some(v_row) = v.next() {
                if let Some((&last, walls)) = v_row.split_last() {
                    for &wall in walls {
                        if wall {
                            buf.extend_from_slice(&['|', ' ', ' ']);
                        } else {
                            buf.extend_from_slice(&[' ', ' ', ' ']);
                        }
                    }
                    if last {
                        buf.push('|');
                    } else {
                        buf.push(' ');
                    }
                }
            }
        }
    }

    /// Given the maze is at `maze_pos`, is there a wall at `pos`?
    fn is_wall(&self, maze_pos: &Position, pos: &Position) -> bool {
        let x = pos.x - maze_pos.x;
        let y = pos.y - maze_pos.y;

        if x < 0 || y < 0 {
            return false;
        }
        if x >= self.char_width as isize || y >= self.char_height as isize {
            return false;
        }

        self.buffer[((y as usize * self.char_width) + x as usize) as usize] != ' '
    }
}

fn draw_o<W: Write>(w: &mut W, pos: &Position) -> std::result::Result<(), std::io::Error> {
    write!(
        w,
        "{}{}",
        cursor::Goto(
            (pos.x + 1).try_into().unwrap(),
            (pos.y + 1).try_into().unwrap()
        ),
        "o"
    )?;
    Ok(())
}

fn main() {
    // Get and lock the stdios so we don't have to get the lock all the time
    let stdout = io::stdout();
    let stdout = stdout.lock();
    let stdin = io::stdin();
    let stdin = stdin.lock();

    let stdout = stdout.into_raw_mode().unwrap();

    {
        let screen = AlternateScreen::from(stdout);
        let mut screen = cursor::HideCursor::from(screen);

        let size = terminal_size().unwrap();
        let size = Size {
            width: size.0.into(),
            height: size.1.into(),
        };

        let o_pos = Position {
            x: size.width as isize / 2,
            y: size.height as isize / 2,
        };
        let mut maze_pos = Position {
            x: o_pos.x + 1,
            y: 0,
        };

        let maze = Maze::new(
            Size {
                width: size.width,
                height: size.height,
            },
            MazeSizeUnit::Chars,
        );

        maze.draw(&mut screen, &size, &maze_pos).unwrap();
        draw_o(&mut screen, &o_pos).unwrap();
        screen.flush().unwrap();

        for c in stdin.keys() {
            write!(screen, "{}", termion::clear::All).unwrap();
            match c {
                Ok(Key::Char('q')) => break,
                Ok(Key::Left) => {
                    let new_maze_pos = maze_pos.right();
                    if !maze.is_wall(&new_maze_pos, &o_pos) {
                        maze_pos = new_maze_pos;
                    }
                }
                Ok(Key::Right) => {
                    let new_maze_pos = maze_pos.left();
                    if !maze.is_wall(&new_maze_pos, &o_pos) {
                        maze_pos = new_maze_pos;
                    }
                }
                Ok(Key::Up) => {
                    let new_maze_pos = maze_pos.down();
                    if !maze.is_wall(&new_maze_pos, &o_pos) {
                        maze_pos = new_maze_pos;
                    }
                }
                Ok(Key::Down) => {
                    let new_maze_pos = maze_pos.up();
                    if !maze.is_wall(&new_maze_pos, &o_pos) {
                        maze_pos = new_maze_pos;
                    }
                }
                _ => {}
            }
            maze.draw(&mut screen, &size, &maze_pos).unwrap();
            draw_o(&mut screen, &o_pos).unwrap();
            screen.flush().unwrap();
        }

        write!(screen, "{}", termion::style::Reset).unwrap();
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn cell_to_char_1x1() {
        let size = MazeSizeUnit::cell_to_char((1, 1));
        assert_eq!(size, (4, 3));
    }

    #[test]
    fn cell_to_char_2x2() {
        let size = MazeSizeUnit::cell_to_char((2, 2));
        assert_eq!(size, (7, 5));
    }

    #[test]
    fn char_to_cell_4x3() {
        let size = MazeSizeUnit::char_to_cell((4, 3));
        assert_eq!(size, (1, 1));
    }
    #[test]
    fn char_to_cell_5x4() {
        let size = MazeSizeUnit::char_to_cell((5, 4));
        assert_eq!(size, (1, 1));
    }

    #[test]
    fn char_to_cell_6x4() {
        let size = MazeSizeUnit::char_to_cell((6, 4));
        assert_eq!(size, (1, 1));
    }

    #[test]
    fn char_to_cell_7x5() {
        let size = MazeSizeUnit::char_to_cell((7, 5));
        assert_eq!(size, (2, 2));
    }

    #[test]
    fn maze_size_chars_to_chars() {
        let size = MazeSizeUnit::Chars.to_chars((5, 5));
        assert_eq!(size, (4, 5));
    }

    #[test]
    fn on_screen_area_fits_inside() {
        let maze_size = Size {
            width: 5,
            height: 5,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: 0, y: 0 };

        let maze = Maze::new(maze_size.clone(), MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(
            area,
            Some(Area {
                width: 0..4,
                height: 0..5
            })
        );
    }

    #[test]
    fn on_screen_same_size() {
        let maze_size = Size {
            width: 10,
            height: 10,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: 0, y: 0 };

        let maze = Maze::new(maze_size.clone(), MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(
            area,
            Some(Area {
                width: 0..10,
                height: 0..9
            })
        );
    }

    #[test]
    fn on_screen_area_off_screen_to_right() {
        let maze_size = Size {
            width: 7,
            height: 7,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: 5, y: 0 };

        let maze = Maze::new(maze_size, MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(
            area,
            Some(Area {
                width: 0..5,
                height: 0..7
            })
        );
    }

    #[test]
    fn on_screen_area_completely_off_screen_to_right() {
        let maze_size = Size {
            width: 7,
            height: 7,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: 15, y: 0 };

        let maze = Maze::new(maze_size, MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(area, None);
    }

    #[test]
    fn on_screen_area_off_screen_to_left() {
        let maze_size = Size {
            width: 7,
            height: 7,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: -3, y: 0 };

        let maze = Maze::new(maze_size, MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(
            area,
            Some(Area {
                width: 3..7,
                height: 0..7
            })
        );
    }

    #[test]
    fn on_screen_area_completely_off_screen_to_left() {
        let maze_size = Size {
            width: 7,
            height: 7,
        };
        let term_size = Size {
            width: 10,
            height: 10,
        };
        let pos = Position { x: -15, y: 0 };

        let maze = Maze::new(maze_size, MazeSizeUnit::Chars);
        let area = maze.on_screen_area(&term_size, &pos);

        assert_eq!(area, None);
    }

    #[test]
    fn maze_buffer_1x1() {
        let mut buf = Vec::<char>::new();
        Maze::draw_to_buffer(
            &vec![vec![true], vec![true]],
            &vec![vec![true, true]],
            &mut buf,
        );
        assert_eq!(
            buf,
            vec![
                '+', '-', '-', '+', //
                '|', ' ', ' ', '|', //
                '+', '-', '-', '+',
            ]
        );
    }

    #[test]
    fn maze_buffer_2x2() {
        let mut buf = Vec::<char>::new();
        Maze::draw_to_buffer(
            &vec![vec![true, true], vec![true, true], vec![true, true]],
            &vec![vec![true, true, true], vec![true, true, true]],
            &mut buf,
        );
        assert_eq!(
            buf,
            vec![
                '+', '-', '-', '+', '-', '-', '+', //
                '|', ' ', ' ', '|', ' ', ' ', '|', //
                '+', '-', '-', '+', '-', '-', '+', //
                '|', ' ', ' ', '|', ' ', ' ', '|', //
                '+', '-', '-', '+', '-', '-', '+',
            ]
        );
    }
}
