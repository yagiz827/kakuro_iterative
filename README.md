# Kakuro Solver (Iterative Backtracking)

A **Kakuro puzzle solver in C++** that uses **iterative backtracking with an explicit stack** instead of recursion. Kakuro is like a crossword with numbers: every run of empty cells must add up to its clue, using the digits 1–9 without repeating a digit within a run.

## How it works

- The board is read from a `.kakuro` file. The first line holds the dimensions, then the grid follows, with clue cells and empty cells.
- The row and column clues (target sums) are extracted and matched to their runs.
- The solver fills the empty cells one by one, pushing `(cell, digit)` pairs onto a **`std::stack`**:
  - a digit is **accepted** if it keeps the row and column valid: no repeated digits, a partial sum that stays below the clue while cells are still empty, and a sum that exactly equals the clue once the run is full
  - if no digit fits, the solver **pops** the stack and tries the next digit in the previous cell (backtracking)
  - the puzzle is solved when the stack holds one entry per empty cell
- The solve time is measured with `std::chrono`, and the solution is written to `visualize.kakuro`.

**Why iterative?** A recursive solver can overflow the call stack on large boards. An explicit stack keeps memory on the heap and makes the algorithm's state easy to inspect. This version was written for a university parallel-computing course, as the sequential baseline before parallelisation.

## Built with

C++ · STL (`stack`, `vector`, `map`) · `chrono` · Visual Studio

## Running

Open `ConsoleApplication2.sln` in Visual Studio. Set the `filename` variable in `main()` to your `.kakuro` board file, then build and run. The program prints the solve time and the solved grid, and writes it to `visualize.kakuro`.
