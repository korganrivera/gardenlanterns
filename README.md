# Garden Lanterns — `lantern_ga.c`

Allele-based genetic algorithm (GA) for optimizing lantern placement inside a garden polygon.  
The program searches for lantern positions that maximize *minimum illuminance (lux)* across the garden while penalizing excessive lantern count.

This is a simulation-driven design tool: you define the garden shape in coordinates, tune physical parameters (lumens, height, grid resolution), and let evolution grind toward a balanced lighting layout.

---

## Overview

This program models lantern placement as a diploid genetic system:

Each lantern “slot” has two alleles, each carrying:
- a coordinate (x, y),
- a dominance marker (`D` or `R`),
- or an OFF state.

Expression rules determine whether each slot produces an active lantern.
Fitness is computed as:

    fitness = minimum_lux_across_grid - 0.02 × number_of_lanterns

The GA evolves populations of candidate layouts toward maximal light uniformity using tournament selection, crossover, and mutation.

---

## Build

Requires GCC or Clang with C99 support.

```bash
gcc -O3 -std=c99 lantern_ga.c -o lantern_ga -lm
```

Optional debug build:

```bash
gcc -std=c99 -O0 -g -Wall -Wextra lantern_ga.c -o lantern_dbg -lm
```

---

## Run

The program prompts for how many lantern slots to optimize (default is 5).

```bash
./lantern_ga
```

Example output:

```
Grid built: 347 interior sample points (NX=100 NY=100)
Running GA for N_slots=5 ... (POP=48, GENS=200)
gen 0: best_score -0.11203 best_minlux 0.0000
gen 50: best_score 0.61241 best_minlux 0.5853
gen 100: best_score 0.68433 best_minlux 0.6231

Best found (min lux): 0.6231 lux | coverage >= 1.00 lux: 72.48% | lanterns used: 4
Lantern positions (feet):
   1: x =  98.12 ft, y =  54.33 ft
   2: x = 120.45 ft, y = 100.34 ft
   3: x =  65.22 ft, y = 140.01 ft
   4: x =  32.00 ft, y =  45.10 ft
Wrote heatmap CSV to heatmap.csv (347 points)
```

---

## Output

### Console
- Minimum lux found
- Percentage of grid points ≥ target lux
- Final lantern coordinates (in feet)

### Heatmap CSV (optional)
If `SAVE_CSV_HEATMAP` is enabled:

`heatmap.csv` columns:

    x_ft,y_ft,lux

Each row is a sampled interior grid point and the lux value after evolution.

This can be visualized in Python, spreadsheets, or GIS tools as a heatmap.

---

## Physical Model

- Illuminance uses inverse-square law:

      lux = lumens / (4πr²)

- Distance `r` is the full 3D distance (horizontal + mounting height)
- Lanterns emit uniformly in all directions
- No occlusion, reflection, or terrain modeling is included

---

## Key Parameters (Top of File)

Edit `#define` values and recompile:

| Macro | Description |
|--------|-------------|
| `POWER_LM` | Lumens per lantern |
| `LANTERN_HEIGHT_FT` | Mount height in feet |
| `TARGET_LUX` | Coverage threshold |
| `GRID_NX / GRID_NY` | Evaluation resolution |
| `POP_SIZE` | GA population |
| `GENERATIONS` | GA iterations |
| `MAX_SLOTS` | Max possible lantern positions |
| `TOURNAMENT_K` | Tournament size |
| `CROSSOVER_RATE` | Genetic crossover rate |
| `MUTATION_RATE_DOM` | Dominant mutation rate |
| `MUTATION_RATE_REC` | Recessive mutation rate |
| `ELITISM` | Survivors per generation |
| `RANDOM_SEED` | RNG seed |
| `SAVE_CSV_HEATMAP` | Write CSV output |

---

## Garden Shape

The garden boundary is currently a hard-coded quadrilateral:

```c
double verts_ft[4][2]
```

Coordinates are provided in feet and converted internally to meters.  
Edit this array to use your own plot layout.

---

## Algorithm Notes

- Grid sampling approximates area lighting coverage
- Dominant alleles override recessive ones
- Alleles may mutate by:
  - moving position,
  - toggling on/off,
  - changing dominance type
- Elitism preserves high performers
- Selection pressure increases spatial uniformity
- The model rewards fewer lamps when coverage is similar

---

## Performance Tips

- Use coarse grids for testing (`50×50`)
- Increase grid size for final runs
- Increase population size for difficult layouts
- Run multiple seeds for stochastic confirmation
- Profile with `-pg` or sanitizers for optimization

---

## Future Enhancements

Suggested expansions:
- CLI flags
- JSON input for polygons
- Directional light patterns
- Obstacle handling
- Multi-threading
- Visualization output formats (PNG, GeoJSON, KML)
- Real fixture photometrics

---

## License

Choose and include a license file (MIT, Apache-2.0, GPL, etc.)

---

## Author

korgan  
GitHub: https://github.com/korganrivera

---