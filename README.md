# Ex1 - 3D Volume Plane Projection

Introduction to Software Engineering - Ariel University

## Problem Description

In this assignment we get a 3D matrix `M` of size `sx * sy * sz` (from an NRRD file like in [3D Slicer](https://www.slicer.org/)) and a plane `p`, and we need to project the plane onto a 2D image `E`.

Because a plane is infinite, we define the area we want to project using 3 or 4 3D points (with `EPS = 0.001`):
- **3 points (a, b, c)**: An isosceles triangle where `|dist(a,b) - dist(a,c)| < 0.001`.
- **4 points (a1, a2, b, c)**: An isosceles trapezoid on the same plane where `|dist(a1,b) - dist(a2,c)| < 0.001` and `dist(a1,a2) <= dist(b,c)`.

## Requirements (SRS)

Based on the [SRS](https://en.wikipedia.org/wiki/Software_requirements_specification) guidelines:
- **Input**:
  - A 3D `.nrrd` file (supports 8, 16, and 32-bit values, raw or gzip).
  - A `.csv` file where each row has 3 or 4 3D points.
- **Validation & Error Handling**:
  - Check if the points in each row are valid (`EPS = 0.001`).
  - If a set of points is not valid, raise an exception (`InvalidPlaneException`).
  - In the video demo, if a row raises an exception, draw an "Error Frame" (a black frame with red text) and continue to the rest of the CSV rows.
- **Output**:
  - An `.mp4` video file containing one frame for each plane in the CSV file in order.
- **Performance**:
  - Needs to run at 30+ FPS (and ideally 100+ FPS) on a `500x500x500` matrix.

## Algorithm Explanation

An isosceles triangle `(a, b, c)` is basically the same as a trapezoid `(a1, a2, b, c)` where `a1 = a` and `a2 = a` (top width is 0). So in the code I convert 3-point inputs into 4 points with `a1 = a2 = a` and run the same algorithm for both:

1. Compute `Dab = dist(a1, b) = dist(a2, c)` and `Dbc = dist(b, c)`.
2. Compute the normalized step vectors `v_ab = (b - a1) / Dab` and `v_ac = (c - a2) / Dab`.
3. Create the 2D output matrix `E` of size `Dab * Dbc`.
4. Loop `y` from `0` to `Dab`:
   - Move `a1` by `v_ab` and `a2` by `v_ac`.
   - Find the distance `d12 = dist(a1, a2)` and normalized vector `dp = (a2 - a1) / d12`.
   - For `x1` from `0` to `d12`, start at `p = a1` and step by `dp`:
     - Round `p` to the closest voxel in `M`. If `p` is outside the bounds of `M`, set `v = 0`.
     - Center the row horizontally using `x = (Dbc / 2) - (d12 / 2) + x1` and set `E[y][x] = v`.

To reach 100+ FPS in Python without slow `for` loops, I used NumPy arrays to calculate all rows and columns at once.

## Data Files (Dataset)

1. **Brain MRI (`MRHead`)**:
   - NRRD File: [MRHead_1.nrrd](https://drive.google.com/file/d/1cqqOgSf8k4e14K5CBs-V_EvMmATX1K-J/view?usp=sharing) (`256 x 256 x 130`, 16-bit).
   - CSV File: [Ex0.csv](https://docs.google.com/spreadsheets/d/1Tfy1yMYhAJjRBA8ZpRp_07pZXqGtxWhWoP02yZW1FQE/edit?usp=sharing) (4 rows: valid triangle, valid rectangle, invalid triangle, and valid trapezoid).
2. **Onion Scan (`onion`)**:
   - NRRD File: `onion.nrrd` (`256 x 256 x 130`, 8-bit).
   - CSV File: `onion.csv` (120 rows of 4-point planes rotating around the onion).
3. **Croissant Scan (`croissant`)**:
   - NRRD File: `croissant.nrrd` (`500 x 500 x 500`, 8-bit) [criossant.nrrd](https://drive.google.com/file/d/1t2SsEvFX8V23J9SjAQkNrU_-OUCBysKW/view?usp=sharing) I uploaded it here not in github because it exceeds the size limit.
   - CSV File: `croissant.csv` (100 rows of 3-point isosceles triangles).

## Testing Strategy

I used `unittest` (`TestPlaneProjection` in `Ex0.py`) to check:
1. Valid 3-point triangles and 4-point trapezoids pass validation.
2. Invalid points (like Row 3 in `Ex0.csv` where `dist(a,b) = 497.39` and `dist(a,c) = 476.86`) raise `InvalidPlaneException`.
3. Points outside the 3D matrix bounds get `v = 0` without crashing.
4. The vectorized projection gives the exact same pixel output as the regular loop algorithm from Page 2.
5. The demo script creates the `.mp4` videos properly and shows an Error Frame for invalid rows.

## How to Run & Results

Install dependencies:
```bash
pip install numpy opencv-python
```

Run the demo videos, unit tests, or FPS benchmark:
```bash
python Ex0.py MRHead_1.nrrd Ex0.csv Ex0_demo.mp4
python Ex0.py onion.nrrd onion.csv onion.mp4
python Ex0.py croissant.nrrd croissant.csv croissant.mp4
python Ex0.py --test
python Ex0.py --benchmark
```

### Results Summary
- **500x500x500 Benchmark**: Runs at **155.7 FPS** (~6.4 ms per frame), which passes the 30–100+ FPS requirement.
- **Brain MRI Video (`Ex0_demo.mp4`)**:
  - Frame 1 (Row 1 - 3 points): Valid triangle projection (`447 x 400`).
  - Frame 2 (Row 2 - 4 points): Valid rectangle slice (`500 x 400`).
  - Frame 3 (Row 3 - 3 points): Invalid triangle -> shows the Error Frame.
  - Frame 4 (Row 4 - 4 points): Valid trapezoid slice (`502 x 400`).
- **Onion Video (`onion.mp4`)**:
  - Projected all 120 4-point planes (`320 x 255` px) at **~1268 FPS**.
- **Croissant Video (`croissant.mp4`)**:
  - Projected all 100 3-point triangle planes (`560 x 512` px) on the `500x500x500` volume at **~180 FPS**.
