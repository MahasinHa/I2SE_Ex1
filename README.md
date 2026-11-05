# Ex1: 3D Volume Plane Projection Engine (`Ex0.py`)

**Course:** Introduction to Software Engineering (OOP / OOD) — Ariel University

**Assignment:** Ex1 — Algorithm: Implementing a Plane Projection of a 3D Volume

## 1. Problem Definition & Project Scope

Given a 3D volumetric raster data matrix $M$ of dimensions $s_x \times s_y \times s_z$ (where $1 \le s_x, s_y, s_z < 2^{10}$ and each voxel $M[x][y][z]$ is an 8-bit, 16-bit, or 32-bit unsigned integer read from a [3D Slicer](https://www.slicer.org/) NRRD file) and a 2D planar region $p$ embedded in 3D space, the objective is to project the voxels of $M$ intersecting $p$ onto a 2D output image matrix $E$.

Because a mathematical plane is infinite while a 2D digital image is finite and rectangular, this project formalizes $p$ as a bounded planar geometric primitive defined by either:

* **3 Points** $(a, b, c)$: An **isosceles triangle** where $\vert{}dist(a,b) - dist(a,c)\vert{} < \epsilon$ ($\epsilon = 0.001$).

* **4 Points** $(a_1, a_2, b, c)$: An **isosceles trapezoid** on the same plane where $\vert{}dist(a_1,b) - dist(a_2,c)\vert{} < \epsilon$, $dist(a_1,a_2) \le dist(b,c)$, and the legs $\vec{a_1 b}$ and $\vec{a_2 c}$ are coplanar.

The project provides:

1. A zero-dependency NRRD v4 parser supporting both `raw` and `gzip` encodings across 8-bit, 16-bit, and 32-bit voxel types.

2. Strict geometric validation (`PlaneValidator`) that raises an `InvalidPlaneException` whenever a 3-point or 4-point specification violates the $\epsilon = 0.001$ constraints.

3. A vectorized NumPy projection engine (`PlaneProjector`) achieving **345+ FPS** on a $500^3$ volume (exceeding both the 30 FPS minimum requirement and the 100 FPS expected target).

4. A CLI video generator (`python Ex0.py <file.nrrd> <file.csv> <file.mp4>`) that renders valid plane projections alongside diagnostic **Error Frames** for invalid CSV rows.

## 2. Software Requirements Specification (SRS)

Following the [Software Requirements Specification](https://en.wikipedia.org/wiki/Software_requirements_specification) standard, the system requirements are organized into functional and non-functional specifications.

### 2.1 Functional Requirements

| ID | Requirement Name | Specification | 
| ----- | ----- | ----- | 
| **FR-1** | **NRRD Volume Loading** | Parse 3D NRRD volumetric files (`NRRD0004`, `raw` or `gzip` encoding, little/big endian, 8/16/32-bit integers) into a contiguous 3D matrix $M[x][y][z]$ with $0 \le x < s_x, 0 \le y < s_y, 0 \le z < s_z$. | 
| **FR-2** | **CSV Plane Parsing** | Read a CSV list of plane coordinate rows, where each row defines either 3 points $(a, b, c)$ or 4 points $(a_1, b, a_2, c)$ / $(a_1, a_2, b, c)$. | 
| **FR-3** | **3-Point Validation** | Verify that 3-point definitions satisfy $\vert{}dist(a,b) - dist(a,c)\vert{} < 0.001$ and are non-collinear. Raise `InvalidPlaneException` on violation. | 
| **FR-4** | **4-Point Validation** | Verify that 4-point definitions satisfy $\vert{}dist(a_1,b) - dist(a_2,c)\vert{} < 0.001$, $dist(a_1,a_2) \le dist(b,c) + 0.001$, coplanarity within $\epsilon = 0.001$, and parallel bases $\vec{a_1 a_2} \parallel \vec{b c}$. Raise `InvalidPlaneException` on violation. | 
| **FR-5** | **Plane Projection & Bounds Handling** | Compute the 2D matrix $E$ of dimensions $D_{ab} \times D_{bc}$ by sampling the nearest voxel (`round(p)`) in $M$, assigning $v = 0$ whenever a 3D point $p$ falls outside $[0, s_x) \times [0, s_y) \times [0, s_z)$. | 
| **FR-6** | **MP4 Demo & Error Frames** | Output an `.mp4` video containing one frame per CSV row in list order. Catch `InvalidPlaneException` during batch demo execution and render a marked black Error Frame so remaining valid planes still render properly. | 

### 2.2 Non-Functional Requirements

| ID | Category | Target Metric | Achieved Result | 
| ----- | ----- | ----- | ----- | 
| **NFR-1** | **Throughput (FPS)** | $\ge 30$ FPS mandatory, $\ge 100$ FPS expected on a $500^3$ byte matrix | **345.3 FPS** ($2.90\text{ ms/frame}$) | 
| **NFR-2** | **Numerical Precision** | Geometric tolerance $\epsilon = 0.001$; nearest-voxel sampling equivalence with reference pseudocode | 100% pixel-identical match (`np.testing.assert_array_equal`) | 
| **NFR-3** | **Modularity & OOD** | Object-oriented separation of concerns (`Volume3D`, `PlaneGeometry`, `PlaneValidator`, `PlaneProjector`, `CSVPlaneReader`, `VideoRenderer`) | Single-responsibility classes with unit tests | 

## 3. Detailed Algorithm Design

### 3.1 Unifying Triangles and Trapezoids

Notice that an isosceles triangle $(a, b, c)$ with apex $a$ is simply a degenerate isosceles trapezoid $(a_1, a_2, b, c)$ where the top base has zero length ($a_1 = a_2 = a$ and $dist(a_1, a_2) = 0 \le dist(b, c)$). By normalizing every valid 3-point triangle into a `PlaneGeometry` object with $a_1 = a$ and $a_2 = a$, the exact same projection kernel handles both 3-point and 4-point planes without branching.

### 3.2 Mathematical Derivation of Separable Vectorization

In the assignment's baseline pseudocode:

1. $D_{ab} = dist(a_1, b) = dist(a_2, c)$ and $D_{bc} = dist(b, c)$.

2. Unit leg vectors: $\vec{v}_{ab} = \frac{b - a_1}{D_{ab}}$ and $\vec{v}_{ac} = \frac{c - a_2}{D_{ab}}$.

3. For each vertical step $y \in \{0, 1, \dots, \text{round}(D_{ab}) - 1\}$ (representing step count $k = y + 1$):
   

   $$
   \vec{a}_1(y) = a_1 + (y + 1)\vec{v}_{ab}, \quad \vec{a}_2(y) = a_2 + (y + 1)\vec{v}_{ac}
   $$

4. Across row $y$, the horizontal span is $d_{12}(y) = \Vert{}\vec{a}_2(y) - \vec{a}_1(y)\Vert{}_2$, and the unit scan direction is:
   

   $$
   \vec{dp}(y) = \frac{\vec{a}_2(y) - \vec{a}_1(y)}{d_{12}(y)}
   $$

5. For each horizontal offset $x_1 \in \{0, 1, \dots, \lfloor d_{12}(y) \rfloor\}$, the output column index $x$ is centered around $D_{bc}/2$:
   

   $$
   x = \text{round}\left(\frac{D_{bc} - d_{12}(y)}{2}\right) + x_1 = x_{\text{start}}(y) + x_1
   $$

   
   Substituting $x_1 = x - x_{\text{start}}(y)$ into the 3D sample position $\vec{p}(y, x_1) = \vec{a}_1(y) + x_1 \vec{dp}(y)$ yields a **separable linear equation** in $y$ and $x$:
   

   $$
   \vec{p}(y, x) = \underbrace{\Big(\vec{a}_1(y) - x_{\text{start}}(y)\vec{dp}(y)\Big)}_{\vec{\text{base}}(y) \text{ — shape } (H, 3)} + x \cdot \vec{dp}(y)
   $$

### 3.3 Why This Achieves 345+ FPS in Python

A naive nested Python `for` loop executes $\sim 200,000$ iterations per frame ($\sim 0.3\text{ s/frame}$, or $3\text{ FPS}$). By expressing $\vec{p}(y, x)$ in separable form, we broadcast a 1D column vector $\vec{\text{base}}(y)$ of length $H$ against a 1D row vector $x \in [0, W)$ for each axis $(p_x, p_y, p_z)$, round to the nearest integer voxel index (`np.rint`), and perform a single 1D C-contiguous stride lookup:

$$
\text{flat\_idx} = p_x \cdot (s_y s_z) + p_y \cdot s_z + p_z
$$


Out-of-bounds coordinates and pixels outside the triangle/trapezoid mask ($x_1 < 0$ or $x_1 > \lfloor d_{12}(y) \rfloor$) are zeroed in one vectorized boolean mask operation.

## 4. Dataset of NRRD Files & Sample CSV Planes

1. **Primary Volumetric Dataset (`MRHead.nrrd` / [MRHead_1.nrrd](https://drive.google.com/file/d/1cqqOgSf8k4e14K5CBs-V_EvMmATX1K-J/view?usp=sharing))**:

   * Source: [3D Slicer](https://www.slicer.org/) Sample Data (`MRHead`).

   * Dimensions: $256 \times 256 \times 130$ voxels (`type: short`, `encoding: gzip`, `endian: little`), voxel intensity range $[0, 279]$.

2. **High-Resolution Benchmark Dataset (Synthetic** $500^3$ **Volume)**:

   * Generated programmatically via `python Ex0.py --benchmark` ($500 \times 500 \times 500 = 125,000,000$ bytes, `uint8`) to benchmark worst-case throughput required by Page 1.

3. **Sample Plane Coordinates (`Ex0.csv` / [Official Ex0.csv Spreadsheet](https://docs.google.com/spreadsheets/d/1Tfy1yMYhAJjRBA8ZpRp_07pZXqGtxWhWoP02yZW1FQE/edit?usp=sharing))**:

   * **Data Row 1 (CSV Line 2)**: `3, 200,200,0, 0,200,400, 400,200,400` — **Valid** 3-point isosceles triangle ($D_{ab} = D_{ac} = 447.2136$, $D_{bc} = 400.0$).

   * **Data Row 2 (CSV Line 3)**: `4, 200,0,0, 200,0,500, 200,400,0, 200,400,500` — **Valid** 4-point rectangle/trapezoid ($D_{ab} = 500.0$, $D_{bc} = 400.0$).

   * **Data Row 3 (CSV Line 4)**: `3, 200,200,-30, -50,200,400, 400,250,400` — **Invalid** 3-point triangle ($dist(a,b) = 497.3932 \ne dist(a,c) = 476.8648$, difference $20.5284 \ge 0.001$). Triggers `InvalidPlaneException` and produces an **Error Frame** in the output video.

   * **Data Row 4 (CSV Line 5)**: `4, 200,50,0, 200,0,500, 200,350,0, 200,400,500` — **Valid** 4-point isosceles trapezoid ($D_{ab} = 502.4938$, top base $= 300.0 \le$ bottom base $= 400.0$).

## 5. Testing Strategy

Our testing strategy (`TestPlaneProjection` inside `Ex0.py`, executed via `python Ex0.py --test`) verifies four layers:

1. **Geometric Validation Unit Tests**:

   * Valid 3-point isosceles triangles and 4-point rectangles/trapezoids pass within $\epsilon = 0.001$.

   * Asymmetric triangles ($\vert{}dist(a,b) - dist(a,c)\vert{} \ge \epsilon$), non-coplanar 4-point sets, collinear points, and inverted trapezoids ($dist(a_1, a_2) > dist(b, c)$) raise `InvalidPlaneException`.

2. **Algorithmic Equivalence Tests**:

   * Compares `PlaneProjector.project_vectorized` against `PlaneProjector.project_reference` (which executes the literal nested loops from Page 2) on both 3-point triangles and 4-point trapezoids and asserts 100% pixel-by-pixel equality.

3. **Boundary & Out-of-Bounds Tests**:

   * Verifies that planes partially or completely outside $[0, s_x) \times [0, s_y) \times [0, s_z)$ safely assign $v = 0$ without index errors.

4. **End-to-End Demo & Fault-Isolation Test**:

   * Runs `python Ex0.py MRHead.nrrd Ex0.csv Ex0_demo.mp4` and confirms that Row 3's invalid triangle is caught and rendered as an Error Frame while Rows 1, 2, and 4 render clean anatomical MRI cross-sections.

## 6. How to Run & Experimental Results

### 6.1 CLI Commands

```
# 1. Run the demo pipeline to produce the MP4 video
python3 Ex0.py MRHead.nrrd Ex0.csv Ex0_demo.mp4

# 2. Run the unit test suite
python3 Ex0.py --test

# 3. Run the 500x500x500 volume FPS benchmark
python3 Ex0.py --benchmark

```

### 6.2 Experimental Benchmark & Demo Results

| Experiment | Input Volume | Plane Geometry | Output Matrix $E$ | Runtime / Throughput | Outcome | 
| ----- | ----- | ----- | ----- | ----- | ----- | 
| **Synthetic** $500^3$ **Benchmark** | $500 \times 500 \times 500$ (`uint8`, 125 MB) | 4-pt Trapezoid ($D_{ab}=502.49, D_{bc}=400$) | $502 \times 400$ px | **2.90 ms/frame (345.3 FPS)** | Exceeds 100 FPS target by 3.45x | 
| **Demo Frame 1 (CSV Line 2)** | `MRHead.nrrd` ($256 \times 256 \times 130$) | 3-pt Isosceles Triangle ($x=0..400, y=200$) | $447 \times 400$ px | **3.58 ms** | Valid triangular MRI projection | 
| **Demo Frame 2 (CSV Line 3)** | `MRHead.nrrd` ($256 \times 256 \times 130$) | 4-pt Rectangle ($x=200, y=0..400, z=0..500$) | $500 \times 400$ px | **2.96 ms** | Valid coronal MRI cross-section | 
| **Demo Frame 3 (CSV Line 4)** | `MRHead.nrrd` ($256 \times 256 \times 130$) | 3-pt Invalid Triangle ($\Delta d = 20.53 \ge 0.001$) | $640 \times 640$ px | **< 0.1 ms** | Caught exception; rendered Error Frame | 
| **Demo Frame 4 (CSV Line 5)** | `MRHead.nrrd` ($256 \times 256 \times 130$) | 4-pt Isosceles Trapezoid ($x=200$, top 300, base 400) | $502 \times 400$ px | **3.49 ms** | Valid trapezoidal MRI cross-section | 
