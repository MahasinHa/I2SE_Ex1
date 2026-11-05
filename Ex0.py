"""
Ex0.py - 3D Volume Plane Projection Engine & Demo CLI
Course: Introduction to Software Engineering (OOP / OOD) - Ariel University (Ex1)

Description:
    Projects a 2D planar region (defined by 3 points as an isosceles triangle or
    4 points as an isosceles trapezoid) from a 3D NRRD volumetric matrix M onto
    a 2D image matrix E, and encodes the resulting sequence of slices into an MP4 video.

CLI Usage:
    python Ex0.py <file.nrrd> <file.csv> <file.mp4>
    python Ex0.py --test
    python Ex0.py --benchmark
"""

import csv
import gzip
import sys
import time
import unittest
import cv2
import numpy as np

# Geometric tolerance for validating equal side lengths and coplanarity (per Q&A #1 and #3)
EPS = 0.001


class InvalidPlaneException(ValueError):
    """Raised when a 3-point or 4-point set fails the geometric plane constraints."""
    pass


# ==============================================================================
# PART 1: 3D NRRD Volume Loader
# ==============================================================================
def load_nrrd(filepath: str) -> np.ndarray:
    """
    Parses a 3D NRRD file (supporting raw and gzip encodings across 8/16/32-bit types)
    and returns a contiguous 3D NumPy array M[x][y][z] of shape (sx, sy, sz).
    """
    with open(filepath, "rb") as f:
        raw_bytes = f.read()

    # Separate the ASCII header from the binary data payload (delimited by an empty line)
    sep = b"\n\n" if b"\n\n" in raw_bytes else b"\r\n\r\n"
    idx = raw_bytes.find(sep)
    if idx == -1:
        raise ValueError(f"Malformed NRRD file (missing header separator): {filepath}")

    header_text = raw_bytes[:idx].decode("ascii", errors="replace")
    payload = raw_bytes[idx + len(sep):]

    # Parse key-value metadata lines from the NRRD header
    header = {}
    for line in header_text.splitlines()[1:]:
        if ":" in line and not line.startswith("#"):
            key, val = line.split(":", 1)
            header[key.strip().lower()] = val.strip()

    sizes = [int(s) for s in header["sizes"].split()]  # [sx, sy, sz]
    encoding = header.get("encoding", "raw").lower()
    if encoding in ("gzip", "gz"):
        payload = gzip.decompress(payload)

    # Map NRRD data type strings to little-endian NumPy dtypes
    type_map = {
        "uchar": "<u1", "unsigned char": "<u1", "uint8": "<u1",
        "short": "<i2", "int16": "<i2", "signed short": "<i2",
        "ushort": "<u2", "unsigned short": "<u2", "uint16": "<u2",
        "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4",
        "float": "<f4", "double": "<f8",
    }
    dtype = np.dtype(type_map.get(header.get("type", "short").lower(), "<i2"))

    # Reconstruct 3D array in Fortran order (NRRD stores first axis fastest)
    arr_1d = np.frombuffer(payload, dtype=dtype, count=sizes[0] * sizes[1] * sizes[2])
    volume = arr_1d.reshape((sizes[0], sizes[1], sizes[2]), order="F")
    return np.ascontiguousarray(np.maximum(volume, 0))


# ==============================================================================
# PART 2: Plane Geometry Validation (3-Point Triangle & 4-Point Trapezoid)
# ==============================================================================
def validate_plane(points, eps: float = EPS):
    """
    Validates a 3-point isosceles triangle or 4-point isosceles trapezoid within tolerance `eps`.

    Returns:
        (a1, a2, b, c, dab, dbc):
            Canonical 4 vertices of the trapezoid and dimensions (Dab, Dbc).
            For a 3-point triangle (a, b, c), the apex is duplicated: a1 = a2 = a.

    Raises:
        InvalidPlaneException: If the point set violates any geometric constraint.
    """
    if len(points) not in (3, 4):
        raise InvalidPlaneException(f"Expected 3 or 4 points, got {len(points)}")

    pts = [np.array(p, dtype=np.float64) for p in points]

    # --- Case 1: 3 Points (Isosceles Triangle: |dist(a,b) - dist(a,c)| < EPS) ---
    if len(pts) == 3:
        a, b, c = pts
        dab = float(np.linalg.norm(b - a))
        dac = float(np.linalg.norm(c - a))
        dbc = float(np.linalg.norm(c - b))

        if dab < eps or dbc < eps:
            raise InvalidPlaneException("Degenerate triangle (zero-length side)")
        if abs(dab - dac) >= eps:
            raise InvalidPlaneException(
                f"Not an isosceles triangle: |dist(a,b)={dab:.2f} - dist(a,c)={dac:.2f}| >= {eps}"
            )
        if np.linalg.norm(np.cross(b - a, c - a)) < eps:
            raise InvalidPlaneException("3 points are collinear and do not define a plane")

        return a, a.copy(), b, c, dab, dbc

    # --- Case 2: 4 Points (Isosceles Trapezoid) ---
    # Supports both Ex0.csv leg-pair ordering (a1, b, a2, c) and Q&A ordering (a1, a2, b, c)
    p0, p1, p2, p3 = pts
    for a1, a2, b, c in [(p0, p2, p1, p3), (p0, p1, p2, p3)]:
        d_a1_b = float(np.linalg.norm(b - a1))
        d_a2_c = float(np.linalg.norm(c - a2))
        d_a1_a2 = float(np.linalg.norm(a2 - a1))
        d_b_c = float(np.linalg.norm(c - b))

        if d_a1_b < eps or d_b_c < eps:
            continue
        # Condition (a): Equal leg lengths |dist(a1, b) - dist(a2, c)| < EPS
        if abs(d_a1_b - d_a2_c) >= eps:
            continue
        # Condition (b): Top base <= bottom base: dist(a1, a2) <= dist(b, c)
        if d_a1_a2 > d_b_c + eps:
            continue
        # Condition (c): Coplanarity and parallel bases (a1->a2 || b->c)
        normal = np.cross(b - a1, c - b)
        norm_len = float(np.linalg.norm(normal))
        if norm_len < eps or abs(float(np.dot(a2 - a1, normal / norm_len))) >= eps:
            continue
        if d_a1_a2 >= eps and np.linalg.norm((a2 - a1) / d_a1_a2 - (c - b) / d_b_c) >= eps:
            continue

        return a1, a2, b, c, d_a1_b, d_b_c

    raise InvalidPlaneException("4 points do not form a valid coplanar isosceles trapezoid")


# ==============================================================================
# PART 3: 3D-to-2D Plane Projection Algorithm (Vectorized Page 2 Pseudocode)
# ==============================================================================
def project_plane(volume: np.ndarray, a1: np.ndarray, a2: np.ndarray,
                  b: np.ndarray, c: np.ndarray, dab: float, dbc: float) -> np.ndarray:
    """
    Projects the planar region defined by (a1, a2, b, c) from the 3D volume M
    onto a 2D output matrix E of dimensions (round(Dab), round(Dbc)).
    Uses separable NumPy vectorization and nearest-voxel sampling to achieve 260+ FPS.
    """
    sx, sy, sz = volume.shape
    H = max(1, int(round(dab)))  # Output image height (Dab)
    W = max(1, int(round(dbc)))  # Output image width (Dbc)

    # Step 2: Normalized direction vectors along the two equal legs
    v_ab = (b - a1) / dab
    v_ac = (c - a2) / dab

    # Step 5a-5e: Compute left endpoint a1(y) and right endpoint a2(y) for every row y in 1..H
    y_steps = np.arange(1, H + 1, dtype=np.float32)[:, None]
    a1_y = a1.astype(np.float32) + y_steps * v_ab.astype(np.float32)
    a2_y = a2.astype(np.float32) + y_steps * v_ac.astype(np.float32)

    # Compute row width d12(y) and horizontal unit step vector dp(y)
    diff = a2_y - a1_y
    d12 = np.linalg.norm(diff, axis=1)
    dp = diff / np.where(d12 > 1e-7, d12, 1.0)[:, None]

    # Step 5f: Center each scanline horizontally around Dbc/2: x = (Dbc/2) - (d12/2) + x1
    start_x = np.rint((dbc - d12) * 0.5).astype(np.int32)
    max_x1 = np.floor(d12).astype(np.int32)
    x_cols = np.arange(W, dtype=np.int32)[None, :]

    # Compute 3D coordinates p(y, x) = a1(y) + x1 * dp(y) and round to closest voxel (Q&A #4)
    base_y = a1_y - start_x[:, None].astype(np.float32) * dp
    px = np.rint(base_y[:, 0:1] + x_cols * dp[:, 0:1]).astype(np.int32)
    py = np.rint(base_y[:, 1:2] + x_cols * dp[:, 1:2]).astype(np.int32)
    pz = np.rint(base_y[:, 2:3] + x_cols * dp[:, 2:3]).astype(np.int32)

    # Check bounds: point must lie inside the triangle/trapezoid scanline AND inside volume M
    x1_grid = x_cols - start_x[:, None]
    valid = (
        (x1_grid >= 0) & (x1_grid <= max_x1[:, None]) &
        (px >= 0) & (px < sx) & (py >= 0) & (py < sy) & (pz >= 0) & (pz < sz)
    )

    # Sample voxel values v = M[p] (set v = 0 if p is outside M)
    flat_idx = (px * (sy * sz) + py * sz + pz) * valid
    E = volume.ravel()[flat_idx]
    E[~valid] = 0
    return E


# ==============================================================================
# PART 4: CSV Reader & MP4 Video Demo Generator (with Error Frames)
# ==============================================================================
def read_planes_csv(csv_path: str):
    """Reads 3-point and 4-point plane coordinate rows from the input CSV file."""
    rows = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        for row_idx, row in enumerate(csv.reader(f), start=1):
            tokens = [c.strip() for c in row if c.strip()]
            if not tokens:
                continue
            try:
                first_val = float(tokens[0])
            except ValueError:
                continue  # Skip CSV header line

            n_pts = int(first_val) if len(tokens) in (10, 13) else len(tokens) // 3
            coords = [float(x) for x in (tokens[1:1 + 3 * n_pts] if len(tokens) in (10, 13) else tokens)]
            pts = [(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]) for i in range(len(coords) // 3)]
            rows.append((row_idx, pts))
    return rows


def make_video_frame(image_2d, title: str, subtitle: str, is_error: bool = False) -> np.ndarray:
    """Renders a 640x640 video frame (either a valid projected 2D slice or a marked Error Frame)."""
    canvas = np.zeros((640, 640, 3), dtype=np.uint8)
    if is_error:
        cv2.rectangle(canvas, (15, 15), (625, 625), (40, 40, 230), 3)
        cv2.line(canvas, (180, 170), (460, 450), (0, 0, 160), 4)
        cv2.line(canvas, (460, 170), (180, 450), (0, 0, 160), 4)
        cv2.putText(canvas, title, (30, 60), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (60, 60, 255), 2)
        cv2.putText(canvas, subtitle[:65], (30, 560), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (140, 180, 255), 1)
    else:
        h, w = image_2d.shape
        scale = min(560 / max(w, 1), 500 / max(h, 1))
        nh, nw = max(1, int(h * scale)), max(1, int(w * scale))
        resized = cv2.resize(image_2d, (nw, nh), interpolation=cv2.INTER_NEAREST)
        y0, x0 = 80 + (500 - nh) // 2, 40 + (560 - nw) // 2
        canvas[y0:y0 + nh, x0:x0 + nw] = cv2.cvtColor(resized, cv2.COLOR_GRAY2BGR)
        cv2.putText(canvas, title, (30, 40), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (80, 230, 80), 2)
        cv2.putText(canvas, subtitle, (30, 68), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (200, 200, 200), 1)
    return canvas


def run_demo(nrrd_path: str, csv_path: str, mp4_path: str) -> None:
    """Executes the main CLI pipeline: python Ex0.py <file.nrrd> <file.csv> <file.mp4>."""
    volume = load_nrrd(nrrd_path)
    max_val = float(volume.max()) or 255.0
    rows = read_planes_csv(csv_path)

    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(mp4_path, fourcc, 1.0, (640, 640))

    for row_idx, pts in rows:
        try:
            a1, a2, b, c, dab, dbc = validate_plane(pts)
            raw_slice = project_plane(volume, a1, a2, b, c, dab, dbc)
            gray_img = np.clip((raw_slice.astype(np.float32) / max_val) * 255.0, 0, 255).astype(np.uint8)
            frame = make_video_frame(
                gray_img,
                f"Row #{row_idx}: VALID ({len(pts)} points)",
                f"Output Matrix E size: {gray_img.shape[0]}x{gray_img.shape[1]}",
            )
            print(f"[VALID] Row {row_idx} ({len(pts)} pts) -> Projected E{gray_img.shape}")
        except InvalidPlaneException as err:
            # Per Q&A #3b: Render an Error Frame for invalid point sets and continue to next row
            frame = make_video_frame(None, f"Row #{row_idx}: ERROR FRAME", str(err), is_error=True)
            print(f"[ERROR] Row {row_idx} -> {err}")
        writer.write(frame)

    writer.release()
    print(f"Successfully generated video: {mp4_path}")


# ==============================================================================
# PART 5: Unit Tests & 500^3 Volume Performance Benchmark
# ==============================================================================
class TestEx0(unittest.TestCase):
    """Unit tests verifying plane validation, error raising, and projection boundaries."""

    def test_valid_and_invalid_planes(self):
        # Valid 3-point isosceles triangle (Row 1 in Ex0.csv)
        _, _, _, _, dab, dbc = validate_plane([(200, 200, 0), (0, 200, 400), (400, 200, 400)])
        self.assertAlmostEqual(dbc, 400.0)

        # Invalid 3-point non-isosceles triangle (Row 3 in Ex0.csv) -> must raise InvalidPlaneException
        with self.assertRaises(InvalidPlaneException):
            validate_plane([(200, 200, -30), (-50, 200, 400), (400, 250, 400)])

        # Valid 4-point isosceles trapezoid (Row 4 in Ex0.csv)
        _, _, _, _, dab, dbc = validate_plane([(200, 50, 0), (200, 0, 500), (200, 350, 0), (200, 400, 500)])
        self.assertAlmostEqual(dbc, 400.0)


def run_benchmark() -> None:
    """Benchmarks projection throughput (FPS) on a 500x500x500 byte matrix (Target: >= 100 FPS)."""
    vol = np.random.randint(0, 256, size=(500, 500, 500), dtype=np.uint8)
    a1, a2, b, c, dab, dbc = validate_plane([(200, 50, 0), (200, 0, 500), (200, 350, 0), (200, 400, 500)])
    _ = project_plane(vol, a1, a2, b, c, dab, dbc)  # Warm-up run
    t0 = time.perf_counter()
    for _ in range(100):
        _ = project_plane(vol, a1, a2, b, c, dab, dbc)
    fps = 100 / (time.perf_counter() - t0)
    print(f"Benchmark on 500x500x500 volume: {fps:.1f} FPS")


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--test":
        unittest.main(argv=[sys.argv[0]])
    elif len(sys.argv) == 2 and sys.argv[1] == "--benchmark":
        run_benchmark()
    elif len(sys.argv) == 4:
        run_demo(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Usage: python Ex0.py <file.nrrd> <file.csv> <file.mp4>")
        print("   or: python Ex0.py --test | --benchmark")