# SuperCell Builder — Arbitrary Zone-Axis Supercells for STEM Simulations

Generate crystallographically correct, orthogonal supercells in **any zone axis** from a standard CIF file — ready for STEM image simulations in [abTEM](https://abtem.github.io/) and compatible tools.

---

## Why this project exists

STEM (Scanning Transmission Electron Microscopy) image simulations require an atomic structure — a "supercell" — whose **z-axis points along the electron beam direction** (the zone axis). Creating such a supercell is deceptively tricky:

- The crystal must be **rotated** into the desired orientation.
- The rotated cell must be **large enough** to fill the simulation box without gaps.
- After rotation, the cell is generally **not orthogonal**, so it must be carved into a clean cuboid.

Existing tools (abTEM, ASE's built-in utilities) handle some of these steps, but they often require deep knowledge of the material's symmetry, and edge cases — especially for low-symmetry or layered structures — frequently produce incorrect or oversized cells. This notebook solves the problem **generically**, with no material-specific knowledge required.

---

## What it does

Given a CIF file and a 3 × 3 transformation matrix describing the desired orientation, the notebook:

1. **Estimates the minimum supercell repetitions** (along the crystal's *a*, *b*, *c* axes) needed so that, after rotation and carving, the final cuboid meets your target dimensions.
2. **Builds the supercell** by tiling the primitive/conventional cell those repetitions.
3. **Applies the transformation matrix** to rotate both the cell vectors and all atomic positions into the desired zone-axis frame.
4. **Carves an orthogonal cuboid** by tiling the rotated cell into a dense atom cloud and cutting out an axis-aligned box of exactly the target size (Lx × Ly × Lz Å).
5. **Centre-crops** to the exact target dimensions, avoiding any edge bias.
6. **Saves the result** as an extended XYZ file, compatible with abTEM and other simulation codes.

---

## How it works — step by step

### Step 1 — Repetition estimation

The notebook reads the unit cell vectors from the CIF and simulates the bounding box that results from tiling by *(na, nb, nc)* repetitions and then applying the transformation. It computes a **sensitivity matrix** — how much each additional repeat along *a*, *b*, or *c* grows each Cartesian dimension of the bounding box — and uses this to find the smallest *(na, nb, nc)* that guarantees the carved cuboid will be at least as large as the target. A greedy refinement loop then tightens the estimate axis by axis.

### Step 2 — Supercell construction

ASE's `*` operator tiles the primitive cell by the computed repetitions along the crystal's own lattice vectors (not Cartesian axes), exactly matching what dedicated tools such as Atomsk would do.

### Step 3 — Orientation transformation

The transformation matrix **A** (a 3 × 3 rotation/permutation matrix obtained from VESTA) is applied to both the cell vectors and the atomic positions via a simple matrix multiplication:

```
new_cell = old_cell @ A.T
new_positions = old_positions @ A.T
```

This places the desired zone axis along z, ready for the electron beam.

### Step 4 — Cuboid carving

The rotated cell is generally a parallelepiped, not a box. To carve a clean orthorhombic cell the notebook:

- Tiles the rotated supercell in **both positive and negative** directions along every lattice vector, creating a dense atom cloud that completely fills the target bounding box on all faces.
- Shifts the cloud so its minimum corner sits at the origin.
- Keeps only atoms inside the axis-aligned box `[0, Lx) × [0, Ly) × [0, Lz)`.

Tiling in both directions is the key fix that prevents the ragged face that appears in naive implementations.

### Step 5 — Centre crop to exact target size

Because the carved box is slightly oversized (a consequence of discrete lattice repetitions), the notebook centre-crops it to the exact requested dimensions, preserving structural symmetry at both faces.

---

## Usage

### Requirements

```
Python ≥ 3.9
ase
abtem
numpy
matplotlib
```

Install with:

```bash
pip install ase abtem numpy matplotlib
```

### Getting the transformation matrix from VESTA

1. Open your CIF file in [VESTA](https://jp-minerals.org/vesta/en/).
2. Go to **Style → Orientation**.
3. Enter your desired **zone axis** (e.g. `[1 1 0]`) and an **upward vector** for the in-plane rotation.
4. Click **Apply** — the 3 × 3 transformation matrix appears above the input fields.
5. Copy the matrix and paste it into the notebook's `text` variable:

```python
text = """
 +0.000000  +1.000000  +0.000000
 +0.000000  +0.000000  +1.000000
 +1.000000  +0.000000  +0.000000
"""
```

### Running the notebook

1. Set `cif_file_path` to point to your CIF file.
2. Paste your transformation matrix from VESTA into the `text` block.
3. Set your target supercell dimensions:

```python
target_Lx = 50.0   # Å
target_Ly = 50.0   # Å
target_Lz = 100.0  # Å  (along the beam / zone axis)
```

4. Run all cells in order. The notebook prints the required repetitions, intermediate diagnostics, and writes the final structure to an extended XYZ file.

---

## Output

| File | Contents |
|---|---|
| `*_{na}x{nb}x{nc}.xyz` | Tiled supercell before rotation |
| `*_transformed.xyz` | Supercell after applying the transformation matrix |
| `*_cuboid.xyz` | Carved orthogonal cell (slightly oversized) |
| `*_cropped_{Lx}x{Ly}x{Lz}.xyz` | Final cell at exact target dimensions |

All files are in extended XYZ format and can be read directly by abTEM:

```python
import abtem, ase.io
atoms = ase.io.read("structure_cropped_50.0x50.0x100.0.xyz")
potential = abtem.Potential(atoms, sampling=0.05)
```

---

## Known limitations

- The transformation matrix must be supplied manually from VESTA.
- Very low-symmetry structures with highly oblique unit cells may require large supercells and significant memory.
- Periodic boundary conditions are set but not enforced across the carved faces; atoms sitting exactly on a boundary may or may not be included depending on floating-point rounding.

---

## Dependencies and related projects

| Library | Role |
|---|---|
| [ASE](https://wiki.fysik.dtu.dk/ase/) | CIF reading, atomic structure manipulation, XYZ I/O |
| [abTEM](https://abtem.github.io/) | STEM image simulation (downstream consumer of this notebook's output) |
| [VESTA](https://jp-minerals.org/vesta/en/) | Visualisation and transformation matrix generation |
| NumPy | All matrix operations |

---

## Contributing

Bug reports and pull requests are welcome. If you find a material for which the repetition estimator produces an incorrect or unnecessarily large supercell, please open an issue with the CIF file and transformation matrix.

---

## License

MIT
