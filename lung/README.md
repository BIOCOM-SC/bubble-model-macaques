# lung/ - the virtual macaque lung

The model's lung: the surfaces of the 7 lobes and 20 segments of a rhesus macaque lung, the 8233
terminal bronchioles of a bronchial tree grown inside it, and the distances between terminals
through the tree. All lengths in mm.

Its construction is described in the Methods of the paper.

Before running the model, unpack the distance file once (it is compressed to fit in the repository):

```sh
xz -dk lung/macaque/D.bin.xz
```

## Files (`macaque/`)

| file | content |
|---|---|
| `lobes/L1.txt` ... `L7.txt` | lobe surfaces (RUL, RML, RLL, LIC, LUL, LML, LLL), triangle meshes |
| `segments/S1.txt` ... `S20.txt` | segment surfaces, same format |
| `terminals.txt` | the terminals of the bronchial tree |
| `D.bin.xz` | distances between terminals through the bronchial tree (`D.bin` once unpacked) |

**Surfaces.** First line: number of triangles; then, for each triangle, its 3 vertices (x y z),
the first vertex repeated (to draw closed triangles), and blank lines.

**`terminals.txt`.** First line: number of terminals (8233); then one line per terminal:

```
x, y, z, segment, d_segment, lobe, d_lobe, diameter, p_initial_infection
```

- `segment` (1-20), `lobe` (1-7): segment and lobe of the terminal.
- `d_segment`, `d_lobe`: distance from the terminal to the segment / lobe wall (nearest vertex of
  the surface). The maximum radius of a lesion growing at the terminal is capped by `d_lobe`
  (model option `-dlobe`, used in the paper) or `d_segment`.
- `diameter`: diameter of the terminal bronchiole.
- `p_initial_infection`: probability that an initial lesion starts at this terminal: the same for
  all terminals of a lobe, summing per lobe to the lobe frequency of the week-3 consolidations
  (`data/processed/lobe_frequency_week3.csv`).

**`D.bin`.** Binary, little-endian doubles: the upper triangle (diagonal included) of the
symmetric 8233 x 8233 matrix of distances through the tree, row by row
(d(0,0), d(0,1), ..., d(0,N-1), d(1,1), ...), N(N+1)/2 values.

`analysis/lung_qc.py` checks this geometry (terminals outside their own lobe mesh, distances to
the walls, terminal diameters; Suppl. Fig. 2).
