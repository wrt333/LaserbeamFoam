# Setting Up a Custom LPBF Simulation with `laserbeamFoam`

> **中文摘要 / Chinese Summary**
>
> 本指南说明如何利用 `laserbeamFoam` 建立您自己的激光粉末床融合（LPBF）仿真案例。
> 从复制教程案例开始，按照以下各节逐步修改材料属性、激光参数、计算域网格以及
> 初始条件，最终运行求解器并在 ParaView 中后处理结果。

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Recommended Starting Point](#2-recommended-starting-point)
3. [Case Directory Structure](#3-case-directory-structure)
4. [Step 1 – Define the Computational Domain and Mesh](#4-step-1--define-the-computational-domain-and-mesh)
5. [Step 2 – Configure Material Properties](#5-step-2--configure-material-properties)
6. [Step 3 – Configure Laser Parameters](#6-step-3--configure-laser-parameters)
7. [Step 4 – Define the Laser Scan Path and Power Schedule](#7-step-4--define-the-laser-scan-path-and-power-schedule)
8. [Step 5 – Set Up Initial Conditions](#8-step-5--set-up-initial-conditions)
9. [Step 6 – Tune Solver Numerics](#9-step-6--tune-solver-numerics)
10. [Step 7 – Run the Simulation](#10-step-7--run-the-simulation)
11. [Step 8 – Post-Process and Visualise Results](#11-step-8--post-process-and-visualise-results)
12. [Optional: Powder Bed Setup with DEM](#12-optional-powder-bed-setup-with-dem)
13. [Common Issues and Tips](#13-common-issues-and-tips)

---

## 1. Prerequisites

Before setting up your own case, make sure you have completed the following:

1. **Install OpenFOAM** (v2506, v2412, or v2512 from [openfoam.com](https://www.openfoam.com)).
2. **Build `laserbeamFoam`** following the instructions in the main
   [`README.md`](../README.md).
3. **Run an existing tutorial** (e.g., `tutorials/laserbeamFoam/LPBF_small`) to
   confirm that the solver compiles and runs correctly.

---

## 2. Recommended Starting Point

The quickest way to create your own case is to **copy an existing tutorial**
and modify it:

```bash
# Copy the LPBF_small tutorial as a starting point
cp -r tutorials/laserbeamFoam/LPBF_small  my_lpbf_case
cd my_lpbf_case
```

For a case **without** a DEM-generated powder bed (e.g., a flat substrate),
the simpler `Plate2D` tutorial is a better starting point:

```bash
cp -r tutorials/laserbeamFoam/Plate2D  my_lpbf_case
cd my_lpbf_case
```

---

## 3. Case Directory Structure

Every `laserbeamFoam` case contains the following directories and files:

```
my_lpbf_case/
├── Allrun               # Script to mesh, set fields, and run the solver
├── Allclean             # Script to remove generated files and reset the case
├── constant/            # Time-independent physical properties
│   ├── LaserProperties      # Laser beam parameters and scan path references
│   ├── transportProperties  # Material thermo-physical properties
│   ├── turbulenceProperties # Turbulence model selection
│   ├── g                    # Gravitational acceleration
│   ├── timeVsLaserPosition  # Laser centre coordinates vs time
│   └── timeVsLaserPower     # Laser power vs time
├── system/              # Mesh and solver configuration
│   ├── blockMeshDict        # Mesh definition
│   ├── controlDict          # Time-stepping and output settings
│   ├── fvSchemes            # Discretisation schemes
│   ├── fvSolution           # Linear-solver settings
│   └── setFieldsDict        # Initial field regions (for flat plate cases)
└── initial/             # Initial field values (copied to 0/ at run time)
    ├── U                    # Velocity
    ├── p_rgh                # Pressure
    ├── T                    # Temperature
    ├── alpha.metal          # Metal volume fraction
    └── TRHS                 # Temperature right-hand side
```

> **Note for powder bed cases:**  `constant/location` contains the DEM particle
> positions; `system/setFieldsDict` is not needed because `setSolidFraction`
> initialises `alpha.metal` from the particle file.

---

## 4. Step 1 – Define the Computational Domain and Mesh

Edit `system/blockMeshDict` to match **your** domain size and desired resolution.

### Key parameters

| Parameter | Description |
|-----------|-------------|
| `vertices` | Coordinates of the eight domain corners |
| `blocks`   | Cell count in each direction `(nx ny nz)` |
| `boundary` | Patch names and types |

### Typical LPBF domain

For a typical LPBF single-track simulation on a flat substrate, the domain
spans ~200 µm in the scan direction, ~150 µm in depth, and ~200 µm in the
transverse direction:

```c++
convertToMeters 1.0;

vertices
(
    (0        0        0      )   // 0
    (200e-6   0        0      )   // 1
    (200e-6   150e-6   0      )   // 2
    (0        150e-6   0      )   // 3
    (0        0        200e-6 )   // 4
    (200e-6   0        200e-6 )   // 5
    (200e-6   150e-6   200e-6 )   // 6
    (0        150e-6   200e-6 )   // 7
);

//  x: scan direction, y: depth, z: transverse
//  Cell count: (nx  ny  nz) — increase for finer resolution
blocks
(
    hex (0 1 2 3 4 5 6 7) (80 60 80) simpleGrading (1 1 1)
);

boundary
(
    leftWall   { type patch; faces ((0 4 7 3)); }
    rightWall  { type patch; faces ((1 2 6 5)); }
    bottomWall { type wall;  faces ((0 3 2 1)); }
    topWall    { type patch; faces ((4 5 6 7)); }
    front      { type patch; faces ((1 5 4 0)); }
    back       { type patch; faces ((3 7 6 2)); }
);
```

### Guidelines

- **Resolution**: Aim for ≥ 5 cells across the expected melt pool width
  (typically 50–150 µm).  Finer meshes improve accuracy but increase
  computation time.
- **Depth**: The domain should be deep enough to contain the build plate plus
  the powder layer (if used).
- **Gas region above powder**: Include at least 50–100 µm of gas above the
  powder surface to allow the vapour plume to develop.

---

## 5. Step 2 – Configure Material Properties

Edit `constant/transportProperties` (also named `phaseProperties` in newer
cases) to specify the thermo-physical properties of **your** metal and the
shielding gas.

```c++
interfaceTrackingScheme isoAdvector;   // or MULES

phases (metal gas);

metal
{
    transportModel  Newtonian;
    nu              5e-07;       // Kinematic viscosity (m²/s)
    rho             8000;        // Density (kg/m³)
    Tsolidus        1658;        // Solidus temperature (K)
    Tliquidus       1723;        // Liquidus temperature (K)
    LatentHeat      2.7e5;       // Latent heat of fusion (J/kg)
    beta            5.0e-6;      // Thermal expansion coefficient (1/K)

    // Polynomial coefficients for thermal conductivity κ(T) (W/m·K)
    // κ = c0 + c1*T + c2*T² + ... (exactly 8 coefficients required)
    // Unused higher-order terms must be set to 0.
    poly_kappa      (10 0.015 0 0 0 0 0 0);

    // Polynomial coefficients for specific heat cp(T) (J/kg·K)
    // cp = c0 + c1*T + c2*T² + ... (exactly 8 coefficients required)
    poly_cp         (520 0.075 0 0 0 0 0 0);
}

gas
{
    transportModel  Newtonian;
    nu              1.48e-05;    // Kinematic viscosity of shielding gas (m²/s)
    rho             1;           // Density of shielding gas (kg/m³)
    Tsolidus        1.0;
    Tliquidus       10.0;
    LatentHeat      1;
    beta            4.0e-5;
    poly_kappa      (0.04 0.0 0 0 0 0 0 0);  // 8 coefficients required
    poly_cp         (520  0.0 0 0 0 0 0 0);   // 8 coefficients required
}

// ---- Phase-change and surface properties ----
sigma               0.07;        // Surface tension (N/m) at the reference temperature
Marangoni_Constant  -0.5e-4;    // dσ/dT (N/m·K)
p0                  100000.0;   // Ambient pressure (Pa)
Tvap                3068.0;     // Vaporisation temperature (K)
Mm                  5.58e-2;    // Molar mass of metal (kg/mol)
LatentHeatVap       7.45e6;     // Latent heat of vaporisation (J/kg)
emS                 0.4;        // Emissivity (solid)
emL                 0.1;        // Emissivity (liquid)
TRef                TRef [0 0 0 1 0 0 0] 300;  // Reference temperature (K)
```

### Common material databases

| Material | `rho` (kg/m³) | `Tsolidus` (K) | `Tliquidus` (K) | `LatentHeat` (J/kg) |
|----------|--------------|----------------|-----------------|---------------------|
| SS316L   | 8000         | 1658           | 1723            | 2.7 × 10⁵           |
| Ti-6Al-4V| 4420         | 1877           | 1933            | 2.86 × 10⁵          |
| AlSi10Mg | 2670         | 833            | 893             | 3.9 × 10⁵           |
| Inconel 625 | 8440      | 1563           | 1623            | 2.27 × 10⁵          |

> Verify these values against the experimental literature for your specific
> alloy grade, as they can vary with composition and processing history.

---

## 6. Step 3 – Configure Laser Parameters

Edit `constant/LaserProperties` to set the beam geometry and optical properties.

```c++
radialPolarHeatSource no;    // 'no' = Gaussian beam; 'yes' = ring/donut beam

// Reference to the position and power schedule files
timeVsLaserPosition
{
    file    "$FOAM_CASE/constant/timeVsLaserPosition";
    outOfBounds clamp;
}

timeVsLaserPower
{
    file    "$FOAM_CASE/constant/timeVsLaserPower";
    outOfBounds clamp;
}

// Incident laser direction (normalised automatically by the solver)
// (0 1 0) = laser beam pointing in the -y direction (downward for typical setup)
V_incident (0 1 0);

// 1/e² beam radius (m) — also called the "spot size" radius
laserRadius 10e-6;           // 10 µm radius (20 µm diameter)

// Number of ray sub-divisions across the beam cross-section
// Higher values = more accurate energy deposition, but slower
N_sub_divisions 1;

// Wavelength of the laser (m) — affects Fresnel absorptivity
wavelength    1.064e-6;      // 1064 nm (Nd:YAG / fibre laser)

// Electron number density of the metal (m⁻³) — used in Fresnel equations
e_num_density 5.83e29;

// Beam radius flavour: 2 = 1/e² radius definition
Radius_Flavour 2.0;

// Set to 'true' when simulating a powder bed; 'false' for flat-plate cases
PowderSim true;
```

### Choosing `laserRadius` and `N_sub_divisions`

- `laserRadius` should be at least **4–5 mesh cells** across to resolve the
  beam accurately.  For a 2.5 µm cell size, a radius of 15–20 µm is
  appropriate.
- Increase `N_sub_divisions` (e.g., to 3–5) for cases with complex surface
  topology (deep keyhole, many reflections) where ray-tracing accuracy is
  critical.

---

## 7. Step 4 – Define the Laser Scan Path and Power Schedule

### `constant/timeVsLaserPower`

A table of `(time  power)` pairs (SI units: seconds and watts).

```
(
    (0        0   )     // laser off at t = 0
    (1e-8   200   )     // ramp up to 200 W very quickly
    (80e-6  200   )     // maintain 200 W until t = 80 µs
    (80.1e-6  0   )     // switch off
    (100e-6   0   )
)
```

### `constant/timeVsLaserPosition`

A table of `(time  (x  y  z))` pairs giving the laser beam centre position
in metres.

```
(
    //  time       ( x-centre    y-centre     z-centre )
    (0             (100e-6    -200e-6    100e-6))   // start position (above domain)
    (80e-6         (100e-6    -200e-6    180e-6))   // end position after scan
    (100e-6        (100e-6    -200e-6    180e-6))
)
```

> **Coordinate convention (after `transformPoints`):**  
> In the tutorial cases, the mesh is rotated so that the laser travels in the
> **z**-direction and the build direction is **y** (downward into the substrate).
> The y-coordinate of the laser centre is set to a value outside the domain
> (e.g., `-200e-6`) because the solver projects the beam onto the free surface.

### Calculating scan speed

Given a track length `L` and a laser-on time `t_on`, the scan speed is simply:

```
v_scan = L / t_on
```

For example, L = 80 µm scanned in 80 µs gives v = 1 m/s.

---

## 8. Step 5 – Set Up Initial Conditions

The `initial/` directory contains the starting values for all solution fields.
These files are copied to `0/` before each run.

### Temperature (`initial/T`)

Set the preheat (ambient) temperature, e.g., 300 K:

```c++
dimensions      [0 0 0 1 0 0 0];   // Kelvin
internalField   uniform 300.0;

boundaryField
{
    bottomWall { type  fixedValue; value uniform 300.0; }  // fixed base temperature
    topWall    { type  zeroGradient; }
    leftWall   { type  zeroGradient; }
    rightWall  { type  zeroGradient; }
    front      { type  zeroGradient; }
    back       { type  zeroGradient; }
}
```

### Metal volume fraction (`initial/alpha.metal`)

For a **flat-plate** case, `alpha.metal` is set to 0 everywhere initially and
then the substrate region is filled by `setFields` using `system/setFieldsDict`.
For a **powder bed** case, `setSolidFraction` fills the field from DEM particle
positions.

```c++
// initial/alpha.metal
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 0;          // start with all gas

boundaryField
{
    back       { type zeroGradient; }
    front      { type zeroGradient; }
    leftWall   { type zeroGradient; }
    rightWall  { type zeroGradient; }
    topWall    { type zeroGradient; }
    bottomWall { type zeroGradient; }
}
```

### `system/setFieldsDict` (flat plate only)

Define a box region that fills the lower portion of the domain with solid metal:

```c++
defaultFieldValues
(
    volScalarFieldValue alpha.metal 0
);

regions
(
    boxToCell
    {
        // Fill the bottom 50 µm with solid metal (substrate)
        // box (xMin yMin zMin) (xMax yMax zMax)
        box (-1 0.0 -1) (1 50e-6 1);
        fieldValues
        (
            volScalarFieldValue alpha.metal 1
        );
    }
);
```

Adjust the `box` coordinates to match your substrate thickness.

---

## 9. Step 6 – Tune Solver Numerics

### `system/controlDict`

```c++
application     laserbeamFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         100e-6;       // Total simulation time (s)

deltaT          1.0e-8;       // Initial time step (s)

writeControl    adjustableRunTime;
writeInterval   5e-6;         // Output every 5 µs

adjustTimeStep  yes;
maxCo           0.2;          // Maximum Courant number (keep ≤ 0.2)
maxAlphaCo      0.2;
maxDeltaT       5e-5;
```

### `system/fvSchemes`

The default schemes from the tutorials are appropriate for most LPBF cases:

```c++
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U)           Gauss linearUpwind grad(U);
    div(phi,alpha)          Gauss interfaceCompression vanLeer 1;
    div(phirb,alpha)        Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
    div(rhophicp,T)         Gauss upwind;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
```

### `system/fvSolution`

Key settings to adjust are the PIMPLE outer correctors and the melting loop:

```c++
MELTING
{
    minTempCorrector  1;
    maxTempCorrector  20;
    epsilonTolerance  1e-6;
    epsilonRelaxation 0.9;
    damperSwitch      true;
}

PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
}
```

Increasing `nOuterCorrectors` (e.g., to 2–3) can improve stability for high
power densities or deep keyhole cases.

---

## 10. Step 7 – Run the Simulation

### Create or update your `Allrun` script

For a flat-plate (no DEM) case:

```bash
#!/bin/bash
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# 1. Initialise field directory
cp -r initial 0

# 2. Generate mesh
runApplication blockMesh

# 3. Set initial metal volume fraction
runApplication setFields

# 4. Run solver in serial
runApplication laserbeamFoam
```

For a **parallel** run (e.g., on 8 cores), replace the last two lines with:

```bash
runApplication decomposePar
runParallel laserbeamFoam
runApplication reconstructPar
```

### Running the case

```bash
chmod +x Allrun
./Allrun
```

Progress is written to `log.blockMesh`, `log.setFields`, and
`log.laserbeamFoam`.  Monitor the residuals in `log.laserbeamFoam`:

```bash
tail -f log.laserbeamFoam
```

### Estimated run times

| Case type      | Mesh (cells) | Typical run time |
|----------------|-------------|-----------------|
| 2-D plate      | ~20 k       | Minutes (serial) |
| 3-D flat plate | ~150 k      | Hours (serial) / Minutes (parallel) |
| 3-D powder bed | ~200–400 k  | Several hours (parallel) |

---

## 11. Step 8 – Post-Process and Visualise Results

### ParaView

Open the case in ParaView:

```bash
paraFoam -builtin &
```

Or generate the `.foam` reader file first:

```bash
touch myCase.foam
paraview myCase.foam &
```

**Useful fields to visualise:**

| Field          | Description |
|----------------|-------------|
| `alpha.metal`  | Metal volume fraction (0 = gas, 1 = solid metal) |
| `T`            | Temperature (K) — identifies melt pool extent |
| `U`            | Velocity (m/s) — shows Marangoni flow patterns |
| `p_rgh`        | Pressure — shows recoil pressure distribution |

**Visualising the melt pool boundary:**  
Apply a `Contour` filter on `T` at the liquidus temperature
(e.g., 1723 K for SS316L) to see the melt pool boundary.

### Visualising laser rays

The solver writes ray VTK files to `VTK/rays_<LASER_NAME>_<TIME>.vtk`.
Open them in ParaView:

1. `File` → `Open...` → select `VTK/rays_<LASER_NAME>.vtk.series`
   (this file is written at the end of the simulation and loads all time steps
   at the correct physical times automatically).
2. Set the colour to black and increase the `Line Width`, or apply the
   `Tube` filter for better visibility.

---

## 12. Optional: Powder Bed Setup with DEM

For a realistic LPBF simulation, a random powder bed can be generated using
[LIGGGHTS®](https://www.cfdem.com) (see [README.md](../README.md) for
installation instructions).

### Workflow

1. **Edit `DEM_small/input.liggghts`** to set your particle size, size
   distribution, and number of particles.
2. **Run the DEM simulation:**
   ```bash
   cd DEM_small && ./Allrun
   ```
3. **Copy the particle positions** to the constant directory:
   ```bash
   cp DEM_small/post/location constant/location
   ```
4. **Run `setSolidFraction`** in the CFD case directory to initialise
   `alpha.metal` from the particle file:
   ```bash
   blockMesh
   setSolidFraction -subDivisions 6
   ```
5. Set `PowderSim true;` in `constant/LaserProperties`.

For detailed DEM instructions, see the
[LPBF_small tutorial README](../tutorials/laserbeamFoam/LPBF_small/README.md).

---

## 13. Common Issues and Tips

### Solver diverges / very small time steps

- Reduce `maxCo` from 0.2 to 0.1.
- Check that `laserRadius` spans at least 4–5 cells.
- Check that the laser position table places the beam **above** the metal
  surface at `t = 0` (the beam is projected onto the surface automatically).
- Ensure material properties (particularly `Tsolidus` and `Tliquidus`) are
  physically consistent.

### Melt pool does not form

- Verify that `timeVsLaserPower` provides non-zero power at the correct times.
- Check that `V_incident` points from the laser source **toward** the surface.
- Confirm that the laser beam centre (from `timeVsLaserPosition`) passes over
  the metal surface during the simulation.

### Poor mass conservation

- Increase `nAlphaSubCycles` in `fvSolution` from 1 to 2 or 3.
- Use `isoAdvector` as `interfaceTrackingScheme` (generally more accurate than
  MULES for 3-D cases).

### Customising the laser for a ring/donut beam

Set `radialPolarHeatSource yes;` in `constant/LaserProperties` to activate a
ring-shaped heat source profile.

### Running on a high-performance cluster

Generate a `system/decomposeParDict` (e.g., using `scotch` decomposition):

```c++
numberOfSubdomains 32;
method scotch;
```

Then run:

```bash
decomposePar
mpirun -np 32 laserbeamFoam -parallel > log.laserbeamFoam 2>&1
reconstructPar
```

---

## Summary Checklist

Use this checklist when setting up a new LPBF case:

- [ ] Copy an existing tutorial case as a starting point
- [ ] Update `system/blockMeshDict` – domain size and mesh resolution
- [ ] Update `constant/transportProperties` – material properties for your alloy
- [ ] Update `constant/LaserProperties` – beam radius, wavelength, `PowderSim`
- [ ] Update `constant/timeVsLaserPower` – laser power schedule
- [ ] Update `constant/timeVsLaserPosition` – laser scan path
- [ ] Update `initial/T` – initial (preheat) temperature
- [ ] Update `system/setFieldsDict` (flat plate) or run DEM + `setSolidFraction` (powder bed)
- [ ] Update `system/controlDict` – `endTime`, `writeInterval`
- [ ] Run `./Allrun` and monitor `log.laserbeamFoam`
- [ ] Visualise results in ParaView

---

*For further background on the solver physics and numerical methods, refer to
the [PDF documentation](LaserbeamFoam___V2.pdf) and the peer-reviewed
publications listed in the main [README.md](../README.md).*
