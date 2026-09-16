# Ubuntu 20.04 Development Handoff

## Purpose

This document transfers the repository understanding established during the initial source audit to a new development session on Ubuntu 20.04. It records confirmed implementation facts, migration risks, the current workspace state, and a safe order for resuming work.

It is a point-in-time handoff, not a replacement for repository rules or source inspection. If this document conflicts with current source code or CMake configuration, the current source is authoritative and the discrepancy should be reported.

## Snapshot

- Handoff date: 2026-09-16
- Audited workspace: `E:\trajopt`
- Git branch at handoff: `main`
- Git commit at handoff: `eafe7d8e9f5a5c5fced902b28cbb810cb2e4d7e1`
- Target migration environment: Ubuntu 20.04
- Work performed: source and documentation audit only; no dependency installation, large experiment, or real-robot execution was performed.

At handoff time, `git status --short` reported the following untracked paths:

```text
?? AGENTS.md
?? RBF_LIO.pdf
?? README.md
?? docs/
```

These files will not be transferred by a normal `git clone` until they are committed. If the migration uses a filesystem copy instead, preserve them explicitly and verify their hashes or contents on the target machine.

## Required Reading Order

Before modifying or building the project, read:

1. [AGENTS.md](../AGENTS.md) — repository-wide rules, invariants, build facts, and protected files.
2. This handoff — migration state, confirmed blockers, and recommended first steps.
3. [README.md](../README.md) — concise project entry and smoke workflow.
4. [DEPLOYMENT.md](DEPLOYMENT.md) — Ubuntu dependencies, build, offline execution, and deployment gaps.
5. [ARCHITECTURE.md](ARCHITECTURE.md) — module boundaries, dependencies, data flow, and data contracts.
6. [ALGORITHM.md](ALGORITHM.md) — mathematical formulation, decision variables, objectives, constraints, and numerical limitations.
7. [EXPERIMENT.md](EXPERIMENT.md) — reproducibility protocol, artifact status, metrics, and proposed experiments.
8. [RBF_LIO.pdf](../RBF_LIO.pdf) — upstream RBF terrain representation reference.

Always cross-check important claims against the current entry points and `CMakeLists.txt` files. Historical `build*` directories, logs, PNG files, and result directories may come from older code.

## Project Scope

The repository implements offline, planning-side trajectory generation for a two-wheel / wheel-legged mechanism on a continuous 2.5D RBF height field.

RBF-LIO is upstream. It may produce the RBF centers and weights consumed here, but this repository does **not** implement:

- RBF-LIO odometry or SLAM;
- LiDAR/IMU fusion;
- online RBF center selection or fitting;
- a robot state estimator;
- ROS communication;
- a low-level controller;
- real-time trajectory execution or replanning.

The primary implementation is C++/Ipopt. Root Python/CasADi code is prototype work and is not equivalent to the C++ pipeline.

## Confirmed Main Pipeline

```text
RBF JSON
  -> terrain_friction_segment
  -> terrain_segmentation_data.txt
  -> rrt_star_region
  -> waypoints_segmented.txt
  -> seg_trajopt_cpp for a rollable segment
     or lift_leg_trajopt_cpp for a lift segment
  -> Ipopt
  -> initial_trajectory.txt / final_trajectory.txt
```

The region-label contract is fixed:

```text
0 = lift / non-rollable under the current friction proxy
1 = rollable
```

Do not change this convention in one producer or consumer without migrating every reader, writer, fixture, and document.

This pipeline is only partially orchestrated. There is no executable that automatically traverses every labeled segment, dispatches the appropriate optimizer, enforces adjacent-segment continuity, stitches the results, time-parameterizes the route, and exports a robot command.

## Main Components

| Component | Confirmed role |
|---|---|
| `terrain_friction_segment` | Reads RBF JSON, evaluates a `512 x 512` terrain grid with CUDA, applies a slope/friction feasibility test, and writes segmentation artifacts. |
| `rrt_star_region` | Runs a two-dimensional region-cost RRT*. Lift regions remain traversable but have higher edge cost. |
| `trajopt_cpp` | Standalone general full-trajectory Ipopt NLP with 10 Chebyshev trajectory blocks and 6 sampled contact-force blocks. It does not consume segmented waypoints. |
| `seg_trajopt_cpp` | Optimizes the first contiguous label-`1` waypoint segment. Wheel `z` is evaluated from the RBF terrain. |
| `lift_leg_trajopt_cpp` | Optimizes the first contiguous label-`0` segment, with piecewise wheel-height coefficients and swing/contact initialization. It falls back to the first rollable segment if no lift segment exists. |
| `core/`, `constraints/`, `geometry/`, `solver/` | Modular Python prototype abstractions. They are not assembled into the C++ mainline. |
| `solve_casadi.py`, `trajopt_sdf.py` | Standalone Python/CasADi research prototypes with hard-coded inputs. |
| `test/terrain_res/` | Smoke-test RBF, segmentation, RRT path, and labeled waypoint artifacts. Their generating commit, command, and environment are not fully recorded. |

`seg_trajopt_cpp` and `lift_leg_trajopt_cpp` compile shared Chebyshev, RBF, JSON-reader, polynomial, and yaw-objective sources directly from the sibling `trajopt_cpp` directory. Preserve this directory relationship.

## RBF Data Contract

The current C++ readers require:

```json
{
  "cur_rbf_grid": [[x0, y0], [x1, y1]],
  "rbf_weight": [w0, w1]
}
```

Confirmed behavior:

- centers and weights must have the same length and matching index order;
- the terrain is a single-valued height field `z = f(x,y)`;
- the Gaussian kernel uses the denominator `2 * sigma^2`;
- C++ `JSONReader` currently hard-codes `sigma = 0.14` rather than reading it from JSON;
- the JSON schema does not define frame ID, units, timestamp, valid region, schema version, or sigma provenance.

The repository documentation interprets positions as metres, yaw as radians, force as newtons, mass as kilograms, and gravity as metres per second squared, but the input JSON itself does not encode these units. This must become explicit before real-robot use.

## Current CLI Contracts

```text
terrain_friction_segment [rbf.json] [output_dir] [mu]

rrt_star_region [data.txt] [start_x start_y goal_x goal_y] [path_out.txt]

trajopt_cpp [rbf.json]

seg_trajopt_cpp [data_dir] [output_dir]

lift_leg_trajopt_cpp [data_dir] [output_dir]
```

Important behavior:

- terrain segmentation defaults to `mu = 0.3` and clamps a supplied value to `[0.1, 1.0]`;
- `seg_trajopt_cpp` expects `rbf.json` and `waypoints_segmented.txt` in `data_dir`, and attempts to read `terrain_segmentation_data.txt` there;
- `lift_leg_trajopt_cpp` expects `rbf.json` and `waypoints_segmented.txt` in `data_dir`;
- segment output files have columns `s xb yb zb psi xl yl zl xr yr zr`;
- `s` is a normalized geometric path parameter, not physical time;
- `final_trajectory.txt` may contain the final Ipopt iterate even when convergence is unsuccessful, so always inspect the return status and logs.

## Ubuntu 20.04 Dependencies

Confirmed project dependencies are:

- C++17 compiler;
- CMake 3.10+ for the three Ipopt projects;
- CMake 3.18+ for the CUDA project;
- Ipopt development library;
- MUMPS, selected in the current entry points with `linear_solver=mumps`;
- `pkg-config` as one Ipopt discovery path;
- CUDA Toolkit and a compatible NVIDIA driver/GPU for terrain segmentation;
- Python 3, NumPy, and Matplotlib for visualization;
- CasADi only for the Python NLP prototypes;
- CuPy only for an optional performance test.

The current CMake files explicitly mention the Ubuntu package `coinor-libipopt-dev`. `libmumps-dev` is also relevant to the configured linear solver. Do not assume a CUDA version or CuPy package variant: select them from the actual GPU, driver, Ubuntu package sources, and CUDA compatibility matrix.

No dependency lock file, container definition, unified installer, or top-level superbuild exists.

## Recommended Fresh Build Layout

Do not reuse or modify copied `build`, `build_lift`, `build_lift_new`, or other historical caches. From the repository root, the proposed fresh layout is:

```bash
cmake -S trajopt_cpp -B _build_ubuntu2004/trajopt_cpp \
  -DCMAKE_BUILD_TYPE=Release
cmake --build _build_ubuntu2004/trajopt_cpp --parallel

cmake -S seg_trajopt_cpp -B _build_ubuntu2004/seg_trajopt_cpp \
  -DCMAKE_BUILD_TYPE=Release
cmake --build _build_ubuntu2004/seg_trajopt_cpp --parallel

cmake -S lift_leg_trajopt_cpp -B _build_ubuntu2004/lift_leg_trajopt_cpp \
  -DCMAKE_BUILD_TYPE=Release
cmake --build _build_ubuntu2004/lift_leg_trajopt_cpp --parallel

cmake -S terrain_friction_segment \
  -B _build_ubuntu2004/terrain_friction_segment \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=<target_arch>
cmake --build _build_ubuntu2004/terrain_friction_segment --parallel
```

These commands were derived from current CMake files but were not executed during the source-audit conversation. Determine `<target_arch>` from the actual target GPU. The CUDA CMake default is `70;80;86`; it may not match the Ubuntu machine.

A useful build order is the standalone CPU optimizer first, then both segment optimizers, then the CUDA project. This separates Ipopt/MUMPS problems from CUDA toolchain problems.

## Recommended Environment Audit

Before installation or build, collect:

```bash
lsb_release -a
uname -m
git status --short
git rev-parse HEAD
g++ --version
cmake --version
pkg-config --modversion ipopt
nvidia-smi
nvcc --version
python3 --version
python3 -c "import numpy, matplotlib; print(numpy.__version__, matplotlib.__version__)"
```

Missing tools are expected on a fresh machine. Record them first; do not silently modify the system environment or select dependency versions without user approval.

## Safe Smoke-Test Strategy

After successful builds, use a new output directory and do not overwrite `test/terrain_res`:

```bash
RUN=/tmp/trajopt_ubuntu_smoke
mkdir -p "$RUN/data" "$RUN/roll" "$RUN/lift" "$RUN/logs"

MPLBACKEND=Agg \
./_build_ubuntu2004/terrain_friction_segment/terrain_friction_segment \
  test/terrain_res/rbf.json "$RUN/data" 0.3 \
  >"$RUN/logs/segmentation.log" 2>&1

cp test/terrain_res/waypoints_segmented.txt "$RUN/data/waypoints_segmented.txt"

./_build_ubuntu2004/seg_trajopt_cpp/seg_trajopt_cpp \
  "$RUN/data" "$RUN/roll" \
  >"$RUN/logs/roll.log" 2>&1

./_build_ubuntu2004/lift_leg_trajopt_cpp/lift_leg_trajopt_cpp \
  "$RUN/data" "$RUN/lift" \
  >"$RUN/logs/lift.log" 2>&1

python3 seg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  "$RUN/roll" --no-show

python3 lift_leg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  "$RUN/lift" --rbf "$RUN/data/rbf.json" --no-show
```

This is a smoke workflow, not a fully regenerated benchmark: it reuses the checked-in labeled waypoints because the current RRT* output routing is not portable. Review Ipopt statuses and logs rather than treating output-file existence as success.

## Confirmed Migration and Development Blockers

### Critical

1. `rrt_star_region` only allows the path text output to be redirected. It still writes `waypoints_segmented.txt` and `rrt_star_path.png` to `/home/yizhe/trajopt/test/terrain_res` and attempts to open the PNG.
2. There is no whole-route segment orchestrator, continuity enforcement, or trajectory stitcher.
3. The generated trajectory has no physical timing, velocity/acceleration contract, or robot-controller message format.
4. No ROS publisher/subscriber, state-feedback loop, command execution layer, safety supervisor, or real-time replanning interface is implemented.
5. Coordinate-frame, timestamp, unit, terrain-snapshot, and start-state contracts are incomplete for robot deployment.

### Important

1. RBF sigma, start/goal values, robot geometry, Ipopt options, objective weights, and many thresholds are hard-coded.
2. The general `trajopt_cpp` and the roll/lift segmented pipeline are parallel implementations, not stages automatically assembled into one solver.
3. Constraints are checked at finite sample/collocation points; continuous-interval feasibility is not guaranteed.
4. Some contact, force-balance, torque-balance, and friction behavior in the general optimizer is represented by soft objectives.
5. Python and C++ models have diverged and should not be assumed numerically equivalent.
6. Existing test outputs have uncertain provenance and are unsuitable as performance claims without regeneration records.

### Convenience

1. There is no unified build/run command, dependency lock, container, or install target.
2. Some programs launch plotting or image-opening commands automatically, complicating SSH and headless use.
3. Output files lack schema version and run metadata.

## Current Documentation Inconsistency

The four detailed documents currently reside under `docs/`:

```text
docs/ARCHITECTURE.md
docs/ALGORITHM.md
docs/EXPERIMENT.md
docs/DEPLOYMENT.md
```

At handoff time, `README.md` still linked to root-level `ARCHITECTURE.md`, `ALGORITHM.md`, `EXPERIMENT.md`, and `DEPLOYMENT.md`. Those links are therefore stale. This handoff task intentionally did not modify `README.md`; repair the links only in a separately authorized change.

## Recommended First Development Tasks

After the Ubuntu environment and baseline builds are verified, the highest-value sequence is:

1. Make every RRT* output path configurable and add a headless/no-open option.
2. Add a single configuration source for RBF sigma, start/goal, robot geometry, solver settings, and objective weights.
3. Add a route orchestrator that enumerates every contiguous labeled segment.
4. Define and enforce segment-boundary state, derivative, yaw, wheel, and contact continuity.
5. Stitch segment outputs and introduce physical time parameterization.
6. Define versioned coordinate-frame, unit, timestamp, and trajectory schemas.
7. Add finite-difference checks for analytic gradients/Jacobians and small deterministic regression tests.
8. Only after the offline contracts are stable, add a trajectory adapter and robot/ROS integration layer.

Do not combine all of these into one initial patch. Establish a reproducible Ubuntu baseline first and preserve the checked-in fixtures and historical outputs.

## Updated Prompt for the Next Conversation

Copy the following prompt into the Ubuntu 20.04 development conversation:

```text
I have migrated the trajopt repository to Ubuntu 20.04. Start with a read-only migration audit; do not modify source code, CMakeLists.txt, test fixtures, existing outputs, or the system environment, and do not install packages until I approve.

Read these files in order:

1. AGENTS.md
2. docs/HANDOFF.md
3. README.md
4. docs/DEPLOYMENT.md
5. docs/ARCHITECTURE.md
6. docs/ALGORITHM.md
7. docs/EXPERIMENT.md

Then cross-check all important facts against the current source entry points and CMakeLists.txt files. Treat current source and CMake as authoritative. Do not reuse copied build*, build_lift*, logs, images, or experiment outputs.

Known project facts to verify:

- This repository is downstream of RBF-LIO; it does not implement RBF-LIO, SLAM, or LiDAR/IMU state estimation.
- The C++/Ipopt code is the mainline; Python/CasADi code is prototype work.
- Region labels are fixed as 0=lift and 1=rollable.
- seg_trajopt_cpp handles only the first rollable segment.
- lift_leg_trajopt_cpp handles only the first lift segment and can fall back to the first rollable segment.
- There is no all-segment orchestration, stitching, time parameterization, ROS interface, or robot controller interface.
- rrt_star_region contains hard-coded /home/yizhe/trajopt output paths for labeled waypoints and PNG output.
- C++ JSONReader currently hard-codes RBF sigma to 0.14.

Perform and report:

1. Confirm Ubuntu release, CPU architecture, repository path, branch, commit, and git status.
2. Confirm that all handoff documents were transferred. Report any untracked files or broken README links.
3. Inventory g++, CMake, pkg-config, Ipopt, MUMPS, CUDA, NVIDIA driver/GPU, Python, NumPy, and Matplotlib. Also report CasADi/CuPy separately as optional prototype dependencies.
4. Determine whether the GPU compute capability matches the CMake default CUDA architectures 70;80;86.
5. List missing dependencies, but do not run apt, pip, CUDA installers, or environment changes without approval.
6. Propose fresh _build_ubuntu2004/ configure and build commands for trajopt_cpp, seg_trajopt_cpp, lift_leg_trajopt_cpp, and terrain_friction_segment.
7. Propose a smoke test using test/terrain_res and a new /tmp output directory. Do not overwrite repository artifacts or run a long Ipopt/CUDA experiment yet.
8. Separate confirmed blockers into Critical, Important, and Convenience.
9. Mark anything that cannot be confirmed from source as UNKNOWN; do not infer standard ROS or robot capabilities.

Finish with: environment compatibility, available/missing dependencies, recommended commands, hard-coded path problems, files that would eventually need changes, and the maximum currently achievable offline workflow without source modification. Wait for my approval before installation, compilation, or code changes.
```

## Handoff Completion Criteria

The Ubuntu migration is ready for active development only after:

- the handoff documents are present and tracked or otherwise preserved;
- the target commit and local modifications are recorded;
- fresh build directories configure successfully;
- Ipopt reports MUMPS availability;
- CUDA architecture and driver/toolkit compatibility are confirmed;
- a non-destructive smoke run writes outputs to a new directory;
- solver statuses and generated file formats are checked;
- no historical fixture or experiment output has been overwritten.
