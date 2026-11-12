#!/usr/bin/env python3
"""Capture faithful simulator replays and render optional 3D MP4 videos.

This is a presentation tool, not a benchmark. See --help and the README.
"""

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np


ROOT = Path(__file__).resolve().parent
PALETTE = np.array([
    [83, 205, 255], [255, 184, 94], [145, 225, 172], [211, 156, 255],
    [255, 132, 157], [241, 223, 125], [123, 178, 255], [128, 222, 212],
], dtype=np.uint8)


def positive_int(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def nonnegative_int(value):
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return result


def finite_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("must be finite")
    return result


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    for command in ("capture", "render", "movie"):
        sub = commands.add_parser(command)
        sub.add_argument("input", type=Path, help="trajectory .npz" if command == "render"
                         else "simulator instance file")
        if command != "render":
            sub.add_argument("--app", type=Path, default=ROOT / "build/app")
            sub.add_argument("--solver", choices=("n2", "barnes-hut", "fmm"), default="n2")
            config = sub.add_mutually_exclusive_group()
            config.add_argument("--theta", type=finite_float)
            config.add_argument("--order", type=int, choices=range(1, 7))
            sub.add_argument("--dt", type=finite_float, required=True)
            sub.add_argument("--steps", type=nonnegative_int, required=True)
            sub.add_argument("--every", type=positive_int, default=10)
            sub.add_argument("--output-dir", type=Path)
        if command != "capture":
            sub.add_argument("--fps", type=positive_int, default=30)
            sub.add_argument("--width", type=positive_int, default=1920)
            sub.add_argument("--height", type=positive_int, default=1080)
            sub.add_argument("--trail-frames", type=nonnegative_int, default=60)
            sub.add_argument("--camera", choices=("fixed", "orbit"), default="fixed")
        if command == "render":
            sub.add_argument("--output", type=Path, help="default: simulation.mp4 beside trajectory")
    return result


def validate_options(args):
    if args.command != "render":
        if args.dt <= 0 or not math.isfinite(args.dt * args.steps):
            raise ValueError("dt must be positive and total simulation time must be finite")
        if args.theta is not None and (args.solver != "barnes-hut" or args.theta < 0):
            raise ValueError("--theta requires barnes-hut and a nonnegative value")
        if args.order is not None and args.solver != "fmm":
            raise ValueError("--order requires fmm")
    if args.command != "capture" and (args.width % 2 or args.height % 2):
        raise ValueError("video width and height must be even for H.264")


def frame_steps(steps, every):
    if steps < 0 or every <= 0:
        raise ValueError("steps must be nonnegative and every must be positive")
    values = list(range(0, steps + 1, every))
    if values[-1] != steps:
        values.append(steps)
    return np.array(values, dtype=np.int64)


def parse_state(text):
    tokens = text.split()
    if not tokens or not tokens[0].isascii() or not tokens[0].isdecimal():
        raise ValueError("state must start with a body count")
    count = int(tokens[0])
    if count == 0:
        raise ValueError("a visualization needs at least one body")
    if len(tokens) != 1 + 10 * count:
        raise ValueError("state does not contain exactly ten fields per body")
    values = np.array(tokens[1:], dtype=np.float64).reshape(count, 10)
    if not np.isfinite(values).all() or not (values[:, 0] > 0).all():
        raise ValueError("state fields must be finite and masses must be positive")
    return values


def validate_trajectory(data):
    positions, masses = data["positions"], data["masses"]
    steps, times, metadata = data["steps"], data["times"], data["metadata"]
    if (positions.ndim != 3 or positions.shape[2] != 3 or positions.shape[0] < 1
            or masses.ndim != 1 or len(masses) < 1 or positions.shape[1] != len(masses)
            or steps.shape != (len(positions),) or times.shape != steps.shape):
        raise ValueError("invalid trajectory array shapes")
    if (positions.dtype != np.float64 or masses.dtype != np.float64
            or times.dtype != np.float64 or steps.dtype != np.int64):
        raise ValueError("trajectory requires float64 coordinates/masses/times and int64 steps")
    if (not all(np.isfinite(a).all() for a in (positions, masses, times))
            or not (masses > 0).all() or steps[0] != 0
            or np.any(steps[1:] <= steps[:-1])):
        raise ValueError("invalid trajectory values or step sequence")
    if (not isinstance(metadata, dict) or metadata.get("version") != 1
            or metadata.get("solver") not in ("n2", "barnes-hut", "fmm")):
        raise ValueError("unsupported trajectory metadata")
    dt = metadata.get("dt")
    if (not isinstance(dt, (int, float)) or not math.isfinite(dt) or dt <= 0
            or not np.array_equal(times, steps * dt)):
        raise ValueError("trajectory times must equal steps multiplied by dt")


def save_trajectory(path, data):
    validate_trajectory(data)
    with Path(path).open("xb") as output:
        np.savez_compressed(output, positions=data["positions"], masses=data["masses"],
                            steps=data["steps"], times=data["times"],
                            metadata_json=np.array(json.dumps(data["metadata"], allow_nan=False)))


def load_trajectory(path):
    with np.load(path, allow_pickle=False) as source:
        data = {key: source[key] for key in ("positions", "masses", "steps", "times")}
        data["metadata"] = json.loads(source["metadata_json"].item())
    validate_trajectory(data)
    return data


def publish(source, destination):
    """Publish a completed file atomically, refusing existing destinations."""
    os.link(source, destination)


def capture(args):
    validate_options(args)
    source, executable = args.input.resolve(), args.app.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("simulator executable is missing or not executable; build it with make first")
    # Validate before reserving an output directory; parse the frozen copy again below.
    parse_state(source.read_text())
    output = args.output_dir or ROOT / "results" / (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + "_visualization")
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    schedule = frame_steps(args.steps, args.every)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    with TemporaryDirectory(prefix=".capture-", dir=output) as temporary:
        work = Path(temporary)
        frozen_input, frozen_app = work / "input.txt", work / "app"
        shutil.copyfile(source, frozen_input)
        shutil.copy2(executable, frozen_app)
        original = parse_state(frozen_input.read_text())
        positions = np.empty((len(schedule), len(original), 3), dtype=np.float64)
        positions[0] = original[:, 1:4]
        final_text = frozen_input.read_text()
        command = [str(frozen_app), "--solver", args.solver]
        if args.theta is not None:
            command += ["--theta", str(args.theta)]
        if args.order is not None:
            command += ["--order", str(args.order)]
        with (output / "capture.stderr.log").open("x") as log:
            for index, step in enumerate(schedule[1:], 1):
                print("Capture {}/{}: step {}".format(index + 1, len(schedule), step), flush=True)
                try:
                    completed = subprocess.run(
                        command + [str(frozen_input), str(args.dt), str(int(step))],
                        cwd=work, env=environment, capture_output=True, text=True, check=False)
                except OSError as error:
                    raise ValueError("capture frame {}, step {}: {}".format(index, step, error)) from error
                log.write("Frame {}, step {}\n{}\n".format(index, step, completed.stderr))
                log.flush()
                try:
                    if completed.returncode:
                        raise ValueError("simulator exited with status {}".format(completed.returncode))
                    state = parse_state(completed.stdout)
                    if state.shape != original.shape or not np.array_equal(state[:, 0], original[:, 0]):
                        raise ValueError("simulator changed body count or masses")
                except ValueError as error:
                    raise ValueError("capture frame {}, step {}: {}\n{}".format(
                        index, step, error, completed.stderr.strip())) from error
                positions[index] = state[:, 1:4]
                final_text = completed.stdout
        metadata = {
            "version": 1, "capture": "faithful-replay", "solver": args.solver,
            "theta": (args.theta if args.theta is not None else 0.5) if args.solver == "barnes-hut" else None,
            "order": (args.order if args.order is not None else 4) if args.solver == "fmm" else None,
            "dt": args.dt, "steps": args.steps, "every": args.every,
            "threads": 1, "input": str(source), "app": str(executable),
        }
        data = {"positions": positions, "masses": original[:, 0].copy(),
                "steps": schedule, "times": schedule * args.dt, "metadata": metadata}
        save_trajectory(work / "trajectory.npz", data)
        (work / "final.txt").write_text(final_text)
        for name in ("input.txt", "final.txt", "trajectory.npz"):
            publish(work / name, output / name)
    print("Trajectory: {}".format(output / "trajectory.npz"), flush=True)
    return output / "trajectory.npz"


def scene_bounds(positions):
    low = positions.min(axis=(0, 1))
    high = positions.max(axis=(0, 1))
    center = low / 2 + high / 2
    span = float(np.max(high - low))
    half = (span if span > 0 else 1.0) * 0.575
    if not np.isfinite(center).all() or not math.isfinite(half):
        raise ValueError("trajectory bounds are too large to render")
    return center, half


def trail_mesh(pv, positions, frame, history, bodies, colors):
    start = max(0, frame - history)
    segments = frame - start
    if segments == 0:
        return None
    points = positions[start:frame + 1, bodies, :].transpose(1, 0, 2).reshape(-1, 3)
    indices = np.arange(len(bodies))[:, None] * (segments + 1) + np.arange(segments)
    lines = np.column_stack((np.full(indices.size, 2), indices.ravel(), indices.ravel() + 1))
    mesh = pv.PolyData(points, lines=lines.ravel())
    rgba = np.empty((len(bodies), segments, 4), dtype=np.uint8)
    rgba[:, :, :3] = colors[bodies, None, :]
    rgba[:, :, 3] = np.linspace(15, 210, segments).astype(np.uint8)
    mesh.cell_data["rgba"] = rgba.reshape(-1, 4)
    return mesh


def verify_video(path, frames, width, height):
    import imageio_ffmpeg

    # Decode every frame in FFmpeg, returning only dimensions and per-frame checksums.
    # This avoids buffering raw movie pixels or leaving decoder pipe handles open.
    completed = subprocess.run(
        [imageio_ffmpeg.get_ffmpeg_exe(), "-v", "error", "-xerror", "-i", str(path),
         "-map", "0:v:0", "-fps_mode", "passthrough", "-f", "framehash", "-hash", "adler32", "-"],
        capture_output=True, text=True, check=False)
    if completed.returncode:
        raise ValueError("video decoding failed: {}".format(completed.stderr.strip()))
    lines = completed.stdout.splitlines()
    if "#dimensions 0: {}x{}".format(width, height) not in lines:
        raise ValueError("encoded video has unexpected dimensions")
    count = sum(bool(line.strip()) and not line.startswith("#") for line in lines)
    if count != frames:
        raise ValueError("encoded video has {} frames; expected {}".format(count, frames))


def render(trajectory, output, args):
    validate_options(args)
    # Keep these optional libraries out of capture and numerical test execution.
    import pyvista as pv

    if sys.platform == "darwin" and not pv.system_supports_plotting():
        raise RuntimeError("the Mac graphics system is unavailable to this process; "
                           "run rendering from a desktop session with graphics access")
    data = load_trajectory(trajectory)
    output = output.resolve()
    if output.suffix.lower() != ".mp4":
        raise ValueError("video output must have an .mp4 extension")
    if output.exists():
        raise FileExistsError("refusing to overwrite {}".format(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    positions = data["positions"]
    center, half = scene_bounds(positions)
    # Translating by one fixed origin improves display precision without changing motion.
    positions = positions - center
    count = len(data["masses"])
    colors = PALETTE[np.arange(count) % len(PALETTE)]
    bodies = np.linspace(0, count - 1, min(count, 16), dtype=int)
    with TemporaryDirectory(prefix=".render-", dir=output.parent) as temporary:
        pending = Path(temporary) / "simulation.mp4"
        plotter = pv.Plotter(off_screen=True, window_size=(args.width, args.height))
        try:
            plotter.set_background("#07111f")
            plotter.enable_anti_aliasing("ssaa")
            cloud = pv.PolyData(positions[0].copy())
            cloud["color"] = colors
            plotter.add_mesh(cloud, scalars="color", rgb=True, style="points",
                             render_points_as_spheres=True, point_size=36 if count <= 16 else 8,
                             show_scalar_bar=False)
            plotter.add_axes(color="#a9b6c9", line_width=2)
            plotter.add_text("N-BODY  /  GRAVITY SIMULATION", position=(0.04, 0.92),
                             viewport=True, font_size=16, color="#dce7f5", font="arial")
            # A bounding sphere fits at every camera azimuth, with no automatic zoom.
            aspect = args.width / args.height
            field = math.radians(30) / 2
            limiting_field = min(field, math.atan(math.tan(field) * aspect))
            radius = float(np.linalg.norm(positions / half, axis=2).max()) * half * 1.15
            radius = radius if radius > 0 else half
            distance = radius / math.sin(limiting_field)
            elevation = math.radians(30)
            plotter.open_movie(str(pending), framerate=args.fps, quality=None,
                               codec="libx264", pixelformat="yuv420p", macro_block_size=1,
                               ffmpeg_params=["-crf", "18", "-movflags", "+faststart"])
            for frame in range(len(positions)):
                cloud.points = positions[frame]
                azimuth = math.radians(45)
                if args.camera == "orbit" and len(positions) > 1:
                    azimuth += 2 * math.pi * frame / (len(positions) - 1)
                eye = distance * np.array([math.cos(elevation) * math.cos(azimuth),
                                           math.cos(elevation) * math.sin(azimuth), math.sin(elevation)])
                plotter.camera_position = (eye, (0, 0, 0), (0, 0, 1))
                plotter.camera.clipping_range = (distance - radius, distance + radius)
                if args.trail_frames:
                    trail = trail_mesh(pv, positions, frame, args.trail_frames, bodies, colors)
                    if trail is not None:
                        plotter.add_mesh(trail, scalars="rgba", rgba=True, preference="cell",
                                         line_width=3, lighting=False, show_scalar_bar=False,
                                         name="trails", reset_camera=False)
                label = "{}  |  {} bodies  |  t = {:.4f}".format(
                    data["metadata"]["solver"], count, data["times"][frame])
                plotter.add_text(label, position="lower_right", font_size=12,
                                 color="#b5c5db", font="arial", name="time")
                plotter.write_frame()
                if frame % 30 == 0 or frame == len(positions) - 1:
                    print("Render {}/{}".format(frame + 1, len(positions)), flush=True)
        finally:
            plotter.close()
        verify_video(pending, len(positions), args.width, args.height)
        publish(pending, output)
    print("Video: {}".format(output), flush=True)
    return output


def main(argv=None):
    options = parser()
    args = options.parse_args(argv)
    try:
        validate_options(args)
        if args.command == "render":
            render(args.input, args.output or args.input.resolve().parent / "simulation.mp4", args)
        else:
            trajectory = capture(args)
            if args.command == "movie":
                render(trajectory, trajectory.parent / "simulation.mp4", args)
    except (ValueError, OSError, ImportError, KeyError, OverflowError, RuntimeError) as error:
        options.exit(1, "Visualization failed: {}\n".format(error))
    except KeyboardInterrupt:
        options.exit(130, "Visualization interrupted; completed trajectory files are preserved.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
