"""Small plotting helpers for :class:`autopas.AutoPas`. Requires matplotlib."""

from __future__ import annotations

import itertools

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Categorical colors in a fixed order, so that a particle type keeps its color from plot to plot.
_CATEGORICAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")
_SEQUENTIAL = "viridis"
_INK = "#52514e"
_GRID = "#e4e3df"


def _style_axes(ax) -> None:
    """Keeps grid and axes recessive."""
    ax.grid(True, color=_GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(_GRID)
    ax.tick_params(colors=_INK, length=0)


def _draw_box(ax, box_min, box_max) -> None:
    """Draws the edges of the domain."""
    corners = list(itertools.product(*zip(box_min, box_max)))
    for a, b in itertools.combinations(corners, 2):
        if sum(x != y for x, y in zip(a, b)) == 1:
            ax.plot(*zip(a, b), color=_INK, linewidth=0.6, alpha=0.6)


def _largest_speed(speed: np.ndarray) -> float:
    """The upper end of a color scale for these speeds, which stays valid if everything is at rest."""
    return max(float(speed.max(initial=0.0)), 1e-12)


def _scatter(ax, sim, color_by: str, size: float, speed_max: float | None = None):
    """Adds the particles to a 3D axes, and a legend or color bar for what the color means."""
    data = sim.particles()
    x, y, z = data["positions"].T
    if color_by == "type":
        artists = []
        types = data["type_ids"]
        for type_id in np.unique(types):
            mask = types == type_id
            color = _CATEGORICAL[int(type_id) % len(_CATEGORICAL)]
            artists.append(ax.scatter(x[mask], y[mask], z[mask], s=size, color=color, label=f"type {type_id}"))
        if len(artists) > 1:
            ax.legend(frameon=False, labelcolor=_INK)
        return artists
    if color_by == "speed":
        speed = np.linalg.norm(data["velocities"], axis=1)
        vmax = speed_max if speed_max is not None else _largest_speed(speed)
        artist = ax.scatter(x, y, z, s=size, c=speed, cmap=_SEQUENTIAL, vmin=0.0, vmax=vmax)
        ax.figure.colorbar(artist, ax=ax, shrink=0.6, label="speed")
        return [artist]
    raise ValueError(f"color_by has to be 'type' or 'speed', not {color_by!r}")


def _new_3d_axes(sim, ax, title: str | None):
    """Sets up the axes with the true proportions of the domain."""
    if ax is None:
        ax = plt.figure(figsize=(6.4, 5.4)).add_subplot(projection="3d")
    box_min, box_max = np.asarray(sim.box_min), np.asarray(sim.box_max)
    _draw_box(ax, box_min, box_max)
    ax.set(xlim=(box_min[0], box_max[0]), ylim=(box_min[1], box_max[1]), zlim=(box_min[2], box_max[2]))
    ax.set_box_aspect(box_max - box_min)
    ax.set(xlabel="x", ylabel="y", zlabel="z")
    if title is not None:
        ax.set_title(title, color=_INK, loc="left")
    return ax


def plot_particles(sim, *, ax=None, color_by: str = "type", size: float = 8.0, title: str | None = None):
    """Plots the current particles and the domain.

    Args:
        sim: The simulation.
        ax: A 3D axes to draw on. By default, creates a new figure.
        color_by: Either "type", for the particle type, or "speed".
        size: Size of the markers.
        title: Title of the plot. Defaults to the iteration.

    Returns:
        The axes.
    """
    ax = _new_3d_axes(sim, ax, f"Iteration {sim.iteration}" if title is None else title)
    _scatter(ax, sim, color_by, size)
    return ax


def plot_history(sim, *, axes=None):
    """Plots what the simulation recorded per iteration: the energies and the temperature.

    The two quantities have different units, so each one gets its own panel. Iterations in which
    the auto tuner was testing a configuration are shaded, as this is where a run is at its slowest.

    Args:
        sim: The simulation.
        axes: Two axes to draw on. By default, creates a new figure.

    Returns:
        The two axes.
    """
    history = sim.history
    if len(history["iteration"]) == 0:
        raise ValueError("There is nothing to plot yet. Run the simulation with a history_interval > 0.")
    if axes is None:
        _, axes = plt.subplots(2, 1, sharex=True, figsize=(6.4, 4.8), constrained_layout=True)
    energy_ax, temperature_ax = axes
    iteration = history["iteration"]
    potential, kinetic = history["potential_energy"], history["kinetic_energy"]

    energy_ax.plot(iteration, kinetic, color=_CATEGORICAL[0], linewidth=2, label="kinetic")
    if np.isnan(potential).all():
        energy_ax.set_title("Energy (the functors provide no potential energy)", color=_INK, loc="left")
    else:
        energy_ax.plot(iteration, potential, color=_CATEGORICAL[1], linewidth=2, label="potential")
        energy_ax.plot(iteration, potential + kinetic, color=_CATEGORICAL[2], linewidth=2, label="total")
        energy_ax.set_title("Energy", color=_INK, loc="left")
    energy_ax.legend(frameon=False, labelcolor=_INK)

    temperature_ax.plot(iteration, history["temperature"], color=_CATEGORICAL[0], linewidth=2)
    temperature_ax.set_title("Temperature", color=_INK, loc="left")
    temperature_ax.set_xlabel("iteration", color=_INK)

    for ax in axes:
        _style_axes(ax)
        ax.fill_between(iteration, 0, 1, where=history["tuning"], transform=ax.get_xaxis_transform(),
                        color=_GRID, alpha=0.6, linewidth=0, step="mid")
    return axes


def animate_particles(sim, frames: int, *, iterations_per_frame: int = 10, color_by: str = "type",
                      size: float = 8.0, interval: int = 50, speed_max: float | None = None) -> FuncAnimation:
    """Advances the simulation while showing it. The simulation stays usable afterwards.

    Args:
        sim: The simulation.
        frames: The number of frames.
        iterations_per_frame: How many iterations to run per frame.
        color_by: Either "type", for the particle type, or "speed".
        size: Size of the markers.
        interval: Time between two frames in milliseconds.
        speed_max: The upper end of the color scale for color_by="speed". By default, follows the fastest particle of
            each frame, so the same color can mean different speeds in different frames.

    Returns:
        The animation. Keep a reference to it, or matplotlib stops it. Use its save() method to get a file, e.g.
        ``animation.save("drop.gif", writer="pillow", fps=12)``.
    """
    ax = _new_3d_axes(sim, None, None)
    artists = _scatter(ax, sim, color_by, size, speed_max)

    def update(_):
        nonlocal artists
        sim(iterations_per_frame)
        data = sim.particles()
        ax.set_title(f"Iteration {sim.iteration}", color=_INK, loc="left")
        if color_by == "speed":
            artists[0]._offsets3d = tuple(data["positions"].T)
            speed = np.linalg.norm(data["velocities"], axis=1)
            artists[0].set_array(speed)
            if speed_max is None:
                artists[0].set_clim(0.0, _largest_speed(speed))
        else:
            for artist, type_id in zip(artists, np.unique(data["type_ids"])):
                artist._offsets3d = tuple(data["positions"][data["type_ids"] == type_id].T)
        return artists

    # Without an init_func, matplotlib draws the first frame an additional time, which would advance the simulation.
    return FuncAnimation(ax.figure, update, init_func=lambda: artists, frames=frames, interval=interval, blit=False,
                         repeat=False)
