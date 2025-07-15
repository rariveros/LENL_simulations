import numpy as np
import matplotlib.pyplot as plt
import qutip
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':

    def rotation_matrix(azimuth_deg, elevation_deg):
        az = np.deg2rad(azimuth_deg)
        el = np.deg2rad(elevation_deg)
        Rz = np.array([[np.cos(az), -np.sin(az), 0],
                       [np.sin(az), np.cos(az), 0],
                       [0, 0, 1]])
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(el), -np.sin(el)],
                       [0, np.sin(el), np.cos(el)]])
        return Rx @ Rz


    azimuth = 235
    elevation = -70
    R = rotation_matrix(azimuth, elevation)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Sphere
    sphere = plt.Circle((0, 0), 0.99, color='lightgray', alpha=0.7, zorder=1)
    ax.add_artist(sphere)

    trajectories = [
        {'name': 'XZ circle', 'color': 'blue',
         'xyz': lambda t: (np.sin(t), np.zeros_like(t), np.cos(t)),
         't_range': (0, 2 * np.pi)},
        {'name': 'YZ circle', 'color': 'red',
         'xyz': lambda t: (np.zeros_like(t), np.sin(t), np.cos(t)),
         't_range': (0, 2 * np.pi)},
        {'name': 'spiral', 'color': 'green',
         'xyz': lambda t: (
             np.sin(t) * np.cos(4 * np.pi * t / np.pi),
             np.sin(t) * np.sin(4 * np.pi * t / np.pi),
             np.cos(t)),
         't_range': (0, np.pi)}
    ]

    for traj in trajectories:
        t = np.linspace(traj['t_range'][0], traj['t_range'][1], 500)
        x, y, z = traj['xyz'](t)
        points = np.vstack((x, y, z))
        rot_points = R @ points
        z_rot = rot_points[2]

        signs = z_rot >= 0
        transitions = np.where(np.diff(signs))[0]
        segments = []
        start_idx = 0

        for idx in transitions:
            seg = np.arange(start_idx, idx + 2)
            segments.append((seg, signs[start_idx]))
            start_idx = idx + 1
        seg = np.arange(start_idx, len(t))
        segments.append((seg, signs[start_idx]))

        for seg, is_front in segments:
            seg_x = rot_points[0][seg]
            seg_y = rot_points[1][seg]
            alpha = 1.0
            zorder = 3 if is_front else 0
            ax.plot(seg_x, seg_y, color=traj['color'], alpha=alpha, linewidth=2, zorder=zorder,
                    label=traj['name'] if is_front and seg[0] == 0 else None)

    # Draw axes as arrows with front/back separation based only on r<R
    axes = [
        {'vec': np.array([1.4, 0, 0]), 'name': 'x'},
        {'vec': np.array([0, 1.3, 0]), 'name': 'y'},
        {'vec': np.array([0, 0, 1.3]), 'name': 'z'}
    ]

    for ax_def in axes:
        vec = ax_def['vec']
        points = np.linspace([0, 0, 0], vec, 200).T  # shape (3,N)
        rot_points = R @ points
        r = np.linalg.norm(points, axis=0)
        signs = r >= 1.0  # True=adelante, False=atrás

        transitions = np.where(np.diff(signs))[0]
        segments = []
        start_idx = 0

        for idx in transitions:
            seg = np.arange(start_idx, idx + 2)
            segments.append((seg, signs[start_idx]))
            start_idx = idx + 1
        seg = np.arange(start_idx, points.shape[1])
        segments.append((seg, signs[start_idx]))

        for seg, is_front in segments:
            seg_x = rot_points[0][seg]
            seg_y = rot_points[1][seg]
            alpha = 1.0
            zorder = 4 if is_front else 0

            if is_front:
                # Draw the front part as an arrow directly:
                start_x, start_y = seg_x[0], seg_y[0]
                dx, dy = seg_x[-1] - start_x, seg_y[-1] - start_y
                ax.arrow(start_x, start_y, dx, dy,
                         length_includes_head=True,
                         head_width=0.07, head_length=0.12,
                         linewidth=3, color='k', alpha=alpha, zorder=zorder)
            else:
                # Draw the behind part as a normal line
                ax.plot(seg_x, seg_y, color='k', alpha=alpha, linewidth=3, zorder=zorder)

    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_aspect('equal')
    ax.axis('off')

    plt.tight_layout()
    plt.show()