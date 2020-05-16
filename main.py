import numpy as np
from lib.geometry import generate_grid, get_ray_lengths_within_cells, get_distance
from lib.anisotropy import aniso_parametrization, bin_az_measurements, fit_aniso, extract_aniso_results
import logging

def generate_synthetic_example(n_stations=65, lat_lims=(45, 50), lon_lims=(10, 20), u_0=0, A=.1, phi_2=60, B=.01, phi_4=20, amplitude_noise=.05):
    """
    Helper function to generate a simple synthetic example.
    Constant anisotropy in entire region.

    :param n_stations: Number of stations in study area
    :param lat_lims: Latitude limits of study area
    :param lon_lims: Longitude limits of study area
    :param A: Amplitude of 2Theta-term
    :param phi_2: Phase of 2Theta-term
    :param B: Amplitude of 4Theta-term
    :param phi_4: Phase of 4Theta-term
    :param u_0: Isotropic velocity
    :returns: station pair locations  and velocities
    """
    from itertools import combinations

    # random station distribution inside grid
    station_x = np.random.uniform(lon_lims[0]-1, lon_lims[1]+1, [n_stations])
    station_y = np.random.uniform(lat_lims[0]-1, lat_lims[1]+1, [n_stations])
    station_coords = list(zip(station_y, station_x))
    station_pairs = list(combinations(station_coords, 2))

    # compute synthetic velocities for station pairs
    vels = []
    for sp in station_pairs:
        dist, az, baz = get_distance(
            point1=sp[0],
            point2=sp[1],
            coord_type='geographic')

        vel_for_sp = aniso_parametrization(
            x=np.deg2rad(np.min([az, baz])),
            A=A,
            B=B,
            u_0=u_0,
            phi_2=np.deg2rad(phi_2),
            phi_4=np.deg2rad(phi_4)
            )
        
        # add noise
        vel_for_sp += np.random.uniform(-amplitude_noise, amplitude_noise)
        
        vels.append(vel_for_sp)

    return station_pairs, vels

if __name__ == '__main__':
    # -- Synthetic or Real data?
    do_synthetic_example = False
    csv_file = 'synth.csv'
    # -- Geometry of grid
    # size of cells in decimal degrees
    grid_point_cellsize_degree = 1.5
    # overlap of cells from 0 to 1.
    grid_point_overlap = .5
    # -- Curve-fitting
    # size of azimuth-bins in degrees
    az_bin_size = 5

    # define grid
    lat_lims = [45, 50]
    lon_lims = [10, 20]
    grid_points, gp_meshgrid = generate_grid(
        lats=lat_lims,
        lons=lon_lims,
        cellsize=grid_point_cellsize_degree,
        grid_overlap=grid_point_overlap
        )

    # get station pair locations and velocities
    if do_synthetic_example:
        station_pairs, vels = generate_synthetic_example()
        output_format = np.empty([len(station_pairs), 5])
        station_pairs = np.array(station_pairs)
        # output_format[:, 0:4] = station_pairs.reshape(station_pairs.shape[0], 4)
        # output_format[:, 4] = np.array(vels)
        # np.savetxt(csv_file, output_format, delimiter=',')

    else:
        real_data = np.loadtxt(csv_file, delimiter=',')
        sps_1 = real_data[:, 0:2]
        sps_2 = real_data[:, 2:4]
        station_pairs = np.array([[sp1, sp2] for sp1, sp2 in zip(sps_1, sps_2)])
        vels = real_data[:, 4]
    
    # get lengths of rays within each cell
    ray_lengths_within_cells = get_ray_lengths_within_cells(
        rays=station_pairs,
        cells=grid_points,
        cell_width=grid_point_cellsize_degree,
        cell_height=grid_point_cellsize_degree,
        rel_loc='mid',
        coord_type='geographic')

    # analyze anisotropy for each cell
    fast_directions = np.zeros(len(grid_points))
    fast_directions[:] = np.nan
    fast_directions_errors = np.zeros(len(grid_points))
    fast_directions_errors[:] = np.nan
    fast_amps_A = np.zeros(len(grid_points))
    fast_amps_A[:] = np.nan
    fast_amps_A_errors = np.zeros(len(grid_points))
    fast_amps_A_errors[:] = np.nan

    for idx, (cell, ray_params) in enumerate(zip(grid_points, ray_lengths_within_cells)):
        # cell has no crossing rays
        if len(ray_params) == 0:
            continue
        # extract ray parameters for cell
        ray_params = np.array(ray_params)
        ray_idxs = ray_params[:, 0]
        ray_geometry_in_cell = ray_params[:, 1]
        ray_lens = [_[0] for _ in ray_geometry_in_cell]
        ray_azs = [_[1] for _ in ray_geometry_in_cell]
        ray_bazs = [_[2] for _ in ray_geometry_in_cell]

        current_vels = np.array(vels)[list(ray_idxs)]
        current_azis = np.min(np.array([ray_azs, ray_bazs]), axis=0)

        # bin velocity measurements into azimuth-bins
        az_bins, vel_bin_medians, vel_bin_stds = bin_az_measurements(
            vels=current_vels,
            bazs=current_azis,
            az_step=az_bin_size
            )

        # fit anisotropic parameters
        if len(vel_bin_medians) < 8:
            logging.info(f"GP {cell[0]:0.2f}/{cell[1]:0.2f}: not enough azimuth-bins w/ data ({len(vel_bin_medians)}) to resolve anisotropy.")
            continue

        params, params_covariance = fit_aniso(
            az_bins=az_bins,
            az_step=az_bin_size,
            vel_bin_medians=vel_bin_medians,
            vel_bin_stds=vel_bin_stds
            )

        phi_2, phi_2_std, rel_A, rel_A_std = extract_aniso_results(params, params_covariance)

        # sanity checks, these will trigger if unstable or overfitted
        if np.isinf(phi_2_std):
            logging.info(f"GP {cell[0]:0.2f}/{cell[1]:0.2f}: sanity check #1 failed. phi_2_std == ∞")
            continue
        if phi_2_std < 0.01:
            logging.info(f"GP {cell[0]:0.2f}/{cell[1]:0.2f}: sanity check #2 failed. phi_2_std < 0.01")
            continue

        fast_directions[idx] = phi_2
        fast_directions_errors[idx] = phi_2_std
        fast_amps_A[idx] = np.abs(rel_A)
        fast_amps_A_errors[idx] = rel_A_std

        ## cell fit visualization
        # import pylab as plt
        # fig, ax = plt.subplots(1, 1)
        # test_az = np.arange(0, 180+1, 1)
        # test_vels = aniso_parametrization(x=test_az, A=rel_A, B=rel_B, u_0=0, phi_2=phi_2, phi_4=phi_4)
        # ax.scatter(current_azis, current_vels, c='#CCCCCC')
        # ax.errorbar(az_bins+az_bin_size/2, vel_bin_medians, fmt='o', yerr=vel_bin_stds, capsize=4, c='k')
        # ax.plot(test_az, test_vels)
        # ax.set_xlim(0, 180)
        # plt.show(fig)
        # plt.close(fig)
        # break

    # visualize results
    import pylab as plt
    fig, axs = plt.subplots(1, 2)
    
    # transform radial to cartesian coordinates
    x_comp = np.multiply(fast_amps_A, np.sin(np.radians(fast_directions)))
    y_comp = np.multiply(fast_amps_A, np.cos(np.radians(fast_directions)))

    xx, yy = gp_meshgrid

    ax = axs[0]
    ax.set_aspect('equal')
    ax.scatter(
        [station_pairs[:, 0, 1], station_pairs[:, 1, 1]], 
        [station_pairs[:, 0, 0], station_pairs[:, 1, 0]],
        marker='^')

    Q = ax.quiver(
        xx,
        yy,
        x_comp.T,
        y_comp.T,
        color='#d62728',
        width=0.008,
        headlength=0,
        headwidth=0,
        headaxislength=0,
        scale=50000,
        pivot='mid'
        )
    
    ax.set_xlim(lon_lims)
    ax.set_ylim(lat_lims)
    
    ax = axs[1]
    ax.set_aspect('equal')
    ax.scatter(
        [station_pairs[:, 0, 1], station_pairs[:, 1, 1]], 
        [station_pairs[:, 0, 0], station_pairs[:, 1, 0]],
        marker='^')

    x_comp_fixed_len = np.multiply(1, np.sin(np.radians(fast_directions)))
    y_comp_fixed_len = np.multiply(1, np.cos(np.radians(fast_directions)))

    Q = ax.quiver(
        xx,
        yy,
        x_comp_fixed_len.T,
        y_comp_fixed_len.T,
        color='#d62728',
        width=0.008,
        headlength=0,
        headwidth=0,
        headaxislength=0,
        scale=10,
        pivot='mid'
        )
    
    ax.set_xlim(lon_lims)
    ax.set_ylim(lat_lims)

    plt.show(fig)
    plt.close(fig)