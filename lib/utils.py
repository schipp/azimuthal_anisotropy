import numpy as np

def bin_az_measurements(vels, bazs, az_step, weights):

    az_bins = np.arange(0, 180+az_step, az_step)
    vel_bin_medians = []
    vel_bin_stds = []
    for az_bin in az_bins:
        in_bin_idx = np.where((az_bin <= bazs) & (bazs < az_bin + az_step))
        vels_in_bin = vels[in_bin_idx]
        weights_in_bin = np.array(weights[in_bin_idx])
        if len(vels_in_bin) > 1:
            # weighted by distances
            # vel_bin_medians.append(weighted_median(vels, dists_in_bin))
            vel_bin_medians.append(np.average(vels_in_bin, weights=weights_in_bin/np.max(weights_in_bin)))
            vel_bin_stds.append(np.nanstd(vels_in_bin))
        else:
            vel_bin_medians.append(np.nan)
            vel_bin_stds.append(np.nan)
    
    az_bins = np.delete(az_bins, np.argwhere(np.isnan(vel_bin_medians)))
    vel_bin_medians = np.delete(vel_bin_medians, np.argwhere(np.isnan(vel_bin_medians)))
    vel_bin_stds = np.delete(vel_bin_stds, np.argwhere(np.isnan(vel_bin_stds)))

    return az_bins, vel_bin_medians, vel_bin_stds

def chunks(lst, n):
    """Yield n chunks and their indices from lst."""
    import numpy as np

    chunk_size = len(lst) // n
    if len(lst) % 2 != 0:
        chunk_size += 1
    for i in range(0, len(lst), chunk_size):
        chunk = lst[i:i + chunk_size]
        chunk_indices = np.arange(i, i + chunk_size)
        yield chunk_indices, chunk