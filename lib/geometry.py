def generate_grid(lats:tuple, lons:tuple, cellsize:float, grid_overlap:float) -> list:
    """
    Generates a grid of `cellsize` large cells in bounds `lats`/`lons` with `grid_overlap` overlap.

    :param lats: (lat_min, lat_max) lower and upper latitude limit for grid points
    :param lons: (lon_min, lon_max) lower and upper longitude limit for grid points
    :param cellsize: size of cells in decimal degrees
    :param grid_overlap: overlap-ratio of cells [0 to 1]
    """
    import numpy as np
    
    dist_lon = np.max(lons) - np.min(lons)
    dist_lat = np.max(lats) - np.min(lats)

    grid_point_non_overlap = 1 - grid_overlap
    gp_spacing_lon = (cellsize / dist_lon) * (np.max(lons) - np.min(lons)) * grid_point_non_overlap
    gp_spacing_lat = (cellsize / dist_lat) * (np.max(lats) - np.min(lats)) * grid_point_non_overlap

    n_gp_lon = int(np.ceil(dist_lon / cellsize) / grid_point_non_overlap) + 1
    n_gp_lat = int(np.ceil(dist_lat / cellsize) / grid_point_non_overlap) + 1

    grid_points_lon_coords = [np.min(lons) + ix * gp_spacing_lon for ix in range(n_gp_lon)]
    grid_points_lat_coords = [np.max(lats) - iy * gp_spacing_lat for iy in range(n_gp_lat)]

    grid_points = [
        [np.max(lats) - iy * gp_spacing_lat, np.min(lons) + ix * gp_spacing_lon] for iy in range(n_gp_lat) for ix in range(n_gp_lon)
        ]

    gp_meshgrid = np.meshgrid(grid_points_lon_coords, grid_points_lat_coords)

    return grid_points, gp_meshgrid


def get_cell_corners(cell:list, width:float, height:float, rel_loc='mid') -> list:
    """
    Returns the corners of a cell (x1, x2, y1, y2).
    
    :param cell: A tuple of lat/Y,lon/X-coordinates of the cell
    :param width: Width of the cell
    :param height: Height of the cell
    :param rel_loc: Which position of the cell does `cell` indicate: center (mid)
    :returns: y1, y2, x1, x2

    :Example:

    """
    y, x = cell

    if rel_loc == 'mid':
        x1 = x - width/2
        x2 = x + width/2
        y1 = y - height/2
        y2 = y + height/2
    elif rel_loc == 'ul':
        x1 = x
        x2 = x + width
        y1 = y - height
        y2 = y
    else:
        raise AttributeError('Invalid rel_loc. Supported: mid, ul')
    
    corners = [
        [y2, x1],
        [y2, x2],
        [y1, x2],
        [y1, x1],
        ]

    return corners


def get_cell_borders(corners:list) -> list:
    """
    Computes the start- and end-points of the four lines delimiting a cell.
    
    :param cell_corners: List of 4 tuples (corner-coordinates in Y,X)
    :returns: List of 4 lists (= 4 cell-borders) containg 2 tuples (= start- and end-points in [Y, X])

    :Example:

    """

    cell_borders = [
        [corners[0], corners[1]],
        [corners[1], corners[2]],
        [corners[2], corners[3]],
        [corners[3], corners[0]],
        ]

    return cell_borders

def is_in_cell(point:list, corners:list) -> bool:
    """
    Checks if a point is within a cell.
    
    :param point: Tuple of lat/Y,lon/X-coordinates
    :param corners: List of corner coordinates
    :returns: Boolean whether point is within cell

    :Example:

    """    
    y1, y2, x1, x2 = corners[2][0], corners[0][0], corners[0][1], corners[2][1]

    if (y1 <= point[0] <= y2) and (x1 <= point[1] <= x2):
        return True
    return False

def get_distance(point1:tuple, point2:tuple, coord_type='geographic') -> float:
    """
    Compute distance between two points
    
    :param point1: Tuple of lat/Y,lon/X-coordinates
    :param point2: Tuple of lat/Y,lon/X-coordinates
    :param coord_type: Type of coordinate system. 'cartesian' or 'geographic'
    :returns: Distance between the points in meters

    :Example:

    """    

    if coord_type == 'geographic':
        from obspy.geodetics import gps2dist_azimuth
        return gps2dist_azimuth(lat1=point1[0], lon1=point1[1], lat2=point2[0], lon2=point2[1])
    elif coord_type == 'cartesian':
        import numpy as np
        return np.sqrt((np.array(point1) - np.array(point2))**2)
    else:
        raise AttributeError('coord_type invalid: Supported: cartesian or geographic')


def get_intersection(line1:list, line2:list) -> list:
    """
    Check if line1 and line2 intersect and returns intersection-point coordinates as list.
    Returns empty list, if no intersection point.
    
    :param line1: List of 2 tuples of lat/Y,lon/X-coordinates
    :param line1: List of 2 tuples of lat/Y,lon/X-coordinates
    :returns: Tuple of intersection-point coordinates [Y, X]. Empty list, if no intersection.

    :Example:

    """  

    from shapely.geometry import LineString
    l1 = LineString(((line1[0][1], line1[0][0]), (line1[1][1], line1[1][0])))
    l2 = LineString(((line2[0][1], line2[0][0]), (line2[1][1], line2[1][0])))
    if l1.intersects(l2):
        point = l1.intersection(l2)
        return [point.y, point.x]
    return []



def get_ray_lengths_within_cells(rays:list, cells:list, cell_width:float, cell_height:float, rel_loc='mid', coord_type='geographic') -> list:
    """
    Main logic of checking all cells for crossing rays and saving the lengths of each ray within cells.
    
    :param rays: List of tuples of tuples of station-pair lat/Y,lon/X-coordinates
    :param cells: List of tuples of lat/Y,lon/X-coordinates of the cell
    :param cell_width: Width of the cell
    :param cell_height: Height of the cell
    :param rel_loc: Which position of the cell does `cell` indicate: center (mid)
    :param coord_type: Type of coordinate system. 'cartesian' or 'geographic'
    :returns: List with length len(cells), each entry containing the ids and lenghts of all crossing rays for each cell.

    :Example:

    """
    # pseudo-code (aka all functions are missing) example of how to determine if a ray crosses
    # a cell and save the length of the segment within the cell (intended as weight later on)
    #
    # cells are a list of X,Y- or lon,lat-coordinates.
    # rays are a list of ray-end-point (aka station-pair) coordinates. 

    # path_lengths an empty list of lists same length as cells
    ray_lengths_within_cells = [[] for _ in cells]
    for cell, ray_lengths_within_cell in zip(cells, ray_lengths_within_cells):
        cell_corners = get_cell_corners(cell=cell, width=cell_width, height=cell_height, rel_loc=rel_loc)
        # determine the four lines that delimit the cell.
        cell_borders = get_cell_borders(corners=cell_corners)
        for ray_idx, ray in enumerate(rays):
            intersection_points = []
            for cell_border in cell_borders:
                # check if and where the line segments intersect
                intersection = get_intersection(line1=cell_border, line2=ray)
                if intersection:
                    intersection_points.append(intersection)
            # no intersections
            if len(intersection_points) == 0:
                # both end-points within cell?
                if is_in_cell(point=ray[0], corners=cell_corners) and is_in_cell(point=ray[1], corners=cell_corners):
                    point1 = ray[0]
                    point2 = ray[1]
                else:
                    # skip to next ray if end-points oustide of cell and no intersection
                    continue
            # 1 intersection: ray ends or starts in the cell
            elif len(intersection_points) == 1:
                point1 = intersection_points[0]
                # which of the two ray-endpoints is in cell
                if is_in_cell(point=ray[0], corners=cell_corners):
                    point2 = ray[0]
                if is_in_cell(point=ray[1], corners=cell_corners):
                    point2 = ray[1]
            # 2 intersections: ray crosses the cell
            elif len(intersection_points) == 2:
                point1 = intersection_points[0]
                point2 = intersection_points[1]
            else:
                print('intersection_points should not be longer than 2 entries. Investigate!')
                continue

            len_in_cell = get_distance(point1=point1, point2=point2, coord_type=coord_type)

            if len_in_cell[0] > 4000E3:
                print(len(intersection_points), intersection_points)
                print('')
                print(ray[0])
                print(ray[1])
                print('')
                print(point1)
                print(point2)
                print('')
            ray_lengths_within_cell.append([ray_idx, len_in_cell])

    return ray_lengths_within_cells

if __name__ == '__main__':
    # Run a random example (random distribution of stations within a grid)
    import numpy as np
    
    # I use lat/lon coordinates for everything (i.e, Y/X)
    # construct grid
    from itertools import product, combinations
    cell_width = .25
    cell_height = .25
    grid_x = np.arange(10, 20, cell_width)
    grid_y = np.arange(40, 50, cell_height)
    grid_points = list(product(grid_y, grid_x))

    # construct rays
    station_x = np.random.uniform(11, 19, [10])
    station_y = np.random.uniform(41, 49, [10])
    station_coords = list(zip(station_y, station_x))
    station_pairs = list(combinations(station_coords, 2))

    ray_lengths_within_cells = get_ray_lengths_within_cells(
        rays=station_pairs,
        cells=grid_points,
        cell_width=cell_width,
        cell_height=cell_height,
        rel_loc='mid',
        coord_type='geographic')
    
    # path_density is simply the total amount of computed segment-lengths within each cell
    path_density = [len(_) for _ in ray_lengths_within_cells]

    # plot example
    import pylab as plt
    from matplotlib.colors import Normalize

    norm = Normalize(vmin=0, vmax=max(path_density))

    fig, axs = plt.subplots(1, 2)
    ax = axs[0]
    ax.set_title('path density')
    ax.set_aspect('equal')

    gps = np.array(grid_points)
    xx, yy = np.meshgrid(grid_x, grid_y)

    pd = np.array(path_density).reshape(len(grid_x), len(grid_y))
    img = ax.pcolormesh(xx-cell_width/2, yy-cell_height/2, pd, norm=norm, cmap='viridis_r')
    
    scs = np.array(station_coords)
    ax.scatter(scs[:,1], scs[:,0], marker='^', facecolors='k', edgecolors='w')
    
    x0, y0, w, h = ax.get_position().bounds
    # cb_ax = fig.add_axes([x0+w+.01*w, y0, .02*w, h])
    cb_ax = fig.add_axes([x0+.025*w, y0+.025*h, .02*w, .95*h])
    plt.colorbar(img, cax=cb_ax)
    

    # visualize the densest cell to check if everythings working
    idx_dense = np.argmax(path_density)
    gp_dense = grid_points[idx_dense]
    rayidx_dense = np.array(ray_lengths_within_cells[idx_dense])[:,0]
    ray_params = np.array(ray_lengths_within_cells[idx_dense])[:,1]
    
    corners_of_dense = get_cell_corners(cell=gp_dense, width=cell_width, height=cell_height)
    borders_of_dense = get_cell_borders(corners=corners_of_dense)
    for border in borders_of_dense:
        ax.plot((border[0][1], border[1][1]), (border[0][0], border[1][0]), c='red', lw=2)

    ax = axs[1]
    ax.set_title('zoom-in: densest cell')
    ax.set_aspect('equal')

    for ray_idx, ray_param in zip(rayidx_dense, ray_params):
        sta1_x = station_pairs[ray_idx][0][1]
        sta1_y = station_pairs[ray_idx][0][0]
        sta2_x = station_pairs[ray_idx][1][1]
        sta2_y = station_pairs[ray_idx][1][0]
        ax.plot(
            [sta1_x, sta2_x], [sta1_y, sta2_y],
            label=f"{ray_param[0]/1000:0.1f}km / {np.min(ray_param[1:]):0.1f}°"
            )
    ax.scatter(scs[:,1], scs[:,0], s=80, marker='^', facecolors='k', edgecolors='w', zorder=10)
    ax.set_xlim(gp_dense[1]-cell_width/2, gp_dense[1]+cell_width/2)
    ax.set_ylim(gp_dense[0]-cell_height/2, gp_dense[0]+cell_height/2)
    ax.set_xticks((gp_dense[1]-cell_width/2, gp_dense[1]+cell_width/2))
    ax.set_yticks((gp_dense[0]-cell_height/2, gp_dense[0]+cell_height/2))
    ax.spines['bottom'].set_color('red')
    ax.spines['top'].set_color('red') 
    ax.spines['right'].set_color('red')
    ax.spines['left'].set_color('red')

    ax.legend(loc=2, bbox_to_anchor=(1, 1), fontsize=6)

    plt.show(fig)
    plt.close(fig)