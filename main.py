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
    l1 = LineString(line1)
    l2 = LineString(line2)
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
                # chech if the line segments intersect
                # get_intersection a function that returns (y, x) if an intersection points is found, and False otherwise
                intersection = get_intersection(line1=cell_border, line2=ray)
                if intersection:
                    intersection_points.append(intersection)
            # no intersections
            if len(intersection_points) == 0:
                # both end-points within cell?
                if is_in_cell(point=ray[0], corners=cell_corners) and is_in_cell(point=ray[1], corners=cell_corners):
                    point1 = ray[0]
                    point2 = ray[1]
                # skip to next ray if end-points oustide of cell and no intersection
                else:
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
                # length of path-segment within cell
                point1 = intersection_points[0]
                point2 = intersection_points[1]
            else:
                print('intersection_points should not be longer than 2 entries. Investigate!')
                continue

            len_in_cell = get_distance(point1=point1, point2=point2, coord_type=coord_type)
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

    fig, ax = plt.subplots(1, 1)

    gps = np.array(grid_points)
    xx, yy = np.meshgrid(grid_x, grid_y)

    pd = np.array(path_density).reshape(len(grid_x), len(grid_y))
    img = ax.pcolormesh(xx-cell_width/2, yy-cell_height/2, pd, norm=norm, cmap='viridis_r')
    
    scs = np.array(station_coords)
    ax.scatter(scs[:,1], scs[:,0], marker='^', facecolors='k', edgecolors='w')
    
    x0, y0, w, h = ax.get_position().bounds
    cb_ax = fig.add_axes([x0+w+.01*w, y0, .02*w, h])
    plt.colorbar(img, cax=cb_ax)
    
    plt.show(fig)
    plt.close(fig)