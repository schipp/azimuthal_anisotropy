# pseudo-code (aka all functions are missing) example of how to determine if a ray crosses
# a cell and save the length of the segment within the cell (intended as weight later on)
#
# cells are a list of X,Y- or lon,lat-coordinates.
cells = get_cells()
# rays are a list of ray-end-point (aka station-pair) coordinates. 
rays = get_rays()
# path_lengths an empty list of lists same length as cells
ray_lengths_within_cells = [[] for _ in cells]
for cell_idx, cell in enumerate(cells):
    # determine the four lines that delimit the cell.
    cell_borders = get_cell_borders(cell)
    for ray in rays:
        intersection_points = []
        for cell_border in cell_borders:
            # chech if the line segments intersect
            # get_intersection a function that returns (x, y) if an intersection points is found, and False otherwise
            intersection = get_intersection(line1=cell_border, line2=ray)
            if intersection:
                intersection_points.append(intersection)
        # no intersections
        if len(intersection_points) == 0:
            # both end-points within cell?
            if in_cell(ray[0]) and in_cell(ray[1]):
                point1 = ray[0]
                point2 = ray[1]
            # skip to next ray if end-points oustide of cell and no intersection
            else:
                continue
        # 1 intersection: ray ends or starts in the cell
        elif len(intersection_points) == 1:
            point1 = intersection_points[0]
            # which of the two ray-endpoints is in cell
            if in_cell(ray[0]):
                point2 = ray[0]
            if in_cell(ray[1]):
                point2 = ray[1]
        # 2 intersections: ray crosses the cell
        elif len(intersection_points) == 2:
            # length of path-segment within cell
            point1 = intersection_points[0]
            point2 = intersection_points[1]
        else:
            continue

        len_in_cell = get_distance(point1=point1, point2=point2)
        ray_lengths_within_cells[cell_idx].append(len_in_cell)
# path_density is simply the number of computed segment-lengths within each cell
path_density = [len(_) for _ in ray_lengths_within_cells]