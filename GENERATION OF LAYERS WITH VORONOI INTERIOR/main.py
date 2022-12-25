# -*- coding: utf-8 -*-
"""
Ege Uğur Aguş 
"""

from collections import defaultdict
import numpy as np
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
import random

def point_from_poly(poly, n, centered = True):
    """
    Returns a list of points that contains 'n' number uniformly distributed 
    points from given polygon 'poly'
    """
    min_x, min_y, max_x, max_y = poly.bounds # Find the min and max values of the boundary polygon
    point_list = [] # Create a list to contain points

    if centered == True:
        
        point_list.append([poly.centroid.x, poly.centroid.y])
        angle = 2*np.pi/n
        
        for i in range(n-1):
            point_list.append([poly.centroid.x + np.cos(i*angle), poly.centroid.y + np.sin(i*angle)])
    else:

        while len(point_list) < n: 
            p = Point([random.uniform(min_x, max_x), random.uniform(min_y, max_y)]) # Pick a random point from uniform dist.
            if poly.contains(p): # Check if the polygon contains the point
                point_list.append([p.x,p.y])
                
    return np.array(point_list) # Return an numpy array


def voronoi_polygons(point_list, diameter):
    """Generate shapely.geometry.Polygon objects corresponding to the
    regions of a scipy.spatial.Voronoi object, in the order of the
    input points. The polygons for the infinite regions are large
    enough that all points within a distance 'diameter' of a Voronoi
    vertex are contained in one of the infinite polygons.

    """
    voronoi = Voronoi(point_list) # Create Voronoi obj
    centroid = voronoi.points.mean(axis=0) # Find 

    # Mapping from (input point index, Voronoi point index) to list of
    # unit vectors in the directions of the infinite ridges starting
    # at the Voronoi point and neighbouring the input point.
    ridge_direction = defaultdict(list)
    for (p, q), rv in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        u, v = sorted(rv)
        if u == -1:
            # Infinite ridge starting at ridge point with index v,
            # equidistant from input points with indexes p and q.
            t = voronoi.points[q] - voronoi.points[p] # tangent
            n = np.array([-t[1], t[0]]) / np.linalg.norm(t) # normal
            midpoint = voronoi.points[[p, q]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centroid, n)) * n
            ridge_direction[p, v].append(direction)
            ridge_direction[q, v].append(direction)

    for i, r in enumerate(voronoi.point_region):
        region = voronoi.regions[r]
        if -1 not in region:
            # Finite region.
            yield Polygon(voronoi.vertices[region])
            continue
        # Infinite region.
        inf = region.index(-1)              # Index of vertex at infinity.
        j = region[(inf - 1) % len(region)] # Index of previous vertex.
        k = region[(inf + 1) % len(region)] # Index of next vertex.
        if j == k:
            # Region has one Voronoi vertex with two ridges.
            dir_j, dir_k = ridge_direction[i, j]
        else:
            # Region has two Voronoi vertices, each with one ridge.
            try:
                dir_j, = ridge_direction[i, j]
                dir_k, = ridge_direction[i, k]
            except:
                continue

        # Length of ridges needed for the extra edge to lie at least
        # 'diameter' away from all Voronoi vertices.
        length = 2 * diameter / np.linalg.norm(dir_j + dir_k)

        # Polygon consists of finite part plus an extra edge.
        finite_part = voronoi.vertices[region[inf + 1:] + region[:inf]]
        extra_edge = [voronoi.vertices[j] + dir_j * length,
                      voronoi.vertices[k] + dir_k * length]
        yield Polygon(np.concatenate((finite_part, extra_edge)))

        
def voronoi_iteration(ncells, diameter, boundary_polygon, centered = True, max_iter_num = 1000, min_th = 0.8, max_th = 1.2):
    """
    Lloyd's algorithm also known as Voronoi iteration
    Allows us to find approximate centers to make voronoi cells uniform
    Find the centroids of the voronoi cells and assign those centroids as new
    voronoi points until the smallest voronoi cell has at least min th of the ideal area and
    the biggest voronoi cell has at most max th of the ideal area   
    """       
    
    points = point_from_poly(boundary_polygon,ncells,centered) # Random points are selected uniformly
    ideal_area = boundary_polygon.area / ncells # Calculate ideal area 
    
    for i in range(max_iter_num): # Start iterations
        centroids = []
        areas = []
        cells = []
        for p in voronoi_polygons(points, diameter):
            x, y = zip(*p.intersection(boundary_polygon).exterior.coords)
            iteration_poly = Polygon(list(zip(x,y)))
            centroids.append([iteration_poly.centroid.x,iteration_poly.centroid.y])
            cells.append(iteration_poly)
            areas.append(iteration_poly.area)       
        points = np.array(centroids)
        if min(areas)/ideal_area > min_th and max(areas) / ideal_area < max_th:
            break
    
    return points, cells   
        
def Main(boundary, ncells, height, thickness, constant_bl = True, centered = True, next_boundaries = []):
    """
    Main function that performs the operations
    """
    
    ax = plt.axes(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    
    nlayers = int((height/thickness)) + 1
    
    # Calculate first layer
    boundary_polygon = Polygon(boundary) # Construct boundary polygon
    diameter = np.linalg.norm(boundary.ptp(axis=0)) # Find maximum radii for converting infinite cells to polygons
    vorpoints, cells = voronoi_iteration(ncells,diameter,boundary_polygon,centered) # Obtain sites and cells

    boundx, boundy = boundary.T
    plt.xlim(round(boundx.min() - 1), round(boundx.max() + 1))
    plt.ylim(round(boundy.min() - 1), round(boundy.max() + 1))

    for cell in cells:
        x, y = zip(*cell.exterior.coords)
        # plt.plot(x, y, 'r-') # 2D plot
        ax.plot3D(x, y, 0,c='red') # 3D plot
    
    # plt.plot(*vorpoints.T, 'b.') # 2D plot
    ax.scatter(vorpoints[:,0], vorpoints[:,1], 0,c='black') # 3D plot

    past_vorpoints = vorpoints

    k = 0.2 # Amplifies the movement of the Voronoi cells at 1 layer

    cell_dir = defaultdict(list)
    
    for c, point in enumerate(vorpoints):
        if c == 0:
            # sign = np.random.choice([-1,1]) # Random initial behavior of the cells
            sign = 1 # Set initial behavior as all the points gather into one point

            aim = list(cells[0].exterior.coords)[0] # This aim can be randomized accoriding to desired properties
            direction = [(aim[0]-vorpoints[0][0])*k*sign, (aim[1]-vorpoints[0][1])*k*sign]
            
            cell_dir[0] = direction
            continue
        
        direction = [(aim[0]-vorpoints[c][0])*k*sign, (aim[1]-vorpoints[c][1])*k*sign]

        cell_dir[c] = direction

    for i in range(1,nlayers):
        
        if constant_bl == False:
            boundary_polygon = Polygon(next_boundaries[i-1]) # Construct boundary polygon
        
        
        ideal_area = boundary_polygon.area / ncells # Calculate ideal area    
        new_vorpoints = [[sum(x) for x in zip(cell_dir[c], past_vorpoints[c])] for c in range(len(past_vorpoints))]

        areas = [cell.area for cell in cells]
        
        for c,m in enumerate(new_vorpoints):
            if boundary_polygon.contains(Point(m)) == False:
                new_vorpoints[c] = point_from_poly(cells[areas.index(max(areas))], 1, centered = False)[0]
                direction = [(aim[0]-new_vorpoints[c][0])*k*sign, (aim[1]-new_vorpoints[c][1])*k*sign]
                cell_dir[c] = direction

        new_vorpoints = np.array(new_vorpoints)

        if i % 4 == 0: # Every ith layer change the direction of the movement of Voronoi points
            for index, key in enumerate(cell_dir):     
                cell_dir[key] =  [-cell_dir[key][0],-cell_dir[key][1]] # Reverse the behaviour of the cells
                sign *= -1
        
        while True:   
            new_cells = []
            big = []
            small = []
            exception = []
            flag = False
            for c, newpolys in enumerate(voronoi_polygons(new_vorpoints, diameter)):
                try:
                    x, y = zip(*newpolys.intersection(boundary_polygon).exterior.coords)
                    new_poly = Polygon(list(zip(x,y)))
                    if new_poly.area > ideal_area * 10:
                        big.append(c)
                        new_cells.append(new_poly)

                        
                    elif new_poly.area < ideal_area / 10:
                        small.append(c)
                        new_cells.append(new_poly)
                    
                    else:
                        new_cells.append(new_poly)

                except: # Polygon does not intersect with the boundary polygon, split biggest cell into 2
                    exception.append(c) 

                    new_vorpoints = new_vorpoints.tolist()
                    new_vorpoints.pop(c)
                    areas = [cell.area for cell in cells]
                    fat_cell_index = areas.index(max(areas))
                    new_vorpoints.insert(c,[past_vorpoints[fat_cell_index][0] + 0.1 , past_vorpoints[fat_cell_index][1] + 0.1 ])
                    new_vorpoints = np.array(new_vorpoints)
                    cell_dir[c] = [(aim[0]-vorpoints[c][0])*k*sign, (aim[1]-vorpoints[c][1])*k*sign]
                    flag = True

            if flag == False:
                if len(new_cells) == ncells:
                    break
                                
        for cell in voronoi_polygons(new_vorpoints, diameter):
            x, y = zip(*cell.intersection(boundary_polygon).exterior.coords)
            ax.plot3D(x, y, i*thickness, c='red') # 3D plot
               
        ax.scatter(new_vorpoints[:,0], new_vorpoints[:,1], i* thickness,c='black') # 3D plot

        cells = new_cells
        past_vorpoints = new_vorpoints
        
    plt.show()
    
    return


boundary = np.array([[-50, -50], [-50, 50], [50,50], [50, -50]])

next_layer1 = np.array([[-60, -50], [-40, 60], [60,40], [60, -40]])

next_layer2 = np.array([[-70, -40], [-30, 70], [70,30], [70, -30]])

next_layer3 = np.array([[-100, -30], [-10, 90], [100,10], [100, -10]])


nb_list =[next_layer1,next_layer2,next_layer3]

ncells = 5

height = 3

thickness = 1


Main(boundary,ncells,height , thickness) # Constant Boundary
# Main(boundary,ncells,height , thickness,constant_bl=False,next_boundaries = nb_list) # Varying Boundary




