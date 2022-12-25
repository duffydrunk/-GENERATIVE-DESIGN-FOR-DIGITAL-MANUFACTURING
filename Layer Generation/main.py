# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 12:48:56 2021

@author: Ege Uğur Aguş
"""
import matplotlib.pyplot as plt
import matplotlib.colors as color
import numpy as np

def file_to_list(path):
    '''
    Parameters
    ----------
    path : str
        Path of the file which a list is going to be created from.

    Returns
    -------
    A numpy array that contains float values of the source file.

    '''
    file = open(path, "r") 
    text = file.read() 
    file.close() 
    text = text.replace("\n", ",") 
    text = text.replace("\t", ",")
    text = text.replace(" ", ",")
    
    return np.fromstring(text, dtype=float, sep=',')

def intersection_point(line1, line2):
    '''
    Parameters
    ----------
    line1 : tuple
        Contains the coordinates of 2 two points from a line
    line2 : tuple
        Contains the coordinates of 2 two points from a line    
    Returns
    -------
    Coordinates of intersection point of the line1 and line2
    by using the formula given at
    "https://mathworld.wolfram.com/Line-LineIntersection.html"
    '''
    x_diff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    y_diff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(x, y):
        return x[0] * y[1] - x[1] * y[0]

    denominator = det(x_diff, y_diff)
    if denominator == 0:
        return False,False

    c = (det(*line1), det(*line2)) 
    x = det(c, x_diff) / denominator
    y = det(c, y_diff) / denominator
    return x, y

def intersection_at_edge(intersection,corner1,corner2):
    '''
    Parameters
    ----------
    intersection : tuple
        Contains the coordinates of a point
    corner1 : tuple
        Contains the coordinates of a point    
    corner2 : tuple
        Contains the coordinates of a point
    Returns
    -------
    True if intersection point is on the line segment of
    corner1 to corner2
    '''
    x1,y1 = corner1         
    x2,y2 = corner2         
    x3,y3 = intersection
    
    if type(x3) == bool or type(y3) == bool:
        return False
    
    corner_to_corner = np.sqrt((x2-x1)**2+(y2-y1)**2)
    corner1_to_intersection = np.sqrt((x3-x1)**2 + (y3-y1)**2)
    corner2_to_intersection = np.sqrt((x2-x3)**2+(y2-y3)**2)
    
    th = 5e-10

    if th >= abs (corner_to_corner - corner1_to_intersection - corner2_to_intersection):
        return True
    else:
        return False
    
def Main(path_bl, path_lines, layer_thickness, number_of_layers):
    '''
    Parameters
    ----------
    path_bl : str
        Path or name of the file that contains the vertices of the boundary polygon
    path_lines : str
        Path or name of the file that contains the information of the line segments 
    layer_thickness : int
        Layer thickness value in mm
    number_of_layers: int
        Number of layers value
    Returns
    -------
    None. Displays the object constructed layer by layer.
    '''
    current_layer = 0

    ax = plt.axes(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    bl_coordinates = file_to_list(path_bl)
    if len(bl_coordinates) == 0:                                            #Check given input file is not empty
        return print("Vertices of the boundary polygon are not defined.")
    
    bl_corners = np.array_split(bl_coordinates, len(bl_coordinates)/2)
    if len(bl_corners) < 3 :                                                #Check if given corner number is larger than 2
        return print("The amount of corners defined is not enough.")
    bl_corners.append(bl_corners[0])
    bl_coordinates = np.append(bl_coordinates,[bl_coordinates[0],bl_coordinates[1]])
    extreme = bl_coordinates.max()*10
    
    line_coordinates = file_to_list(path_lines)
    if len(line_coordinates) == 0:
        for i in range(number_of_layers):
            ax.plot3D(bl_coordinates[::2],bl_coordinates[1::2],current_layer,c='black')
            current_layer += layer_thickness
        return
    
    line_coordinates = np.array_split(line_coordinates, len(line_coordinates)/5)
      
    for k in range(len(line_coordinates)):
        spacing = 1/(len(line_coordinates)) 
        line_color = np.array([[0+(spacing)*k,1,1]])  # Define a distinct value for every line in HSV
        line_color = color.hsv_to_rgb(line_color) # Covnert color to RGB from HSV
        line_coordinates[k] = np.append(line_coordinates[k],line_color) # Add information to the specific line 
    
    for i in range(number_of_layers):
        ax.plot3D(bl_coordinates[::2],bl_coordinates[1::2],current_layer,c='black') # Select BL color as Black since our line's never going to 
                                                                                    #get black as color
        for line in line_coordinates:
            intersection_points =[]
            x_founded = []
            y_founded = []
            draw = line.copy()
            line_color, draw = draw[-3:], draw[:-3]
            
            for i in range(1,len(bl_corners)):
                line1 = ((bl_corners[i-1][0],bl_corners[i-1][1]),(bl_corners[i][0],bl_corners[i][1]))
                line2 = ((draw[0],draw[1]),(draw[2],draw[3]))
                x,y = intersection_point(line1,line2)
                intersection_points.append((x,y))               
        
            intersection_points.append(intersection_points[0])
            for i in range(1,len(intersection_points)):
                if intersection_at_edge(intersection_points[i-1],bl_corners[i-1],bl_corners[i]):
                    x_founded.append(intersection_points[i-1][0])
                    y_founded.append(intersection_points[i-1][1])
                    
            if len(x_founded) > 2:
                x_founded , y_founded = list(map(list, zip(*list(dict.fromkeys(list(zip(x_founded,y_founded)))))))
                arrange_values = list(zip(x_founded,y_founded))
                arrange_values.sort(key = lambda x:x[0])
                x_founded , y_founded = list(map(list,zip(*arrange_values)))                
            
            if len(x_founded) == 2:
                if any(([x_founded[0],y_founded[0]] == x).all() for x in bl_corners) and \
                    any(([x_founded[1],y_founded[1]] == x).all() for x in bl_corners):
                    pass
                else:
                    ax.plot3D(x_founded,y_founded,current_layer,c=line_color)

            elif len(x_founded) > 2:
                x_founded.append(x_founded[0])
                y_founded.append(y_founded[0])
                for counter1 in range(1,len(x_founded)-1):
                    multiple_segment_points = []
                    multiple_x_points = []
                    multiple_y_points = []
                    mid_point_coor = ((x_founded[counter1-1]+x_founded[counter1])/2,(y_founded[counter1-1]+y_founded[counter1])/2)
                    line3 = ((mid_point_coor[0],mid_point_coor[1]),(extreme,mid_point_coor[1]))
                    
                    for j in range(1,len(bl_corners)):
                        line4 = ((bl_corners[j-1][0],bl_corners[j-1][1]),(bl_corners[j][0],bl_corners[j][1]))
                        x,y = intersection_point(line3,line4)      
                        multiple_segment_points.append((x,y))   

                    multiple_segment_points.append([multiple_segment_points[0]])
                    for m in range(1,len(multiple_segment_points)):
                        if intersection_at_edge(multiple_segment_points[m-1],bl_corners[m-1],bl_corners[m]):
                            multiple_x_points.append(multiple_segment_points[m-1][0])
                            multiple_y_points.append(multiple_segment_points[m-1][1])  
                    
                    points_at_right = [x for x in multiple_x_points if x > mid_point_coor[0]] 
                    if len(points_at_right) % 2 == 1: # Mid point is inside the polygon
                        ax.plot3D([x_founded[counter1-1],x_founded[counter1]],[y_founded[counter1-1],y_founded[counter1]],current_layer,c=line_color)                                                   

            teta  = np.arctan2(line[3]-line[1],line[2]-line[0])
            y_displacement = np.cos(teta)*(line[4]*layer_thickness)
            x_displacement = np.sin(teta)*(line[4]*layer_thickness)
            for i in range(2):
                line[i*2] += x_displacement
                line[i*2+1] -= y_displacement

        current_layer += layer_thickness
 
