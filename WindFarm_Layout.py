import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.spatial.distance import cdist


def create_parallelogram(vertices):
    """
    Create the vertices of a parallelogram based on given dimensions.

    Parameters:
    vertices (list of tuples): List of vertices (x, y)

    Returns:
    np.array: Vertices of the parallelogram
    """
    return np.array(vertices)

def is_point_in_parallelogram(point, parallelogram):
    """
    Check if a given point is inside or on the border of the parallelogram.

    Parameters:
    point (tuple): Point to check (x, y)
    parallelogram (np.array): Vertices of the parallelogram

    Returns:
    bool: True if the point is inside or on the border of the parallelogram, False otherwise
    """
    from matplotlib.path import Path
    path = Path(parallelogram)
    return path.contains_point(point) or path.contains_point(point, radius=-1e-9)

def distance_point_to_line(point, line_start, line_end):
    """
    Calculate the perpendicular distance from a point to a line defined by two points.

    Parameters:
    point (tuple): Coordinates of the point (x, y)
    line_start (tuple): Coordinates of the start point of the line (x, y)
    line_end (tuple): Coordinates of the end point of the line (x, y)

    Returns:
    float: Perpendicular distance from the point to the line
    """
    # px, py = point
    # x1, y1 = line_start
    # x2, y2 = line_end
    # num = abs((y2 - y1) * px - (x2 - x1) * py + x2 * y1 - y2 * x1)
    # den = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    
    px, py = point
    angle = 65.3
    # Defining equation of first line y = dx + e
    d = np.tan(np.radians(angle))
    e = py - d * px

    # Defining equation of second line y =mx + c
    x1, y1 = line_start
    x2, y2 = line_end
    
    m = (y2 -y1)/(x2-x1)
    c = y1 - m *x1

    # Equating both lines - dx + e = mx + c
    x = (c-e)/(d-m)
    y = m*x +c
    final = [ x, y ]
    distance = math.dist(point, final)

    return distance

def generate_points(parallelogram, t_dist):
    """
    Generate points within the parallelogram.

    Parameters:
    parallelogram (np.array): Vertices of the parallelogram
    t_dist (float): Distance used to divide edges to determine number of points

    Returns:
    tuple: x coordinates of the points, y coordinates of the points, FirstColPoints, FirstRowPoints
    """
    x_coords = []
    y_coords = []

    A, B, C, D = parallelogram

    NbrPoints = 0

    # Calculate FirstColPoints
    distance_AB = np.linalg.norm(B - A)
    FirstColPoints = int(distance_AB // t_dist)


    # Calculate FirstRowPoints
    distance_AD = np.linalg.norm(D - A)
    FirstRowPoints = int(distance_AD // t_dist)
    x_spacing = distance_AD/(FirstRowPoints)


    # Generate points along the line A-B
    for i in range(FirstColPoints +1):
        t = i / (FirstColPoints)
        x = (1 - t) * A[0] + t * B[0]
        y = (1 - t) * A[1] + t * B[1]
        x_coords.append(x)
        y_coords.append(y)

        NbrPoints += 1

    # Generate points for the next columns
    ChangeinY = np.tan(np.radians(-19)) * t_dist
    NumberToNew = -int(t_dist // ChangeinY) + 1

    counter = 1 # Counting number of points in column
    next = 0 # Switch to decrease number of points

    for j in range(2,FirstRowPoints+1):
        # print("Calculating col. ", j)
        # Testing if we decrease the amount of points in the column
        if j % NumberToNew == 0:
            # WE DECREASE ON THE NEXT
            
            next = 1
            
            # Calculate the distance from the first point of the column to the line B-C
            first_point_column = (A[0] + (j-1)*x_spacing, A[1])
            distance_to_BC = distance_point_to_line(first_point_column, B, C)
            for i in range(FirstColPoints - counter +1):
                # t = i / (FirstColPoints - 1)
                # x = x_coords[i] + (j-1) * t_dist

                # y_spacing = distance_to_BC/(FirstRowPoints - counter)
                # y =  A[1] + (i) * y_spacing 

                if i == 0: # Defining first point of each col
                    
                    x = A[0] + (j-1) * x_spacing
                    y = A[1]

                else:
                    # y_spacing = distance_to_BC/(FirstRowPoints - counter)
                    # y =  A[1] + (i) * y_spacing 
                    
                    y_spacing = distance_to_BC/(FirstColPoints - counter)
                    
                    x = (j-1) * x_spacing + i * np.cos(np.radians(65.3)) * (y_spacing)

                    y =  np.sin(np.radians(65.3)) * (y_spacing * (i))

                if is_point_in_parallelogram((x, y), parallelogram):
                    x_coords.append(x)
                    y_coords.append(y)
                    NbrPoints += 1

        elif next == 1:
            # WE DECREASE NOW

            # Calculate the distance from the first point of the column to the line B-C
            first_point_column = (A[0] + (j-1)*x_spacing, A[1])
            distance_to_BC = distance_point_to_line(first_point_column, B, C)
            for i in range(FirstColPoints  - counter + 1):
                # t = i / (FirstColPoints - 1)

                if i == 0: # Defining first point of each col
                    x = (j-1) * x_spacing
                    y = 0

                else:
                    # y_spacing = distance_to_BC/(FirstRowPoints - counter)
                    # y =  A[1] + (i) * y_spacing 
                    
                    y_spacing = distance_to_BC/(FirstColPoints - counter)

                    x = (j-1) * x_spacing + i * np.cos(np.radians(65.3)) * (y_spacing)

                    y =  np.sin(np.radians(65.3)) * (y_spacing * (i))

                if is_point_in_parallelogram((x, y), parallelogram):
                    x_coords.append(x)
                    y_coords.append(y)
                    NbrPoints += 1

            counter += 1 # Decrease the number of points by one for the next column
            next = 0 # Reset back to first column after decrease in points
            
        else: 
            
            # Calculate the distance from the first point of the column to the line B-C
            first_point_column = (A[0] + (j-1)*x_spacing, A[1])
            distance_to_BC = distance_point_to_line(first_point_column, B, C)
            for i in range(FirstColPoints  - counter +1):
                # t = i / (FirstColPoints - 1)
                # x = x_coords[i] + (j-1) * t_dist
                if i == 0: # Defining first point of each col
                    x = (j-1) * x_spacing
                    y = 0

                else:
                    # y_spacing = distance_to_BC/(FirstRowPoints - counter)
                    # y =  A[1] + (i) * y_spacing 
                    
                    y_spacing = distance_to_BC/(FirstColPoints - counter)
                    
                    x = (j-1) * x_spacing + i * np.cos(np.radians(65.3)) * (y_spacing)

                    y =  np.sin(np.radians(65.3)) * (y_spacing * (i))


                if is_point_in_parallelogram((x, y), parallelogram):
                    x_coords.append(x)
                    y_coords.append(y)
                    NbrPoints += 1
        

    # Calculate FirstColPoints
    distance_CD = np.linalg.norm(C - D)
    LastColPoints = int(distance_CD // t_dist)

    # Generate points along the line C-D
    for i in range(LastColPoints):
        t = i / (LastColPoints - 1)
        x = (1 - t) * C[0] + t * D[0]
        y = (1 - t) * C[1] + t * D[1]
        x_coords.append(x)
        y_coords.append(y)

        NbrPoints += 1
    
    print("Total Points: ", NbrPoints)


    return x_coords, y_coords, FirstColPoints, FirstRowPoints, distance_to_BC, NbrPoints

def plot_parallelogram_and_points(parallelogram, x_coords, y_coords):
    """
    Plot the parallelogram and the points within it.

    Parameters:
    parallelogram (np.array): Vertices of the parallelogram
    x_coords (list): x coordinates of the points
    y_coords (list): y coordinates of the points
    """
    # Plot the parallelogram
    plt.figure(figsize=(10, 6))
    parallelogram_path = np.vstack([parallelogram, parallelogram[0]])  # Close the loop
    plt.fill(parallelogram_path[:, 0], parallelogram_path[:, 1], 'b', alpha=0.5)
    

    # Plot the points
    plt.scatter(x_coords, y_coords, c='red', marker='o')

    # Add North reference
    plt.annotate('N', xy=(2, 25), xytext=(2, 23),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=4, headlength=6),
                 fontsize=12, ha='center', va='bottom')
    
    plt.title("Windfarm Optimised Layout - Lagelander")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.show()


# def filter_points(x_coords, y_coords, t_dist):
#     """
#     Filter points to remove those that are closer than t_dist to the previous point.

#     Parameters:
#     x_coords (list): x coordinates of the points
#     y_coords (list): y coordinates of the points
#     t_dist (float): Minimum distance between points

#     Returns:
#     tuple: Filtered x and y coordinates
#     """
#     filtered_x_coords = []
#     filtered_y_coords = []

#     for i in range(0, len(x_coords)):
#         final = [x_coords[i],y_coords[i]]
#         point = [x_coords[i-1],y_coords[i-1]]


#         distance = np.abs(math.dist(point, final))
#         # distance = np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
#         if distance >= t_dist:
#             filtered_x_coords.append(x_coords[i])
#             filtered_y_coords.append(y_coords[i])

#     return filtered_x_coords, filtered_y_coords


def filter_points(x_coords, y_coords, t_dist):
    """
    Filter points to remove those that are closer than t_dist to the previous point.

    Parameters:
    x_coords (list): x coordinates of the points
    y_coords (list): y coordinates of the points
    t_dist (float): Minimum distance between points

    Returns:
    tuple: Filtered x and y coordinates
    """
    coords = [x_coords,y_coords]

    coords = np.column_stack((x_coords, y_coords))

    for i in range(0, len(x_coords)):
        point = np.array([x_coords[i],y_coords[i]])   
        # coords[(cdist(coords[:,:],point[None]) > t_dist).ravel()]
        coords = coords[np.linalg.norm(coords[:,:] - point, axis=1) > t_dist]
    
    filtered_x_coords = coords[:, 0]
    filtered_y_coords = coords[:, 1]

    return filtered_x_coords, filtered_y_coords


# Define the vertices of the parallelogram based on given dimensions and angles
A = np.array([0, 0])
B = np.array([27.5 * np.cos(np.radians(65.3)), 27.5 * np.sin(np.radians(65.3))])
C = np.array([B[0] + 43.9 * np.cos(np.radians(-19)), B[1] + 43.9 * np.sin(np.radians(-19))])
D = np.array([A[0] + 46.7, A[1]])
parallelogram = create_parallelogram([A, B, C, D])

def calculate_parallelogram_area(A,B,C,D):
    first_triangle = np.abs(1/2*B[0]*B[1])
    second_triangle = np.abs(1/2*(C[0]-B[0])*(C[1]-B[1]))
    third_triangle = np.abs(1/2*(D[0]-C[0])*(D[1]-C[1]))
    square = np.abs(C[1]*(D[0]-B[0]))

    area = first_triangle + second_triangle + third_triangle + square

    return area

# Input function for t_dist
t_dist = float(input("Enter the value for t_dist: "))

# Generate points
x_coords, y_coords, FirstColPoints, FirstRowPoints, distance_to_BC, NbrPoints = generate_points(parallelogram, t_dist)

# print(x_coords, y_coords, len(x_coords), len(y_coords))
# Filter points based on the distance
filtered_x_coords, filtered_y_coords = filter_points(x_coords, y_coords, t_dist)

area_parallelogram = calculate_parallelogram_area(A,B,C,D)
# Print the number of points
print(f"FirstColPoints: {FirstColPoints}")
print(f"FirstRowPoints: {FirstRowPoints}")
print("Total Area of Parellelogram: ", area_parallelogram)
print("Energy Density: ", NbrPoints*30/area_parallelogram)

# Plot the parallelogram and the points within it
plot_parallelogram_and_points(parallelogram, x_coords, y_coords)

# Plot the parallelogram and the points within it
# plot_parallelogram_and_points(parallelogram, filtered_x_coords, filtered_y_coords)
