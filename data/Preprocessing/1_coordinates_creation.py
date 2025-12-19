import csv

# Definition of the spatial domain.
# The objective is to build a regular grid covering the Mont-Blanc area.
# Latitude and longitude boundaries define a rectangular region, and the
# resolution determines the spacing between grid points.
LAT_MIN, LAT_MAX = 45.7426, 45.9226
LON_MIN, LON_MAX = 6.7452, 6.9852
RESOLUTION = 0.0045  # Degrees. This yields 2160 points.


# The program generates a CSV file containing every (lat, lon) coordinate
# of the grid, sampled at the specified resolution. The output is intended
# for use in subsequent meteorological computations.
with open("mont_blanc_grid_2160.csv", "w", newline="") as f:
    writer = csv.writer(f)
    
    # Write header row for clarity and downstream usability
    writer.writerow(["lat", "lon"])
    
    # Iteration over latitude and longitude in nested loops
    # to generate the full grid.
    lat = LAT_MIN
    while lat <= LAT_MAX:
        lon = LON_MIN
        while lon <= LON_MAX:
            # Each coordinate is rounded for compactness and file readability
            writer.writerow([round(lat, 6), round(lon, 6)])
            lon += RESOLUTION
        lat += RESOLUTION

print("File mont_blanc_grid_2160.csv created !")
