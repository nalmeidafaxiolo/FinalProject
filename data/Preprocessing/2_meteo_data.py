import csv
import requests
import time
import os

# File paths and API configuration
GRID_FILE = "mont_blanc_grid_2160.csv"       # CSV file with latitude-longitude grid points
OUTPUT_FILE = "mont_blanc_weather_2025.csv"  # CSV file where we will save the weather data
API_BASE = "https://archive-api.open-meteo.com/v1/archive"  # Open-Meteo historical data API endpoint

# Dates for the period we want to retrieve weather data
START_DATE = "2025-05-01"
END_DATE   = "2025-05-07"


# Hourly weather variables to retrieve from API
HOURLY_VARS = ",".join([
    "temperature_2m",       # Temperature at 2 meters above ground (°C)
    "relative_humidity_2m", # Relative humidity at 2 meters (%)
    "pressure_msl",         # Mean sea-level atmospheric pressure (hPa)
    "surface_pressure",     # Surface pressure at the point (hPa)
    "wind_speed_10m",       # Wind speed at 10 meters above ground (m/s)
    "shortwave_radiation",  # Incoming shortwave radiation (W/m²)
    "rain",                 # Rainfall in mm
    "snowfall",             # Snowfall in mm
    "snow_depth",           # Snow depth in meters
    "is_day"                # Binary indicator whether it is daytime
])

PAUSE = 0.25  # Pause in seconds between API requests to avoid overloading the server


# Function: load_grid
# Purpose: read the grid CSV and return a list of (latitude, longitude) tuples
def load_grid(filename):
    grid = []  # Empty list to store all points
    with open(filename, newline="") as f:
        reader = csv.DictReader(f)  # Use DictReader to get named columns
        for row in reader:
            # Convert the string values to float and store as a tuple
            grid.append((float(row["lat"]), float(row["lon"])))
    return grid  # Return the complete list of points


# Function: fetch_point
# Purpose: request historical weather data for a single grid point
# - Uses the Open-Meteo API
# - Implements retry logic if too many requests are made (HTTP 429)
# - Returns JSON data or None if there is an error
def fetch_point(lat, lon):
    # Build API request parameters
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": START_DATE,
        "end_date": END_DATE,
        "hourly": HOURLY_VARS,
        "timezone": "UTC"  # All times are returned in UTC
    }

    try:
        # Send GET request to the API
        r = requests.get(API_BASE, params=params, timeout=15)

        # If we get a "Too Many Requests" response, pause and retry
        if r.status_code == 429:
            print("429 Too Many Requests → pausing 5 seconds…")
            time.sleep(5)
            return fetch_point(lat, lon)  # Recursive retry

        # If another error occurs, print first 100 chars of the response
        if r.status_code != 200:
            print(f"API Error {r.status_code}: {r.text[:100]}")
            return None

        return r.json()  # Return the API response as JSON

    except Exception as e:
        # Catch any network or unexpected errors
        print(f"Network exception: {e}")
        return None


# Function: save_point
# Purpose: save hourly data for a single grid point to the CSV
# - Takes the CSV writer, latitude, longitude, and JSON data
# - Loops over each hour and writes a row with all variables
def save_point(writer, lat, lon, data):
    hours = data["hourly"]["time"]  # List of timestamps for each hour
    for i in range(len(hours)):
        # Initialize a row with coordinates and timestamp
        row = {
            "lat": lat,
            "lon": lon,
            "time": hours[i]
        }
        # Add all other variables for this hour
        for var in data["hourly"]:
            if var != "time":  # Skip the time key since it is already included
                row[var] = data["hourly"][var][i]
        writer.writerow(row)  # Write the row to the CSV


# Main workflow
# - Loads grid points
# - Fetches weather data for each point
# - Saves data to CSV
# - Supports resuming if process was previously interrupted
def main():
    # Step 1: load grid points
    print("Loading grid...")
    grid = load_grid(GRID_FILE)
    print(f"➡️ {len(grid)} points found.")

    start_idx = 0          # Index to start processing grid points
    processed_last = False # Flag indicating whether we are resuming

    # Step 2: check if output CSV already exists (resume functionality)
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, newline="") as f:
            reader = list(csv.DictReader(f))
            if reader:
                # Find the last processed point
                last_lat = float(reader[-1]["lat"])
                last_lon = float(reader[-1]["lon"])
                print(f"Resuming from {last_lat}, {last_lon}")

                for i, (plat, plon) in enumerate(grid):
                    if plat == last_lat and plon == last_lon:
                        start_idx = i + 1  # Start after last processed point
                        processed_last = True
                        break

    # Step 3: set file mode: append if resuming, write if starting fresh
    write_mode = "a" if processed_last else "w"

    # Step 4: open CSV for writing data
    with open(OUTPUT_FILE, write_mode, newline="") as f:
        header = ["lat", "lon", "time"] + HOURLY_VARS.split(",")
        writer = csv.DictWriter(f, fieldnames=header)
        if write_mode == "w":  # Only write header if starting fresh
            writer.writeheader()

        print(f"Starting at point {start_idx}/{len(grid)}")

        try:
            # Step 5: loop through all grid points
            for idx in range(start_idx, len(grid)):
                lat, lon = grid[idx]
                print(f"{idx+1}/{len(grid)} → lat={lat}, lon={lon}")

                # Fetch data from API
                data = fetch_point(lat, lon)
                if data:
                    # Save data to CSV
                    save_point(writer, lat, lon, data)

                # Pause briefly to avoid overloading API
                time.sleep(PAUSE)

        except KeyboardInterrupt:
            # Allow manual interruption and resume later
            print("\n Manual stop, resume possible")
            return

    # Step 6: end of process
    print("Completed: file generated →", OUTPUT_FILE)


# Standard Python idiom for running the script directly

if __name__ == "__main__":
    main()
