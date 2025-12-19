#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Physical constants

typedef struct {
    double sigma;
    double emissivity;
    double cw;
    double rho_water;
    double cp_air;
    double Rv;
    double R;
    double g;
    double lapse_rate;
    double M_air;
    double L_sv;
    double cs_snow;
    double L_f;
    double Rs;
} Constants;

const Constants C = {
    .sigma = 5.67e-8, //Stephan-Boltzmann constant
    .emissivity = 1.0, // of the snow
    .cw = 4186.0, // specific heat of water
    .rho_water = 1000.0, // water density
    .cp_air = 1005.0,  // specific heat of air
    .Rv = 461.5, // water vapor constant
    .R = 8.314,  // ideal gas constant
    .g = 9.81, // gravitational acceleration 
    .lapse_rate = 0.0065, //lapse rate (K/m)
    .M_air = 0.02896, // dry molare air mass (kg/mol)
    .L_sv = 2.834e6,  // latent heat of sublimation 
    .cs_snow = 2105.0,  // specific heat of snow (J/kg/K)
    .L_f = 3.34e5, // latent heat of fusion (J/kg)
    .Rs = 287.058 //  specific gas constant for dry air (J/(kg·K))
};


// defining the structure of the variables of the input file

typedef struct {
    double lat, lon;
    char time[20];  // format YYYY-MM-DDTHH:MM (19 chars + '\0')
    double temperature_2m;
    double relative_humidity_2m;
    double pressure_msl;
    double surface_pressure;
    double wind_speed_10m;
    double shortwave_radiation;
    double rain;
    double snowfall;
    double snow_depth;
} PointData;

// defining the structure of the varaibles of the output file 
typedef struct {
    double lat, lon;
    double altitude;
    char time[32];

    double snow_temp;

    // Radiations
    double Rn_polluted;
    double Rn_pure;

    // Fluxes
    double Li, Lo, H, LE, M, Qs;

    // melt enregie
    double Qm_polluted;
    double Qm_pure;

    // melt in mm/h
    double melt_polluted_mmh;
    double melt_pure_mmh;

} Result;

// grid resolution 
const double RESOLUTION = 0.0045;
const double meters_per_deg_lat = 111320.0; // meters per degree latitude

// global totals 
double total_volume_pure_m3 = 0.0;
double total_volume_polluted_m3 = 0.0;

// These are function declarations; the actual definitions appear later in the code.

void skip_header(FILE *);
int  read_point(FILE *, PointData *);

double altitude(const PointData *, const Constants *);
double Rn(double sw, double albedo);
double Li(double T, const Constants *);
double Lo(double snow_temp, const Constants *);
double es(double T);
double ea(double rh, double es);
double vapor_flux(double wind, double T, double ea, double es, const Constants *);
double air_density(double T, double P, double ea, const Constants *);
double sensible(double wind, double T, double T0, double rho, const Constants *);
double latent(double Me, const Constants *);
double rain_heat(double rain, double T, double snow_temp, const Constants *);
double snow_heat(double snowfall, double T, double snow_temp, const Constants *);
double total_melt_energy(double Rn, double Li, double Lo, double H, double LE, double M, double Qs);
double melt(double Qm, const Constants *);
double snow_temp_linear(double altitude);



//IMPLEMENTATION

// Function to skip CSV header
void skip_header(FILE *f) {
    char buf[1024];
    fgets(buf, sizeof(buf), f);
}

// Function to read a data point from CSV
int read_point(FILE *f, PointData *p) {
    char line[512];
    if (!fgets(line, sizeof(line), f)) return 0;
    return sscanf(line, "%lf,%lf,%19[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                  &p->lat, &p->lon, p->time,
                  &p->temperature_2m, &p->relative_humidity_2m,
                  &p->pressure_msl, &p->surface_pressure,
                  &p->wind_speed_10m, &p->shortwave_radiation,
                  &p->rain, &p->snowfall, &p->snow_depth) == 12;
}

// Function to calculate the altitude based on the temperature and the pressure at a point
double altitude(const PointData *p, const Constants *C) {
    double T = p->temperature_2m + 273.15;
    double exponent = -C->R * C->lapse_rate / (C->M_air * C->g);
    double ratio = p->surface_pressure / p->pressure_msl;
    return (T / C->lapse_rate) * (pow(ratio, exponent) - 1.0);
}

// Net shortwave radiation
double Rn(double sw, double albedo) { return sw * (1.0 - albedo); }

// Ingoing longwave radiation
double Li(double T, const Constants *C) {
    double Tk = T + 273.15;
    return (0.7 + 0.005*T) * C->sigma * pow(Tk, 4);
}

// Outgoing longwave radiation
double Lo(double snow_temp, const Constants *C) {
    double Ts = snow_temp + 273.15;
    return C->emissivity * C->sigma * pow(Ts, 4);
}

// Saturation vapor pressure
double es(double T) { return 0.611 * exp(17.27 * T / (T + 237.3)); }
// Actual vapor pressure
double ea(double rh, double es) { return rh * es; }

// Surface vapor transport flux
double vapor_flux(double wind, double T, double ea, double es, const Constants *C) {
    double Tk = T + 273.15;
    return 0.003 * wind * (es - ea)*1000 / (C->Rv * Tk);
}
// We used the coefficient 0.003 because using 0.03 produced unrealistic values.
// We opted to adjust the coefficient based on Mahat, Tarboton, and Molotch (2013),
// which provides more appropriate empirical coefficients for our model and parameter ranges.

// Air density
double air_density(double T, double P, double ea, const Constants *C) {
    double Tk = T + 273.15;
    return ((P*100 - ea*1000) / (C->Rs * Tk)) + ((ea*1000) / (C->Rv * Tk));
}

// Sensible heat flux
double sensible(double wind, double T, double T0, double rho, const Constants *C) {
    return 0.003 * wind * rho * C->cp_air * (T - T0);
}

// Latent heat flux
double latent(double Me, const Constants *C) { return C->L_sv * Me; }

// Advected heat from rain
double rain_heat(double rain, double T, double snow_temp, const Constants *C) {
    double R = rain * 0.001 / 3600.0;
    double T_rain = (T > snow_temp ? T : snow_temp);
    return R * C->rho_water * C->cw * (T_rain - snow_temp);
}

// Advected heat from snow precipitation
double snow_heat(double snowfall, double T, double snow_temp, const Constants *C) {
    double S = snowfall * 0.001 / 3600.0;
    double Ts = (T < snow_temp ? T : snow_temp);
    return S * C->rho_water * C->cs_snow * (Ts - snow_temp);
}

// Global energy balance
double total_melt_energy(double Rn, double Li, double Lo, double H, double LE, double M, double Qs) {
    return Rn + Li - Lo + H + LE + M + Qs;
}

// Meltrate speed of snow
double melt(double Qm, const Constants *C) {
    if (Qm <= 0) return 0;
    double melt_m_per_s = Qm / (C->rho_water * C->L_f);
    return melt_m_per_s * 3.6e6;   // conversion m/s to mm/h
}

// Linear snow temperature approximation based on altitude
double snow_temp_linear(double altitude) {

    const double z_min  = 818.24;   // minimum altitude of the dataset
    const double z_max  = 4747.62;  // maximum altitude of the dataset(the summit of Mont Blanc)
    const double Ts_min = 0.0;      // T snow at z_min
    const double Ts_max = -12.0;    // T snow at z_max (begining of may)

    // clamp on borders
    if (altitude <= z_min) return Ts_min;
    if (altitude >= z_max) return Ts_max;

    // linear interpolation
    double frac = (altitude - z_min) / (z_max - z_min);
    return Ts_min + frac * (Ts_max - Ts_min);
}

// MAIN 

int main(int argc, char **argv) {
    // Check if the user provided an input CSV file as a command-line argument

    printf("=== Snow Energy Balance Model ===\n");

    FILE *in = fopen("../data/mont_blanc_weather_2025_day.csv", "r");
    if (!in) {
        perror("Opening error CSV"); // Print an error if the file cannot be opened
        return 1;
    }

    printf("Input file opened.\n");

    // Create results folder if it doesn't exist
    const char *resultsFolder = "../results";

    //  Checking for errors when creating the result file
    FILE *out = fopen("../results/results_melt.csv", "w");
    if (!out) {
        perror("Creation error of results_melt.csv"); // Print an error if the file cannot be created
        fclose(in); // Close input file before exiting
        return 1;
    }   

    printf("Output file results_melt.csv created.\n");

    // Creation of result file header
    fprintf(out,
        "lat,lon,time,altitude,snow_temp,"
        "Rn_polluted,Rn_pure,"
        "Li,Lo,H,LE,M,Qs,"
        "Qm_polluted,Qm_pure,"
        "melt_polluted_mmh,melt_pure_mmh\n");
    

    skip_header(in);

    // Calculation of all the physical functions for every point at every time
    printf("Starting point-by-point calculations...\n");
    PointData p;
    int count = 0;
    while (read_point(in, &p)) {

        Result r;
        r.lat = p.lat;
        r.lon = p.lon;
        strcpy(r.time, p.time);

        double esv = es(p.temperature_2m);
        double eav = ea(p.relative_humidity_2m, esv);
        double rho = air_density(p.temperature_2m, p.surface_pressure, eav, &C);
        double Mev = vapor_flux(p.wind_speed_10m, p.temperature_2m, eav, esv, &C);

        r.altitude = altitude(&p, &C);

        // Computing snow temperature dynamically based on altitude
        double snow_temp = snow_temp_linear(r.altitude);
        r.snow_temp = snow_temp;

        // Shortwave net radiation for polluted and pure snow (different albedos)
        r.Rn_polluted = Rn(p.shortwave_radiation, 0.7);
        r.Rn_pure    = Rn(p.shortwave_radiation, 0.85);

        // Common fluxes
        r.Li = Li(p.temperature_2m, &C);
        r.Lo = Lo(snow_temp, &C);
        r.H  = sensible(p.wind_speed_10m, p.temperature_2m, snow_temp, rho, &C);
        r.LE = latent(Mev, &C);
        r.M  = rain_heat(p.rain, p.temperature_2m, snow_temp, &C);
        r.Qs = snow_heat(p.snowfall, p.temperature_2m, snow_temp, &C);

        // Melt energy
        r.Qm_polluted = total_melt_energy(r.Rn_polluted, r.Li, r.Lo, r.H, r.LE, r.M, r.Qs);
        r.Qm_pure    = total_melt_energy(r.Rn_pure,    r.Li, r.Lo, r.H, r.LE, r.M, r.Qs);

        // Melt (mm/h)
        r.melt_polluted_mmh = melt(r.Qm_polluted, &C);
        r.melt_pure_mmh    = melt(r.Qm_pure,    &C);

        // Grid Cell Area
        double lat_rad = p.lat * M_PI / 180.0;
        double meters_per_deg_lon = meters_per_deg_lat * cos(lat_rad);
        double area_point_m2 = RESOLUTION * meters_per_deg_lat * RESOLUTION * meters_per_deg_lon;

        // Volume calculation
        double dt_h = 1.0; // hourly timestep
        double vol_pure_m3    = (r.melt_pure_mmh    / 1000.0) * dt_h * area_point_m2;
        double vol_polluted_m3 = (r.melt_polluted_mmh / 1000.0) * dt_h * area_point_m2;

        total_volume_pure_m3    += vol_pure_m3;
        total_volume_polluted_m3 += vol_polluted_m3;

        // Write CSV line
        fprintf(out,
            "%lf,%lf,%s,%lf,%lf,"
            "%lf,%lf,"
            "%lf,%lf,%lf,%lf,%lf,%lf,"
            "%lf,%lf,"
            "%lf,%lf\n",
        
            r.lat, r.lon, r.time, r.altitude, r.snow_temp,
            r.Rn_polluted, r.Rn_pure,
            r.Li, r.Lo, r.H, r.LE, r.M, r.Qs,
            r.Qm_polluted, r.Qm_pure,
            r.melt_polluted_mmh, r.melt_pure_mmh
        );
        

        count++;
    }

    printf("Finished ! %d treated lines.\n", count);

    fclose(in);
    fclose(out);

    printf("Closed input and output files.\n");

    // Percentage increase 
    double vol_diff = total_volume_polluted_m3 - total_volume_pure_m3;
    double pct_increase = (vol_diff / total_volume_pure_m3) * 100.0;
    // Compare to olympic pool
    const int OLYMPIC_POOL_VOLUME = 2500; // m³
    double num_pools = vol_diff / OLYMPIC_POOL_VOLUME;
    
    
    // Prepare LaTeX file path
    char texPath[512];
    snprintf(texPath, sizeof(texPath), "%s/melt_volume_results.tex", resultsFolder);

   // Write global snow melt volume results to a LaTeX file

    FILE *tex = fopen(texPath, "w");
    if (tex == NULL) {
     perror("Error while creating melt_volume_results.tex");
    } else {

        fprintf(tex,
            "\\documentclass{article}\n"
            "\\usepackage{siunitx}\n"
            "\\usepackage{geometry}\n"
            "\\geometry{margin=2.5cm}\n"
            "\\begin{document}\n\n"

            "\\section*{Integrated Snow Melt Volumes}\n\n"

            "This document summarizes the total snow melt volumes obtained from the "
            "energy balance model applied over the Mont-Blanc area during the first "
            "week of May 2025. Two surface conditions are compared: fresh snow and "
            "polluted snow with reduced albedo.\n\n"


            "\\begin{itemize}\n"
            "  \\item Total melted volume (fresh snow): "
            "\\SI{%.2f}{\\cubic\\meter}\n"
            "  \\item Total melted volume (polluted snow): "
            "\\SI{%.2f}{\\cubic\\meter}\n"
            "  \\item Additional melt due to pollution: "
            "\\SI{%.2f}{\\cubic\\meter}\n"
            "\\end{itemize}\n\n"

            "The excess melt caused by snow pollution corresponds to "
            "\\textbf{%.2f Olympic-size swimming pools}, assuming a reference volume "
            "of \\SI{2500}{\\cubic\\meter} per pool.\n\n"

            "Overall, the total snow melt is \\textbf{%.2f\\%% higher} under polluted "
            "snow conditions compared to fresh snow.\n\n"

            "\\end{document}\n",
            total_volume_pure_m3,
            total_volume_polluted_m3,
            vol_diff,
            num_pools,
            pct_increase
        );

    fclose(tex);
    printf("melt_volume_results.tex completed.\n");
    }

    printf("=== Model finished successfully ===\n");

    return 0;
}
