#include <iostream>
#include <string>
#include <array>
#include <windows.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>

extern "C" 
{
    #include "SpiceUsr.h"
}

// kilometers per astronomic unit
double AU = 149597870.7;

// I wrote these backwards once and spent hours 
// trying to debug what was going wrong...
double rad2deg(double x)
{
    return x * 180 / pi_c();
}

double deg2rad(double x)
{
    return x * pi_c() / 180;
}

// --- --- --- STRUCT AND CLASS DEFINITIONS --- --- --- [START]
struct Vec3 // extremely self-explanatory
{
    double x, y, z;

    Vec3()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vec3(double xp, double yp, double zp)
    {
        x = xp;
        y = yp;
        z = zp;
    }

    Vec3(std::vector<double> vec)
    {
        x = vec[0];
        y = vec[1];
        z = vec[2];
    }

    Vec3(std::array<double, 3Ui64> vec)
    {
        x = vec[0];
        y = vec[1];
        z = vec[2];
    }

    Vec3(double RA, double DEC) // this is equivalent to a spherical2cartezian() function
    {
        x = cos(deg2rad(DEC)) * cos(deg2rad(RA));
        y = cos(deg2rad(DEC)) * sin(deg2rad(RA));
        z = sin(deg2rad(DEC));
    }

    Vec3 operator+(const Vec3& other) const 
    {
        return { x + other.x, y + other.y, z + other.z };
    }

    Vec3 operator-(const Vec3& other) const
    {
        return { x - other.x, y - other.y, z - other.z };
    }

    Vec3 operator*(double scalar) const 
    {
        return { x * scalar, y * scalar, z * scalar };
    }

    Vec3 operator/(double scalar) const
    {
        return { x / scalar, y / scalar, z / scalar };
    }

    Vec3& operator+=(const Vec3& other) 
    {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vec3 cross(const Vec3& other)
    {
        return Vec3(y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x);
    }

    double dot(const Vec3& other)
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 normalized()
    {
        return Vec3(x, y, z) / Vec3(x, y, z).mag();
    }

    double mag()
    {
        return sqrt(x * x + y * y + z * z);
    }

    void printout()
    {
        std::cout << "Vec3(" << x << ", " << y << ", " << z << ")\n";
    }
};

// represents a physical gravitational field generating body (Sun + planets usually)
// used for orbit propagation
struct Body 
{
    std::string name;
    double GM;
    Vec3 pos;

    Body(std::string namep, double GMp)
    {
        name = namep;
        GM = GMp;
        pos = Vec3();
    }
};

// represents an object which does not generate a notable gravitational field and is
// propagated via Yoshida8 instead of retrieving through SPICE
struct MinorPlanet 
{
    Vec3 pos;
    Vec3 vel;

    MinorPlanet(Vec3 p, Vec3 v)
    {
        pos = p;
        vel = v;
    }
};

// ground observatories with MPC obscodes
class Observatory {
public:
    std::string name;
    std::string code;
    double l; // long
    double s; // sin
    double c; // cos

    Observatory()
    {
        name = "None";
        code = "---";
        l = 0;
        s = 0;
        c = 1;
    }

    Observatory(const std::string& name, const std::string& code,
        double l, double c, double s) : name(name), code(code), l(l), s(s), c(c) {}

    std::string toString() const
    {
        return code + " (" + name + ")";
    }
};

// Obs80 (MPC 80-column format) observations
class Observation
{
public:
    std::string obs_str;
    std::string perm;
    std::string prov;
    SpiceDouble et;
    double RA;
    double DEC;
    double mag;
    std::string obs_code;

    Observation(const std::string& line)
    {
        if (line.length() < 79)
        {
            throw std::runtime_error("Invalid observation line length.");
        }

        obs_str = line;

        perm = line.substr(0, 5);
        prov = line.substr(5, 7);
        std::string date_str = line.substr(15, 16);
        std::string RA_str = line.substr(32, 11);
        std::string DEC_str = line.substr(44, 11);
        std::string mag_str = line.substr(65, 4);
        obs_code = line.substr(77, 3);

        int year, month;
        double day_frac;
        std::istringstream ds(date_str);
        ds >> year >> month >> day_frac;

        int day_int = static_cast<int>(day_frac);
        double decimal_day = day_frac - day_int;
        double secs = decimal_day * 86400.0;

        char utc_str[64];
        std::snprintf(utc_str, sizeof(utc_str),
            "%04d-%02d-%02dT00:00:00", year, month, day_int);

        SpiceDouble base_et;
        str2et_c(utc_str, &base_et);
        et = base_et + secs;

        std::istringstream ra_stream(RA_str);
        double ra_h, ra_m, ra_s;
        ra_stream >> ra_h >> ra_m >> ra_s;
        RA = ra_h * 15.0 + ra_m * 15.0 / 60.0 + ra_s * 15.0 / 3600.0;

        // Parse DEC
        int dec_sign = (line[44] == '-') ? -1 : 1;
        std::istringstream dec_stream(DEC_str);
        double dec_d, dec_m, dec_s;
        dec_stream >> dec_d >> dec_m >> dec_s;
        DEC = dec_sign * (std::abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0);

        // Parse magnitude
        try {
            mag = std::stod(mag_str);
        }
        catch (...) {
            mag = 0.0;
        }
    }

    std::string toString() const {
        return obs_str;
    }
};

// Return type: [body_index][position_or_velocity][x/y/z]
using StateMatrix = std::array<std::array<std::array<double, 3>, 2>, 9>; // I don't even know why I did it this way
// --- --- --- STRUCT AND CLASS DEFINITIONS --- --- --- [END]

// --- --- --- MISC UTILS --- --- --- [START]
void loadAllKernels(const std::string& directory) 
{
    std::string search_path = directory + "\\*.*";
    WIN32_FIND_DATAA fd;
    HANDLE hFind = ::FindFirstFileA(search_path.c_str(), &fd);

    if (hFind == INVALID_HANDLE_VALUE) {
        std::cerr << "Unable to open directory: " << directory << '\n';
        return;
    }

    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            std::string filepath = directory + "\\" + fd.cFileName;
            furnsh_c(filepath.c_str());
            // std::cout << "Loaded kernel: " << filepath << '\n';
        }
    } while (::FindNextFileA(hFind, &fd));

    ::FindClose(hFind);
}

std::string trim(const std::string& str) // trim trailing/leading whitespace, nothing fancy, stole it from somewhere
{
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

int digitsBeforeDecimal(double value) {
    if (value == 0.0) return 1;
    value = std::fabs(value); // Make it positive
    return static_cast<int>(std::floor(std::log10(value))) + 1;
}

double etToJD(SpiceDouble et)
{
    return 2451545.0 + et / 86400.0;
}

SpiceDouble JDToEt(double JD)
{
    return (JD - 2451545.0) * 86400.0;
}

std::string alignDecimal(const std::string& prefix, double value, int precision = 5) {
    std::ostringstream oss;

    // Width for just the number, not including the prefix
    const int number_field_width = 4 + 1 + precision; // 4 digits before decimal, 1 dot, rest after
    oss << prefix;

    // Format number with fixed precision and right-align in number_field_width
    std::ostringstream numstream;
    numstream << std::fixed << std::setprecision(precision) << value;
    std::string num_str = numstream.str();

    int pad_spaces = number_field_width - static_cast<int>(num_str.length());
    if (pad_spaces > 0) {
        oss << std::string(pad_spaces, ' ');
    }

    oss << num_str;
    return oss.str();
}

// this is the fun and simpler part of the equation
std::vector<double> stateVector2Kepler(Vec3 r_equ, Vec3 v_equ, SpiceDouble et, double mu = 1.3271244004193938E+11)
{
    // convert from equatorial to ecliptic frame
    double rot[3][3];
    pxform_c("J2000", "ECLIPJ2000", et, rot);

    double r_equ_arr[3] = { r_equ.x, r_equ.y, r_equ.z };
    double v_equ_arr[3] = { v_equ.x, v_equ.y, v_equ.z };

    double r_arr[3] = { 0, 0, 0 };
    double v_arr[3] = { 0, 0, 0 };

    mxv_c(rot, r_equ_arr, r_arr);
    mxv_c(rot, v_equ_arr, v_arr);

    Vec3 r = Vec3(r_arr[0], r_arr[1], r_arr[2]);
    Vec3 v = Vec3(v_arr[0], v_arr[1], v_arr[2]);

    // now calculate orbital params
    double r_mag = r.mag();
    double v_mag = v.mag();

    Vec3 h = r.cross(v);
    double h_mag = h.mag();

    double inclination = rad2deg(acos(h.z / h_mag));

    Vec3 k = Vec3(0, 0, 1);
    Vec3 n = k.cross(h);
    double n_mag = n.mag();

    double omega;
    double eccentricity;
    double arg_periapsis;
    double true_anomaly;
    double sma;
    double mean_anomaly;

    if (n_mag != 0)
    {
        omega = rad2deg(acos(n.x / n_mag));
        if (n.y < 0)
        {
            omega = 360 - omega;
        }
    }
    else
    {
        omega = 0;
    }

    Vec3 e_vec = (v.cross(h) - r * mu / r_mag) * (1 / mu);
    eccentricity = e_vec.mag();

    if (n_mag != 0)
    {
        if (eccentricity != 0)
        {
            arg_periapsis = rad2deg(acos(n.dot(e_vec) / (n_mag * eccentricity)));
            if (e_vec.z < 0)
            {
                arg_periapsis = 360 - arg_periapsis;
            }
        }
        else
        {
            arg_periapsis = 0;
        }
    }
    else
    {
        arg_periapsis = 0;
    }

    if (eccentricity != 0)
    {
        true_anomaly = rad2deg(acos(e_vec.dot(r) / (eccentricity * r_mag)));
        if (r.dot(v) < 0)
        {
            true_anomaly = 360 - true_anomaly;
        }
    }
    else
    {
        true_anomaly = rad2deg(acos(r.normalized().dot(v.normalized())));
    }

    double specific_energy = v_mag * v_mag / 2 - mu / r_mag;
    if (abs(eccentricity - 1) > 1e-8)
    {
        sma = -mu / (2 * specific_energy);
    }
    else
    {
        sma = 999999;
    }

    if (eccentricity < 1)
    {
        double E = 2 * atan(tan(deg2rad(true_anomaly) / 2) * sqrt((1 - eccentricity) / (1 + eccentricity)));
        if (E < 0)
        {
            E = E + 2 * pi_c();
        }

        mean_anomaly = rad2deg(E - eccentricity * sin(E));
    }
    else if(eccentricity > 1)
    {
        double F = 2 * atanh(tan(deg2rad(true_anomaly) / 2) * sqrt((eccentricity - 1) / (eccentricity + 1)));
        mean_anomaly = rad2deg(eccentricity * sinh(F) - F);
    }
    else
    {
        mean_anomaly = -1.0; // random val.
    }

    return std::vector<double> {sma, eccentricity, inclination, omega, arg_periapsis, true_anomaly, mean_anomaly};

}

std::vector<double> cartezian2spherical (Vec3 v)
{
    double d = v.mag();
    double RA = rad2deg(atan2(v.y, v.x));
    double DEC = rad2deg(asin(v.z / d));

    if (RA < 0)
    {
        RA = RA + 360;
    }

    return std::vector<double> {RA, DEC};
}

std::string RATohms(double ra_deg) {
    double total_hours = ra_deg / 15.0;
    int hours = static_cast<int>(total_hours);
    double remainder_minutes = (total_hours - hours) * 60.0;
    int minutes = static_cast<int>(remainder_minutes);
    double seconds = (remainder_minutes - minutes) * 60.0;

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setw(2) << minutes << ":"
        << std::fixed << std::setprecision(1) << std::setw(4) << seconds;
    return oss.str();
}

std::string DECToDms(double dec_deg) {
    char sign = (dec_deg >= 0) ? '+' : '-';
    double abs_deg = std::fabs(dec_deg);
    int degrees = static_cast<int>(abs_deg);
    double remainder_minutes = (abs_deg - degrees) * 60.0;
    int minutes = static_cast<int>(remainder_minutes);
    double seconds = (remainder_minutes - minutes) * 60.0;

    std::ostringstream oss;
    oss << sign << std::setfill('0') << std::setw(2) << degrees << ":"
        << std::setw(2) << minutes << ":"
        << std::fixed << std::setprecision(1) << std::setw(4) << seconds;
    return oss.str();
}

std::string doubleToResidualStr(double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << std::abs(value);
    // why is the sign printed at the end of the number on MPECs, wtf
    oss << (value >= 0 ? '+' : '-');
    return oss.str();
}

// I could've made a KeplerOrbit class or whatever but this works too
void printOrbitalElements(std::vector<double> orbel)
{
    std::cout << "a: " << orbel[0] / AU << " e: " << orbel[1] << " i: " << orbel[2] << " Node: " << orbel[3] << " Peri: " << orbel[4] << " M: " << orbel[6] << "\n";
}

std::vector<double> flattenResiduals(const std::pair<std::vector<double>, std::vector<double>>& residual_pair) {
    const std::vector<double>& ra_res = residual_pair.first;
    const std::vector<double>& dec_res = residual_pair.second;
    std::vector<double> flat;
    flat.reserve(ra_res.size() + dec_res.size());

    for (std::size_t i = 0; i < ra_res.size(); ++i) {
        flat.push_back(ra_res[i]);
        flat.push_back(dec_res[i]);
    }

    return flat;
}

// --- --- --- MISC UTILS --- --- --- [END]

// --- --- --- OBSERVATIONS AND OBSERVERS HANDLING --- --- --- [START]
std::vector<Observatory> readObsCodes(const std::string& filepath = "data/ObsCodes.dat")
{
    std::cout << "Reading ObsCodes... " << std::flush;
    std::vector<Observatory> obscodes;
    int skips = 0;

    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        std::cerr << "Could not open file: " << filepath << "\n";
        return obscodes;
    }

    std::string line;
    std::getline(infile, line);

    while (std::getline(infile, line))
    {
        try
        {
            std::string code_str = line.substr(0, 3);
            std::string l_str = line.substr(6, 7);
            std::string c_str = line.substr(13, 8);
            std::string s_str = line.substr(21, 9);
            std::string name_str = trim(line.substr(30));

            obscodes.push_back(Observatory(name_str, code_str, strtod(l_str.c_str(), NULL), strtod(c_str.c_str(), NULL), strtod(s_str.c_str(), NULL)));
        }
        catch (...)
        {
            ++skips;
            continue;
        }
    }

    std::cout << "Done, skipped " << skips << " obscodes.\n";
    return obscodes;
}

std::vector<Observation> readObsFile(const std::string& filename = "primary.obs") {
    std::cout << "Reading observations file: " << filename << " ... " << std::flush;
    std::vector<Observation> observations;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return observations;
    }

    std::string line;
    while (std::getline(infile, line)) {
        try {
            Observation obs(line);
            observations.push_back(obs);
        }
        catch (...) {
            // Could not parse line, skip
        }
    }

    std::cout << "Done.\n\n";
    return observations;
}

// lat - lon to an geocentric Earth-fixed vector (hence ECEF)
Vec3 geodeticToECEF(double lat_deg, double lon_deg)
{
    double a = 6378.137;
    double f = 1 / 298.257223563; // flattening. honestly not even needed.
    double e2 = f * (2 - f);

    double lat = deg2rad(lat_deg);
    double lon = deg2rad(lon_deg);

    double slat = sin(lat);
    double N = a / sqrt(1 - e2 * slat * slat);

    double clat = cos(lat);

    double x = N * clat * cos(lon);
    double y = N * clat * sin(lon);
    double z = N * (1 - e2) * slat;

    return Vec3(x, y, z);
}

// give obscode and time, get heliocentric position (equatorial J2000)
Vec3 getObserverPos(Observatory obsv, SpiceDouble et)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

    double lat_rad = atan2(obsv.s, obsv.c);
    double lat_deg = lat_rad * 180 / pi_c();

    Vec3 ecef = geodeticToECEF(lat_deg, obsv.l);
    SpiceDouble rot[3][3];
    pxform_c("ITRF93", "J2000", et, rot);

    SpiceDouble ecef_sd[3] = { ecef.x, ecef.y, ecef.z };
    SpiceDouble observer_j2000[3];
    mxv_c(rot, ecef_sd, observer_j2000);

    return earth_pos + Vec3(ecef_sd[0], ecef_sd[1], ecef_sd[2]);
}

// give obscode and time, get heliocentric velocity (equatorial J2000)
Vec3 getObserverVel(Observatory obsv, SpiceDouble et)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_vel = Vec3(state[3], state[4], state[5]);

    double lat_rad = atan2(obsv.s, obsv.c);
    double lat_deg = lat_rad * 180 / pi_c();

    Vec3 ecef = geodeticToECEF(lat_deg, obsv.l);
    SpiceDouble xform[6][6];
    sxform_c("ITRF93", "J2000", et, xform);

    SpiceDouble ecef_mxvform[6] = { ecef.x, ecef.y, ecef.z, 0, 0, 0 };
    SpiceDouble mxvg_out[6];
    mxvg_c(xform, ecef_mxvform, 6, 6, mxvg_out);

    Vec3 obsv_geo_vel = Vec3(mxvg_out[3], mxvg_out[4], mxvg_out[5]);

    Vec3 obsv_vel = earth_vel + obsv_geo_vel;
    return obsv_vel;
}

// give time, observatory and heliocentric object position, get RA and DEC
std::vector<double> getRADEC(SpiceDouble et, Observatory obsv, Vec3 p)
{
    Vec3 obsv_pos = getObserverPos(obsv, et);
    Vec3 rel_pos = p - obsv_pos;

    return cartezian2spherical(rel_pos);
}

std::vector<double> getGeocentricRADEC(SpiceDouble et, Vec3 p)
{
    SpiceDouble state[6], lt;
    spkezr_c("EARTH", et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
    Vec3 earth_pos = Vec3(state[0], state[1], state[2]);

    Vec3 rel_pos = p - earth_pos;

    return cartezian2spherical(rel_pos);
}
// --- --- --- OBSERVATIONS AND OBSERVERS HANDLING --- --- --- [END]

// --- --- --- ORBIT PROPAGATION --- --- --- [START]
// I do regret doing this the way I did it, but it doesn't bite, it can stay
StateMatrix getSolarSystemStates(SpiceDouble et)
{
    const char* bodies[9] = {
        "SUN",                 // index 0
        "MERCURY BARYCENTER",  // index 1
        "VENUS BARYCENTER",    // index 2
        "EARTH BARYCENTER",    // index 3
        "MARS BARYCENTER",     // index 4
        "JUPITER BARYCENTER",  // index 5
        "SATURN BARYCENTER",   // index 6
        "URANUS BARYCENTER",   // index 7
        "NEPTUNE BARYCENTER"   // index 8
    };

    StateMatrix states{};

    for (int i = 0; i < 9; ++i) 
    {
        SpiceDouble state[6];
        SpiceDouble lt;

        spkezr_c(bodies[i], et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);

        for (int j = 0; j < 3; ++j) 
        {
            states[i][0][j] = state[j];     // Position (km)
            states[i][1][j] = state[j + 3]; // Velocity (km/s)
        }
    }

    return states;
}

// not running numbers on Schwardschild metrics
Vec3 gravAccel(MinorPlanet& mp, std::vector<Body>& bodies)
{
    Vec3 accel{};

    for (Body& body : bodies)
    {
        double dist = (body.pos - mp.pos).mag();
        Vec3 grav_dir = (body.pos - mp.pos) / dist;
        double grav_mag = body.GM / (dist * dist);

        accel = accel + grav_dir * grav_mag;
    }

    return accel;
}

// symplectic 8th order orbit integrator
// RKN or whatever also works probably
// I like this because I already use this a lot in other projects too, 
// so I don't have to spend brain power again implementing something else.
void stepYoshida8(MinorPlanet& mp, std::vector<Body>& bodies, SpiceDouble date0_et, double dt)
{
    constexpr double w1 = 0.311790812418427e0;
    constexpr double w2 = -0.155946803821447e1;
    constexpr double w3 = -0.167896928259640e1;
    constexpr double w4 = 0.166335809963315e1;
    constexpr double w5 = -0.106458714789183e1;
    constexpr double w6 = 0.136934946416871e1;
    constexpr double w7 = 0.629030650210433e0;
    constexpr double w0 = 1.65899088454396;

    constexpr double ds[15] = { w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7 };
    constexpr double cs[16] = { 0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165 };

    SpiceDouble et = date0_et;

    for (int i = 0; i < 15; ++i) 
    {
        mp.pos += mp.vel * (cs[i] * dt);

        // Advance SPICE time
        et += cs[i] * dt;

        // Update body positions
        for (Body& body : bodies) 
        {
            SpiceDouble state[6], lt;
            spkezr_c(body.name.c_str(), et, "J2000", "NONE", "SOLAR SYSTEM BARYCENTER", state, &lt);
            body.pos = { state[0], state[1], state[2] };
        }

        Vec3 accel = gravAccel(mp, bodies);
        mp.vel += accel * (ds[i] * dt);
    }

    mp.pos += mp.vel * (cs[15] * dt);
}

// orbit propagator main function
std::vector<Vec3> propagate(Vec3 p0, Vec3 v0, SpiceDouble date_init, SpiceDouble date_final, double dt)
{
    MinorPlanet mp = MinorPlanet(p0, v0);

    const char* body_names[9] = {
        "SUN",
        "MERCURY BARYCENTER",
        "VENUS BARYCENTER",
        "EARTH BARYCENTER",
        "MARS BARYCENTER",
        "JUPITER BARYCENTER",
        "SATURN BARYCENTER",
        "URANUS BARYCENTER",
        "NEPTUNE BARYCENTER"
    };

    const double body_GMs[9] = {
        1.3271244004193938E+11,
        2.2031780000000021E+04,
        3.2485859200000006E+05,
        4.0350323550225981E+05,
        4.2828375214000022E+04,
        1.2671276480000021E+08,
        3.7940585200000003E+07,
        5.7945486000000080E+06,
        6.8365271005800236E+06
    };

    std::vector<Body> bodies;
    StateMatrix system_state = getSolarSystemStates(date_init);

    for (int i = 0; i < 8; i++)
    {
        bodies.push_back(Body(body_names[i], body_GMs[i]));
        bodies[i].pos = Vec3(system_state[i][0]);
    }

    double time_interval = date_final - date_init;

    if (dt <= 0)
    {
        dt = time_interval / 4;
        while (dt > 10 * 86400)
        {
            dt = dt / 2;
        }
    }

    int N_cycles = (int)(time_interval / dt + 1);
    SpiceDouble date_final_actual = date_init + N_cycles * dt;

    for (int cycle = 0; cycle < N_cycles; cycle++)
    {
        SpiceDouble cycle_date = date_init + cycle * dt;
        stepYoshida8(mp, bodies, cycle_date, dt);
    }

    // final little correction
    double extra_time = N_cycles * dt - time_interval;
    mp.pos = mp.pos - mp.vel * extra_time;

    return std::vector<Vec3> {Vec3(mp.pos), Vec3(mp.vel)};
}

// very little difference with propagate(), only in handling of dt (time step) value
std::vector<Vec3> backpropagate(Vec3 p0, Vec3 v0, SpiceDouble date_init, SpiceDouble date_final, double dt)
{
    MinorPlanet mp = MinorPlanet(p0, v0);

    const char* body_names[9] = {
        "SUN",
        "MERCURY BARYCENTER",
        "VENUS BARYCENTER",
        "EARTH BARYCENTER",
        "MARS BARYCENTER",
        "JUPITER BARYCENTER",
        "SATURN BARYCENTER",
        "URANUS BARYCENTER",
        "NEPTUNE BARYCENTER"
    };

    const double body_GMs[9] = {
        1.3271244004193938E+11,
        2.2031780000000021E+04,
        3.2485859200000006E+05,
        4.0350323550225981E+05,
        4.2828375214000022E+04,
        1.2671276480000021E+08,
        3.7940585200000003E+07,
        5.7945486000000080E+06,
        6.8365271005800236E+06
    };

    std::vector<Body> bodies;
    StateMatrix system_state = getSolarSystemStates(date_init);

    for (int i = 0; i < 8; i++)
    {
        bodies.push_back(Body(body_names[i], body_GMs[i]));
        bodies[i].pos = Vec3(system_state[i][0]);
    }

    double time_interval = date_final - date_init;

    if (dt >= 0)
    {
        dt = time_interval / 4;
        while (dt < -10 * 86400)
        {
            dt = dt / 2;
        }
    }

    int N_cycles = (int)(time_interval / dt + 1);
    SpiceDouble date_final_actual = date_init + N_cycles * dt;

    for (int cycle = 0; cycle < N_cycles; cycle++)
    {
        SpiceDouble cycle_date = date_init + cycle * dt;
        stepYoshida8(mp, bodies, cycle_date, dt);
    }

    // final little correction
    double extra_time = N_cycles * dt - time_interval;
    mp.pos = mp.pos - mp.vel * extra_time;

    return std::vector<Vec3> {Vec3(mp.pos), Vec3(mp.vel)};
}
// --- --- --- ORBIT PROPAGATION --- --- --- [END]

// this guesses an initial velocity assuming observations 1 and 2 are the same distance from the Solar System Barycenter
Vec3 guessCircularV0(Vec3 p0, std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map)
{
    double r = p0.mag();
    Observation o1 = obs_all[1];

    // We know the heliocentric distance because we already have p0 located. We also know the 
    // direction of the object from the observer. We also know the position of the observer.
    // Using these, locate the object in heliocentric space.

    // Fixing the distance from SSB gives us a spherical surface.
    // Fixing the direction gives us a semi-infinite line segment.
    // Their intersection should be the object location.
    // 
    // See https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    //
    // Now we find d.

    Vec3 u = Vec3(o1.RA, o1.DEC);
    Vec3 o = getObserverPos(obscode_map[o1.obs_code], o1.et);

    double u_dot_oc = u.dot(o);
    double o_mag = o.mag();

    double disc = u_dot_oc * u_dot_oc - (o_mag * o_mag - r * r);
    double d_1 = -u_dot_oc + sqrt(disc);
    double d_2 = -u_dot_oc - sqrt(disc);

    double d = d_1;

    if (d_1 < 0)
    {
        d = d_2;
    }

    Vec3 p1 = o + u * d; // assumed position at second observation

    // now we find velocity
    // Bill Gray (author of find_orb and astcheck) has a good explanation for what I'm doing so I won't bother repeating here
    // https://www.projectpluto.com/herget.htm
    Vec3 p_mid = (p0 + p1) * 0.5;

    // get an acceleration guess at midpoint
    double mu = 1.3271244004193938E+11; // Sun standrad gravitational parameter
    Vec3 a_mid = p_mid.normalized() * (- mu) / (r * r);

    // time it takes to go from observation 1 to 2
    double t_diff = o1.et - obs_all[0].et;

    Vec3 v0 = (p1 - p0) / t_diff - a_mid * t_diff / 2;
    return v0;
}

std::pair<std::vector<double>, std::vector<double>> getResiduals(Vec3 p0, Vec3 v0, std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map)
{
    Observation o1 = obs_all[0];

    Vec3 p0_temp = p0;
    Vec3 v0_temp = v0;
    SpiceDouble et_temp = o1.et;

    std::vector<double> RA_res;
    std::vector<double> DEC_res;
    for (int idx_o = 0; idx_o < obs_all.size(); idx_o++)
    {
        Observation o = obs_all[idx_o];

        std::vector<Vec3> prop_res = propagate(p0_temp, v0_temp, et_temp, o.et, -1);
        Vec3 p_check = prop_res[0];
        Vec3 v_new = prop_res[1];
        std::vector<double> RADEC_prop = getRADEC(o.et, obscode_map[o.obs_code], p_check);

        double RA_prop = RADEC_prop[0];
        double DEC_prop = RADEC_prop[1];

        double RA_diff = o.RA - RA_prop;
        if (RA_diff > 180.0) RA_diff -= 360.0;
        else if (RA_diff < -180.0) RA_diff += 360.0;

        double RA_err = RA_diff * 3600;
        double DEC_err = (o.DEC - DEC_prop) * 3600;

        RA_res.push_back(RA_err);
        DEC_res.push_back(DEC_err);

        p0_temp = p_check;
        v0_temp = v_new;
        et_temp = o.et;
    }

    return std::pair<std::vector<double>, std::vector<double>> {RA_res, DEC_res};
}

// gets observations and state vectors, returns RA-DEC error score
double getErrorRADEC(Vec3 p0, Vec3 v0, std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map)
{
    std::pair<std::vector<double>, std::vector<double>> residuals = getResiduals(p0, v0, obs_all, obscode_map);
    std::vector<double> RA_res = residuals.first;
    std::vector<double> DEC_res = residuals.second;

    double err_val = 0;
    for (int iter = 0; iter < RA_res.size(); iter++)
    {
        err_val += RA_res[iter] * RA_res[iter] + DEC_res[iter] * DEC_res[iter];
    }

    return err_val;
}

std::vector<Vec3> initialGuess(std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map, double R1)
{
    Observation o1 = obs_all[0];
    double best_R = 2.0 * AU; // default-ish value
    if (R1 < 0)
    {
        // test for best fitting observer - object distance
        double R1_ts[63] = { 0.1 * AU, 0.2 * AU, 0.3 * AU, 0.5 * AU, 0.75 * AU, 1 * AU,
                            1.25 * AU, 1.5 * AU, 1.75 * AU, 2 * AU, 2.25 * AU, 2.5 * AU, 
                            2.75 * AU, 3 * AU, 3.25 * AU, 3.5 * AU, 3.75 * AU, 4 * AU,
                            5 * AU, 6 * AU, 7 * AU, 8 * AU, 9 * AU, 10 * AU, 11 * AU, 12 * AU,
                            13 * AU, 14 * AU, 15 * AU, 16 * AU, 18 * AU, 20 * AU, 22 * AU, 24 * AU,
                            26 * AU, 28 * AU, 29 * AU, 30 * AU, 31 * AU, 32 * AU, 33 * AU, 34 * AU,
                            35 * AU, 36 * AU, 37 * AU, 38 * AU, 39 * AU, 40 * AU, 41 * AU, 42 * AU,
                            43 * AU, 44 * AU, 45 * AU, 46 * AU, 47 * AU, 48 * AU, 49 * AU, 50 * AU,
                            51 * AU, 52 * AU, 53 * AU, 54 * AU, 55 * AU};
        double R_err_prev = 1E+15; // arbitrary large number

        for (int idx_Rt = 0; idx_Rt < 63; idx_Rt++)
        {
            double R1_t = R1_ts[idx_Rt];
            Vec3 p0 = getObserverPos(obscode_map[o1.obs_code], o1.et) + Vec3(o1.RA, o1.DEC) * R1_t;
            Vec3 v0 = guessCircularV0(p0, obs_all, obscode_map);

            // check current error
            double err_val = getErrorRADEC(p0, v0, obs_all, obscode_map);

            if (err_val < R_err_prev)
            {
                best_R = R1_t;
                R_err_prev = err_val;
            }
        }
    }
    else
    {
        best_R = R1;
    }

    Vec3 p0 = getObserverPos(obscode_map[o1.obs_code], o1.et) + Vec3(o1.RA, o1.DEC) * best_R;
    Vec3 v0 = guessCircularV0(p0, obs_all, obscode_map);

    return std::vector<Vec3> {p0, v0};
}

std::vector<double> computeGradient(Vec3 p0, Vec3 v0, std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map)
{
    double err0 = getErrorRADEC(p0, v0, obs_all, obscode_map);
    std::vector<double> grad = {0, 0, 0, 0, 0, 0};

    double eps_vel = 1;
    double eps_pos = 1;

    grad[0] = (getErrorRADEC(p0, v0 + Vec3(1, 0, 0) * eps_vel, obs_all, obscode_map) - err0) / eps_vel;
    grad[1] = (getErrorRADEC(p0, v0 + Vec3(0, 1, 0) * eps_vel, obs_all, obscode_map) - err0) / eps_vel;
    grad[2] = (getErrorRADEC(p0, v0 + Vec3(0, 0, 1) * eps_vel, obs_all, obscode_map) - err0) / eps_vel;

    grad[3] = (getErrorRADEC(p0 + Vec3(1, 0, 0) * eps_pos, v0, obs_all, obscode_map) - err0) / eps_pos;
    grad[4] = (getErrorRADEC(p0 + Vec3(0, 1, 0) * eps_pos, v0, obs_all, obscode_map) - err0) / eps_pos;
    grad[5] = (getErrorRADEC(p0 + Vec3(0, 0, 1) * eps_pos, v0, obs_all, obscode_map) - err0) / eps_pos;

    return grad;
}

std::vector<std::array<double, 6>> computeJacobian(Vec3 p0, Vec3 v0, std::vector<Observation> obs_all,
    std::unordered_map<std::string, Observatory>& obscode_map,
    double eps_v = 1e-2, double eps_p = 1e3) 
{
    int N = obs_all.size();
    std::vector<std::array<double, 6>> J(2 * N);  // 2 residuals per observation

    std::pair<std::vector<double>, std::vector<double>> base_pair = getResiduals(p0, v0, obs_all, obscode_map);
    std::vector<double> base_residuals = flattenResiduals(base_pair);

    for (int j = 0; j < 3; ++j) 
    {
        // perturb velocity
        Vec3 dv(0, 0, 0);
        if (j == 0) dv.x = eps_v;
        else if (j == 1) dv.y = eps_v;
        else dv.z = eps_v;

        std::pair<std::vector<double>, std::vector<double>> plus_pair_v = getResiduals(p0, v0 + dv, obs_all, obscode_map);
        std::pair<std::vector<double>, std::vector<double>> minus_pair_v = getResiduals(p0, v0 - dv, obs_all, obscode_map);

        std::vector<double> res_plus_v = flattenResiduals(plus_pair_v);
        std::vector<double> res_minus_v = flattenResiduals(minus_pair_v);

        for (int i = 0; i < 2 * N; ++i) 
        {
            J[i][j] = (res_plus_v[i] - res_minus_v[i]) / (2.0 * eps_v);
        }

        // perturb position
        Vec3 dp(0, 0, 0);
        if (j == 0) dp.x = eps_p;
        else if (j == 1) dp.y = eps_p;
        else dp.z = eps_p;

        std::pair<std::vector<double>, std::vector<double>> plus_pair_p = getResiduals(p0 + dp, v0, obs_all, obscode_map);
        std::pair<std::vector<double>, std::vector<double>> minus_pair_p = getResiduals(p0 - dp, v0, obs_all, obscode_map);

        std::vector<double> res_plus_p = flattenResiduals(plus_pair_p);
        std::vector<double> res_minus_p = flattenResiduals(minus_pair_p);

        for (int i = 0; i < 2 * N; ++i) 
        {
            J[i][j + 3] = (res_plus_p[i] - res_minus_p[i]) / (2.0 * eps_p);
        }
    }

    return J;
}

std::array<double, 6> computeGaussNewtonStep(std::vector<std::array<double, 6>>& J, std::vector<double>& r) 
{
    int M = J.size();

    double JTJ[6][6] = {};
    std::array<double, 6> JTr = {};

    for (int i = 0; i < M; ++i) 
    {
        for (int j = 0; j < 6; ++j) 
        {
            JTr[j] += J[i][j] * r[i];
            for (int k = 0; k < 6; ++k) 
            {
                JTJ[j][k] += J[i][j] * J[i][k];
            }
        }
    }

    // JTJ * dx = -JTr
    double A[6][7]; // augmented matrix

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            A[i][j] = JTJ[i][j];
        }
        A[i][6] = -JTr[i];
    }

    // Gaussian elimination
    for (int i = 0; i < 6; ++i) 
    {
        int pivot = i;
        for (int j = i + 1; j < 6; ++j) 
        {
            if (std::fabs(A[j][i]) > std::fabs(A[pivot][i])) 
            {
                pivot = j;
            }
        }
        for (int j = 0; j < 7; ++j) 
        {
            std::swap(A[i][j], A[pivot][j]);
        }

        double div = A[i][i];
        for (int j = i; j < 7; ++j) 
        {
            A[i][j] /= div;
        }

        for (int j = 0; j < 6; ++j) 
        {
            if (j != i) 
            {
                double f = A[j][i];
                for (int k = i; k < 7; ++k) 
                {
                    A[j][k] -= f * A[i][k];
                }
            }
        }
    }

    std::array<double, 6> dx;
    for (int i = 0; i < 6; ++i) 
    {
        dx[i] = A[i][6];
    }

    return dx;
}

// the actual orbit determination algorithm
std::vector<Vec3> determineOrbit(std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map,
    double R1, int maxiter)
{
    std::cout << "Orbit determination...\n";

    Observation o1 = obs_all[0];
    Observation o2 = obs_all[(int)(obs_all.size() / 2)];
    Observation o3 = obs_all[obs_all.size() - 1];

    SpiceDouble date_init = o1.et;
    SpiceDouble date_final = o2.et;
    SpiceDouble date_check = o3.et;

    std::vector<Vec3> initialGuessSv = initialGuess(obs_all, obscode_map, R1);
    Vec3 p0 = initialGuessSv[0];
    Vec3 v0 = initialGuessSv[1];

    double err_val = getErrorRADEC(p0, v0, obs_all, obscode_map);

    std::vector<double> orbital_elems_init = stateVector2Kepler(p0, v0, o1.et);
    std::cout << "Initial guess:\n";
    printOrbitalElements(orbital_elems_init);
    std::cout << "Initial pos. [km]  : "; p0.printout();
    std::cout << "Initial vel. [km/s]: "; v0.printout();
    std::cout << "\nIterating to refine solution";

    double adjust_factor = 1;

    bool good_fit = false;
    int retry_count = 0;
    int MAX_RETRY = maxiter;
    while (!good_fit && retry_count <= MAX_RETRY)
    {
        // check current error
        err_val = getErrorRADEC(p0, v0, obs_all, obscode_map);

        if (err_val < 10)
        {
            good_fit = true;
            break;
        }

        // get gradient in state-vector-space
        std::vector<double> grad = computeGradient(p0, v0, obs_all, obscode_map);

        std::pair<std::vector<double>, std::vector<double>> residuals = getResiduals(p0, v0, obs_all, obscode_map);
        std::vector<double> r = flattenResiduals(residuals);
        std::vector<std::array<double, 6>> J = computeJacobian(p0, v0, obs_all, obscode_map);
        std::array<double, 6> dx = computeGaussNewtonStep(J, r);

        Vec3 v0_new = v0 + Vec3(dx[0], dx[1], dx[2]) * adjust_factor;
        Vec3 p0_new = p0 + Vec3(dx[3], dx[4], dx[5]) * adjust_factor;

        double err_val_new = getErrorRADEC(p0_new, v0_new, obs_all, obscode_map);

        if (err_val_new < err_val)
        {
            v0 = v0_new;
            p0 = p0_new;
            adjust_factor *= 1.8;
        }
        else
        {
            adjust_factor *= 0.1;
        }

        retry_count += 1;
        std::cout << ".";
    }

    std::cout << "\nDone iterating.\n";
    return std::vector<Vec3> {p0, v0};
}

std::pair<std::vector<Vec3>, double> propagateToNextMidnight(Vec3 p0, Vec3 v0, SpiceDouble t0)
{
    double JD_0 = etToJD(t0);
    double JD_1 = std::floor(JD_0) + 0.5 + std::ceil(JD_0 - std::floor(JD_0 + 0.5)); // go to next mindight
    SpiceDouble t1 = JDToEt(JD_1);

    return { propagate(p0, v0, t0, t1, -1), JD_1 };
}

std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> generateEphemeris(SpiceDouble t0, Vec3 p0, Vec3 v0, std::vector<SpiceDouble> ts)
{
    std::vector<double> RA_result;
    std::vector<double> DEC_result;

    for (int idx_tdiff = 0; idx_tdiff < ts.size(); idx_tdiff++)
    {
        if (ts[idx_tdiff] > t0)
        {
            std::vector<Vec3> sv_prop = propagate(p0, v0, t0, ts[idx_tdiff], -1);
            Vec3 p_t = sv_prop[0];
            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p_t);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);
        }
        else if (ts[idx_tdiff] == t0) // the date is epoch date
        {
            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p0);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);
        }
        else
        {
            std::vector<Vec3> sv_prop = backpropagate(p0, v0, t0, ts[idx_tdiff], 1);
            Vec3 p_t = sv_prop[0];
            std::vector<double> RADEC = getGeocentricRADEC(ts[idx_tdiff], p_t);
            RA_result.push_back(RADEC[0]);
            DEC_result.push_back(RADEC[1]);
        }
    }

    // convert RA and DEC from decimal degrees to hh:mm:ss and sDD:mm:ss
    std::vector<std::string> RA_hms;
    std::vector<std::string> DEC_Dms;
    std::vector<std::string> dates;

    for (int idx_result = 0; idx_result < RA_result.size(); idx_result++)
    {
        SpiceChar utcstr[20];
        et2utc_c(ts[idx_result], "ISOC", 0, 20, utcstr);
        dates.push_back(utcstr);
        RA_hms.push_back(RATohms(RA_result[idx_result]));
        DEC_Dms.push_back(DECToDms(DEC_result[idx_result]));
    }

    return std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> {dates, RA_hms, DEC_Dms};
}

std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> generateEphemerisMPEC(Vec3 p0, Vec3 v0, SpiceDouble epoch_et)
{
    // get today's date with clock at midnight
    std::time_t now = std::time(nullptr);
    std::tm utc_tm{};
    errno_t err = gmtime_s(&utc_tm, &now);
    if (err) {
        throw std::runtime_error("gmtime_s failed");
    }

    utc_tm.tm_hour = 0;
    utc_tm.tm_min = 0;
    utc_tm.tm_sec = 0;

    std::time_t utc_midnight = _mkgmtime(&utc_tm);

    char utc_str[32];
    std::strftime(utc_str, sizeof(utc_str), "%Y-%m-%dT%H:%M:%S", &utc_tm);

    SpiceDouble et_today;
    utc2et_c(utc_str, &et_today);

    std::vector<SpiceDouble> target_ets;
    target_ets.push_back(et_today - 30 * 86400.0);
    target_ets.push_back(et_today - 15 * 86400.0);
    target_ets.push_back(et_today - 7 * 86400.0);
    target_ets.push_back(et_today - 1 * 86400.0);
    target_ets.push_back(et_today);
    target_ets.push_back(et_today + 1 * 86400.0);
    target_ets.push_back(et_today + 7 * 86400.0);
    target_ets.push_back(et_today + 15 * 86400.0);
    target_ets.push_back(et_today + 30 * 86400.0);

    return generateEphemeris(epoch_et, p0, v0, target_ets);
}

// styles the output in the form of a Minor Planet Center MPEC for familiarity (and therefore easy reading)
void printQuasiMPEC(std::vector<double> orbital_elements,
    std::vector<Observation> observations,
    std::unordered_map<std::string, Observatory>& observatories,
    double epoch_JD, double orbital_period, Vec3 p0, Vec3 v0,
    std::ostream& out = std::cout)
{
    out << "This quasi-MPEC is NOT an official Minor Planet Electronic Circular!\n";
    out << "It is merely styled as such for familiarity to minor planet astronomers.\n\n";

    // print 80-column observations, also get observatory keys
    std::vector<std::string> obscodes = {};
    out << "Observations:\n";
    for (int idx_obs = 0; idx_obs < (int)observations.size(); idx_obs++)
    {
        Observation obs = observations[idx_obs];
        out << obs.obs_str << "\n";

        if (std::find(obscodes.begin(), obscodes.end(), obs.obs_code) == obscodes.end()) // if obscode NOT present in vector
        {
            obscodes.push_back(obs.obs_code);
        }
    }

    out << "\n";

    // print observatories details
    out << "Observatories:\n";
    for (int idx_obscode = 0; idx_obscode < (int)obscodes.size(); idx_obscode++)
    {
        out << obscodes[idx_obscode] << " " << observatories[obscodes[idx_obscode]].name << "\n";
    }

    out << "\n";

    // print orbital elements
    std::string objname;
    out << "Orbital elements:\n";
    if (observations[0].prov.empty())
    {
        objname = observations[0].perm;
    }
    else
    {
        objname = observations[0].prov;
    }
    out << objname << "\n";

    SpiceChar utcstr[20];
    int utclen = 20;
    et2utc_c(JDToEt(epoch_JD), "ISOC", 0, utclen, utcstr); // convert epoch JD to ISO

    double n = 360.0 / orbital_period * 86400.0; // deg/day
    double q = (1 - orbital_elements[1]) * orbital_elements[0] / AU; // perihelion dist.

    std::pair<std::vector<double>, std::vector<double>> residuals = getResiduals(p0, v0, observations, observatories);
    std::vector<double> RA_res = residuals.first;
    std::vector<double> DEC_res = residuals.second;

    std::vector<std::string> residual_elements;
    for (int idx_res = 0; idx_res < (int)RA_res.size(); idx_res++)
    {
        std::string residual_elem = observations[idx_res].obs_str.substr(17, 2) + observations[idx_res].obs_str.substr(20, 2) + observations[idx_res].obs_str.substr(23, 2) +
            " " + observations[idx_res].obs_code + " " + doubleToResidualStr(RA_res[idx_res]) + " " + doubleToResidualStr(DEC_res[idx_res]);
        residual_elements.push_back(residual_elem);
    }

    out << "Epoch: " << utcstr << " = JDT " << std::fixed << std::setprecision(1) << epoch_JD << "                   MPFT\n";
    out << alignDecimal("M ", orbital_elements[6]) << "              (2000.0)\n";
    out << alignDecimal("n ", n, 5) << "        " << alignDecimal("Peri. ", orbital_elements[4]) << "\n";
    out << alignDecimal("a ", orbital_elements[0] / AU, 5) << "        " << alignDecimal("Node  ", orbital_elements[3]) << "\n";
    out << alignDecimal("e ", orbital_elements[1], 5) << "        " << alignDecimal("Incl. ", orbital_elements[2]) << "\n";
    out << alignDecimal("P ", orbital_period / 86400.0 / 365, 2) << "\n";

    out << "Residuals in seconds of arc\n";
    const int cols = 3;
    const int col_width = 25;
    const int rows = (int)((residual_elements.size() + cols - 1) / cols);

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int idx = r + rows * c;
            if (idx < (int)residual_elements.size())
                out << std::left << std::setw(col_width) << residual_elements[idx];
            else
                out << std::left << std::setw(col_width) << "";
        }
        out << '\n';
    }

    // print out ephemeris
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> ephemeris = generateEphemerisMPEC(p0, v0, JDToEt(epoch_JD));

    out << "\nEphemeris:\n";
    out << objname << "                a,e,i = " << std::fixed << std::setprecision(2) << orbital_elements[0] / AU << ", "
        << std::fixed << std::setprecision(2) << orbital_elements[1] << ", "
        << std::fixed << std::setprecision(0) << orbital_elements[2] << "                   " <<
        "q = " << std::fixed << std::setprecision(4) << q << "\n";

    out << "Date                     R. A. (2000) Decl.\n";
    out << std::get<0>(ephemeris)[0] << "    " << std::get<1>(ephemeris)[0] << " " << std::get<2>(ephemeris)[0] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[1] << "    " << std::get<1>(ephemeris)[1] << " " << std::get<2>(ephemeris)[1] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[2] << "    " << std::get<1>(ephemeris)[2] << " " << std::get<2>(ephemeris)[2] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[3] << "    " << std::get<1>(ephemeris)[3] << " " << std::get<2>(ephemeris)[3] << "\n";
    out << std::get<0>(ephemeris)[4] << "    " << std::get<1>(ephemeris)[4] << " " << std::get<2>(ephemeris)[4] << "\n";
    out << std::get<0>(ephemeris)[5] << "    " << std::get<1>(ephemeris)[5] << " " << std::get<2>(ephemeris)[5] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[6] << "    " << std::get<1>(ephemeris)[6] << " " << std::get<2>(ephemeris)[6] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[7] << "    " << std::get<1>(ephemeris)[7] << " " << std::get<2>(ephemeris)[7] << "\n";
    out << "...\n";
    out << std::get<0>(ephemeris)[8] << "    " << std::get<1>(ephemeris)[8] << " " << std::get<2>(ephemeris)[8] << "\n";
    out << "\n";
    out << "M. P. F. T. Software         (C/) Copyleft MPFT Authors              Quasi-MPEC\n";
}


void printHelpMsg()
{
    std::cout << " === MPFT HELP ===\n";
    std::cout << "MPFT is an orbit determination software for minor planets in heliocentric orbit.\n\n";
    std::cout << "It requires SPICE kernels, an observation file in Minor Planet Center's 80-column format, ";
    std::cout << "and an ObsCodes file which holds the observatory codes used by the MPC.\n\n";
    std::cout << "The minimal invocation is merely 'mpft' - this assumes default locations for all files, which are:\n\n";
    std::cout << "SPICE Files:\n";
    std::cout << "data/SPICE/naif0012.tls\n";
    std::cout << "data/SPICE/de440.bsp\n";
    std::cout << "data/SPICE/pck00011.tpc\n";
    std::cout << "data/SPICE/earth_000101_250316_241218.bpc\n\n";
    std::cout << "Minor Planet Center Observatory Codes:\n";
    std::cout << "data/ObsCodes.dat\n\n";
    std::cout << "Your input observations (astrometry) in Obs80 format (MPC's 80-column format):\n";
    std::cout << "primary.obs\n\n";
    std::cout << "Optional arguments are as follows (can be entered in any order):\n\n";
    std::cout << "-obs <Obs80_input_file> -obscode <obscode_file> -spice <spice_folder_path> -R <initial_object_dist_to_observer_guess (AU)> -maxiter <max_orbit_refinement_iterations> -out <output_file_path>\n\n";
    std::cout << "If an initial distance guess is not given, the program tries several guesses and tries to pick the best one.\n";
    std::cout << "Default number of maximum iterations is 20. Output file will be named quasi-MPEC.txt by default.\n\n";
    std::cout << "MPFT was developed by H. A. Guler.\n\n";
    std::cout << "MPFT is licensed under GNU General Public License version 2.0 (GPL-2.0 License)\n\n";
}

int main(int argc, char *argv[])
{
    std::cout << "MPFT v0.5.0\n\n";

    // default parameters
    std::string obs_path = "primary.obs";
    int param_R = -1;
    std::string obscode_path = "data/ObsCodes.dat";
    std::string spice_path = "data/SPICE/";
    int param_maxiter = 20;
    std::string out_path = "quasi-MPEC.txt";

    // handle command line arguments
    // there is a more compact version of doing this but this is easier for my brain
    int argtype = 0;
    for (int idx_cmd = 1; idx_cmd < argc; idx_cmd++)
    {
        if (!strcmp(argv[idx_cmd], "-obs") || !strcmp(argv[idx_cmd], "-observations")) // also the default argument
        {
            argtype = 0;
        }
        else if (!strcmp(argv[idx_cmd], "-R"))
        {
            argtype = 1;
        }
        else if (!strcmp(argv[idx_cmd], "-obscode"))
        {
            argtype = 2;
        }
        else if (!strcmp(argv[idx_cmd], "-spice"))
        {
            argtype = 3;
        }
        else if (!strcmp(argv[idx_cmd], "-maxiter"))
        {
            argtype = 4;
        }
        else if (!strcmp(argv[idx_cmd], "-out") || !strcmp(argv[idx_cmd], "-output"))
        {
            argtype = 5;
        }
        else if (!strcmp(argv[idx_cmd], "-h") || !strcmp(argv[idx_cmd], "--help")) // two dashes because people are more used to it
        {
            printHelpMsg();
            return 0;
        }
        else 
        {
            switch (argtype)
            {
            case 0:
                obs_path = argv[idx_cmd];
                break;
            case 1:
                param_R = strtod(argv[idx_cmd], NULL) * AU;
                break;
            case 2:
                obscode_path = argv[idx_cmd];
                break;
            case 3:
                spice_path = argv[idx_cmd];
                break;
            case 4:
                param_maxiter = atoi(argv[idx_cmd]);
                break;
            case 5:
                out_path = argv[idx_cmd];
                break;
            }
        }
    }

    std::cout << "Loading SPICE kernels... ";
    /*
    std::string naif_path = spice_path + "/naif0012.tls";
    std::string de440_path = spice_path + "/de440.bsp";
    std::string pck_path = spice_path + "/pck00011.tpc";
    std::string earthfile_path = spice_path + "/earth_000101_250316_241218.bpc";

    // load SPICE kernels

    furnsh_c(naif_path.c_str());
    furnsh_c(de440_path.c_str());
    furnsh_c(pck_path.c_str());
    furnsh_c(earthfile_path.c_str());
    */
    loadAllKernels(spice_path);
    std::cout << "Done\n";

    // read and parse MPC ObsCodes
    std::vector<Observatory> obscodes = readObsCodes(obscode_path);
    std::unordered_map<std::string, Observatory> obscode_map;

    // mapify obscodes by... well... code
    for (int i = 0; i < obscodes.size(); i++)
    {
        obscode_map.insert({ obscodes[i].code, obscodes[i] });
    }

    std::vector<Observation> observations = readObsFile(obs_path);

    // DEBUG part
    // Debug commands are best written here, if you want to test out a function or sth.
    // END DEBUG part

    // main part, good luck
    std::vector<Vec3> state_vec = determineOrbit(observations, obscode_map, param_R, param_maxiter);
    Vec3 p0 = state_vec[0];
    Vec3 v0 = state_vec[1];

    // propagate the solution to a nearby midnight so that we have a good epoch
    std::pair<std::vector<Vec3>, double> final_solution = propagateToNextMidnight(p0, v0, observations[0].et);
    Vec3 p0_final = final_solution.first[0];
    Vec3 v0_final = final_solution.first[1];
    double epoch_JD = final_solution.second;

    // post-processing equivalent of an orbit fit problem
    std::cout << "Calculating final orbital parameters...\n";
    std::vector<double> orbital_elements = stateVector2Kepler(p0_final, v0_final, JDToEt(epoch_JD));
    // orbital_elements = {sma, eccentricity, inclination, omega, arg_periapsis, true_anomaly, mean_anomaly}
    // sma in km, angles in degrees
    
    // 1.3271244004193938E+11 = Sun GM
    double orbital_period = 2 * pi_c() * sqrt(orbital_elements[0] * orbital_elements[0] * orbital_elements[0] / 1.3271244004193938E+11);
    double orbit_traverse_fraction = orbital_elements[6] / 360;

    double periapsis_JD = epoch_JD - orbit_traverse_fraction * orbital_period / 86400.0;

    std::cout << "\n\n";
    std::cout << "================================================================================\n";
    printQuasiMPEC(orbital_elements, observations, obscode_map, epoch_JD, orbital_period, p0, v0);
    std::cout << "================================================================================\n";

    // now write to file
    std::cout << "\nWriting output quasi-MPEC to " << out_path << "... ";
    std::ofstream outfile(out_path);
    if (outfile.is_open())
    {
        printQuasiMPEC(orbital_elements, observations, obscode_map, epoch_JD, orbital_period, p0, v0, outfile);
        outfile.close();
    }
    std::cout << "Done. Program end.\n";
}
