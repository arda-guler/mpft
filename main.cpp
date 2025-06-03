#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>

// #include "Eigen/Dense"

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
std::string trim(const std::string& str) // trim trailing/leading whitespace, nothing fancy, stole it from somewhere
{
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// this is the fun and simpler part of the equation
std::vector<double> stateVector2Kepler(Vec3 r, Vec3 v, double mu = 1.3271244004193938E+11)
{
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

/*
std::vector<std::complex<double>> solvePolynomial(const std::vector<double>& coefficients) {
    int degree = coefficients.size() - 1;

    if (degree < 1 || coefficients[0] == 0.0) {
        return {};
    }

    // Normalize coefficients to make polynomial monic
    std::vector<double> normalized_coeffs(degree + 1);
    double leading = coefficients[0];
    for (int i = 0; i <= degree; ++i) {
        normalized_coeffs[i] = coefficients[i] / leading;
    }

    // Build companion matrix
    Eigen::MatrixXd companion = Eigen::MatrixXd::Zero(degree, degree);
    for (int i = 1; i < degree; ++i) {
        companion(i - 1, i) = 1.0;  // subdiagonal
    }

    for (int i = 0; i < degree; ++i) {
        companion(degree - 1, i) = -normalized_coeffs[i + 1];
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd eigenvalues = solver.eigenvalues();

    std::vector<std::complex<double>> roots(degree);
    for (int i = 0; i < degree; ++i) {
        roots[i] = eigenvalues[i];
    }

    return roots;
}
*/

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

// I could've made a KeplerOrbit class or whatever but this works too
void printOrbitalElements(std::vector<double> orbel)
{
    std::cout << "a: " << orbel[0] / AU << " e: " << orbel[1] << " i: " << orbel[2] << " Node: " << orbel[3] << " Peri: " << orbel[4] << " M: " << orbel[5] << "\n";
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

    std::cout << "Done.\n";
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
        dt = time_interval / 16;
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
// --- --- --- ORBIT PROPAGATION --- --- --- [END]

/*
std::vector<double> get_rho2s(Observation o1, Observation o2, Observation o3,
    std::unordered_map<std::string, Observatory> &obscode_map)
{
    double mu = 1.3271244004193938E+11;

    Vec3 rhat_1 = Vec3(o1.RA, o1.DEC);
    Vec3 rhat_2 = Vec3(o2.RA, o2.DEC);
    Vec3 rhat_3 = Vec3(o3.RA, o3.DEC);

    Vec3 R_1 = getObserverPos(obscode_map[o1.obs_code], o1.et);
    Vec3 R_2 = getObserverPos(obscode_map[o2.obs_code], o2.et);
    Vec3 R_3 = getObserverPos(obscode_map[o3.obs_code], o3.et);

    double t3mt1 = o3.et - o1.et;
    Vec3 rhatprime_2 = (rhat_3 - rhat_1) / t3mt1;
    Vec3 rhatprimeprime_2 = (rhat_3 - rhat_2 * 2 + rhat_1) / (t3mt1 * t3mt1);

    double tau_1 = o1.et - o2.et;
    double tau_3 = o3.et - o2.et;
    double tau = o3.et - o1.et;

    Vec3 p_1 = rhat_2.cross(rhat_3);
    Vec3 p_2 = rhat_1.cross(rhat_3);
    Vec3 p_3 = rhat_1.cross(rhat_2);

    double D_0 = rhat_1.dot(p_1);
    double D_11 = R_1.dot(p_1);
    double D_12 = R_1.dot(p_2);
    double D_13 = R_1.dot(p_3);
    double D_21 = R_2.dot(p_1);
    double D_22 = R_2.dot(p_2);
    double D_23 = R_2.dot(p_3);
    double D_31 = R_3.dot(p_1);
    double D_32 = R_3.dot(p_2);
    double D_33 = R_3.dot(p_3);

    double A = 1 / D_0 * (-D_12 * tau_3 / tau + D_22 + D_32 * tau_1 / tau);
    double B = 1 / (6 * D_0) * (D_12 * (tau_3 * tau_3 - tau * tau) * tau_3 / tau + D_32 * (tau * tau - tau_1 * tau_1) * tau_1 / tau);
    double E = R_2.dot(rhat_2);

    double R_2sq = R_2.dot(R_2);

    double a = -(A * A + 2 * A * E + R_2sq);
    double b = -2 * mu * B * (A + E);
    double c = -mu * mu * B * B;

    std::vector<double> coeffs = { 1, 0, a, 0, 0, b, 0, 0, c };

    std::cout << "COEFFS: ";
    for (int aa = 0; aa < coeffs.size(); aa++)
    {
        std::cout << coeffs[aa] << " ";
    }
    std::cout << "\n";

    std::vector<std::complex<double>> rho_2s = solvePolynomial(coeffs);

    std::cout << "RHO: ";
    for (int aa = 0; aa < rho_2s.size(); aa++)
    {
        std::cout << rho_2s[aa] << " ";
    }
    std::cout << "\n";

    std::vector<double> rho_2s_sanitized = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < rho_2s.size(); i++)
    {
        rho_2s_sanitized[i] = abs(rho_2s[i].real());
    }

    return rho_2s_sanitized;
}
*/

// the actual orbit determination algorithm
std::vector<Vec3> determineOrbit(std::vector<Observation> obs_all, std::unordered_map<std::string, Observatory>& obscode_map,
    double R1, double deltaR)
{
    double mu = 1.3271244004193938E+11;

    std::cout << "Orbit determination...\n";

    Observation o1 = obs_all[0];
    Observation o2 = obs_all[(int)(obs_all.size() / 2)];
    Observation o3 = obs_all[obs_all.size() - 1];

    SpiceDouble date_init = o1.et;
    SpiceDouble date_final = o2.et;
    SpiceDouble date_check = o3.et;

    // std::vector<double> rho_2s = get_rho2s(o1, o2, o3, obscode_map);
    
    // double R1 = 3 * AU; // *std::max_element(rho_2s.begin(), rho_2s.end());
    // double deltaR = 0;
    double R2 = R1 + deltaR;

    Vec3 p0 = getObserverPos(obscode_map[o1.obs_code], o1.et) + Vec3(o1.RA, o1.DEC) * R1;
    Vec3 pf = getObserverPos(obscode_map[o2.obs_code], o2.et) + Vec3(o2.RA, o2.DEC) * R2;

    Vec3 v0 = Vec3(-p0.y / p0.mag(), p0.x / p0.mag(), 0) * sqrt(mu / p0.mag());

    std::vector<double> orbital_elems_init = stateVector2Kepler(p0, v0);
    std::cout << "Initial guess:\n";
    printOrbitalElements(orbital_elems_init);
    std::cout << "\nBeginning iterations...\n";

    double adjust_factor = 1;
    double adjust_pfactor = 1;

    bool good_fit = false;
    int retry_count = 0;
    int MAX_RETRY = 100;
    while (!good_fit && retry_count <= MAX_RETRY)
    {
        // check current error
        double err_val = 0;
        for (int idx_o = 0; idx_o < obs_all.size(); idx_o++)
        {
            Observation o = obs_all[idx_o];

            // p0.printout();

            std::vector<Vec3> prop_res = propagate(p0, v0, o1.et, o.et, -1);
            Vec3 p_check = prop_res[0];
            std::vector<double> RADEC_prop = getRADEC(o.et, obscode_map[o.obs_code], p_check);

            double RA_prop = RADEC_prop[0];
            double DEC_prop = RADEC_prop[1];

            double RA_err = o.RA - RA_prop;
            double DEC_err = o.DEC - DEC_prop;

            err_val += sqrt(RA_err * RA_err + DEC_err * DEC_err);
        }

        // adjust p0, v0
        if (err_val < 10 / 3600)
        {
            good_fit = true;
        }
        else
        {
            std::vector<double> adjust_vals = { 0, 0, 0, 0, 0, 0 };
            std::vector<Vec3> adjust_vecs =
            {
                Vec3(0.01,  0,      0),
                Vec3(-0.01, 0,      0),
                Vec3(0,     0.01,   0),
                Vec3(0,     -0.01,  0),
                Vec3(0,     0,      0.01),
                Vec3(0,     0,      -0.01),
            };

            for (int idx_dirad = 0; idx_dirad < 6; idx_dirad++)
            {
                adjust_vecs[idx_dirad] = adjust_vecs[idx_dirad] * adjust_factor;
            }

            // adjust vel
            for (int idx_dir = 0; idx_dir < 6; idx_dir++)
            {
                Vec3 v0_1 = v0 + adjust_vecs[idx_dir];

                adjust_vals[idx_dir] = 0;
                for (int idx_oj = 0; idx_oj < obs_all.size(); idx_oj++)
                {
                    Observation oj = obs_all[idx_oj];
                    if (oj.et != date_init)
                    {
                        std::vector<Vec3> prop_res = propagate(p0, v0_1, o1.et, oj.et, -1);
                        Vec3 p_check = prop_res[0];
                        std::vector<double> RADEC_prop_1 = getRADEC(oj.et, obscode_map[oj.obs_code], p_check);

                        double RA_prop_1 = RADEC_prop_1[0];
                        double DEC_prop_1 = RADEC_prop_1[1];

                        double RA_err_1 = oj.RA - RA_prop_1;
                        double DEC_err_1 = oj.DEC - DEC_prop_1;

                        adjust_vals[idx_dir] += sqrt(RA_err_1 * RA_err_1 + DEC_err_1 * DEC_err_1);
                    }
                }
            }

            std::vector<double>::iterator it = std::min_element(adjust_vals.begin(), adjust_vals.end());
            int idx_min = std::distance(std::begin(adjust_vals), it);
            if (adjust_vals[idx_min] < err_val)
            {
                v0 = v0 + adjust_vecs[idx_min];
                adjust_factor *= 2;
            }
            else
            {
                adjust_factor *= 0.1;
            }

            // adjust pos
            std::vector<double> adjust_pvals = { 0, 0, 0, 0, 0, 0 };
            std::vector<Vec3> adjust_pvecs =
            {
                Vec3(0.01,  0,      0),
                Vec3(-0.01, 0,      0),
                Vec3(0,     0.01,   0),
                Vec3(0,     -0.01,  0),
                Vec3(0,     0,      0.01),
                Vec3(0,     0,      -0.01),
            };

            for (int idx_pdirad = 0; idx_pdirad < 6; idx_pdirad++)
            {
                adjust_pvecs[idx_pdirad] = adjust_pvecs[idx_pdirad] * adjust_pfactor;
            }

            for (int idx_pdir = 0; idx_pdir < 6; idx_pdir++)
            {
                Vec3 p0_1 = p0 + adjust_pvecs[idx_pdir];

                adjust_pvals[idx_pdir] = 0;
                for (int idx_ok = 0; idx_ok < obs_all.size(); idx_ok++)
                {
                    Observation ok = obs_all[idx_ok];
                    if (ok.et != date_init)
                    {
                        std::vector<Vec3> prop_res = propagate(p0_1, v0, o1.et, ok.et, -1);
                        Vec3 p_check = prop_res[0];
                        std::vector<double> RADEC_prop_2 = getRADEC(ok.et, obscode_map[ok.obs_code], p_check);

                        double RA_prop_2 = RADEC_prop_2[0];
                        double DEC_prop_2 = RADEC_prop_2[1];

                        double RA_err_2 = ok.RA - RA_prop_2;
                        double DEC_err_2 = ok.DEC - DEC_prop_2;

                        adjust_pvals[idx_pdir] += sqrt(RA_err_2 * RA_err_2 + DEC_err_2 * DEC_err_2);
                    }
                }
            }

            std::vector<double>::iterator pit = std::min_element(adjust_pvals.begin(), adjust_pvals.end());
            int idx_pmin = std::distance(std::begin(adjust_pvals), pit);
            if (adjust_pvals[idx_pmin] < err_val)
            {
                p0 = p0 + adjust_pvecs[idx_pmin];
                adjust_pfactor *= 2;
            }
            else
            {
                adjust_pfactor *= 0.1;
            }
        }

        retry_count += 1;
    }

    std::cout << "Done iterating.\n";
    return std::vector<Vec3> {p0, v0};
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
    std::cout << "-obs <Obs80_input_file> -obscode <obscode_file> -spice <spice_folder_path> -R <initial_object_dist_to_observer_guess (AU)>\n\n";
    std::cout << "MPFT was developed by H. A. Guler.\n\n";
    std::cout << "MPFT is licensed under GNU General Public License version 2.0 (GPL-2.0 License)\n\n";
}

int main(int argc, char *argv[])
{
    std::cout << "MPFT v0.1.0\n";

    // default parameters
    std::string obs_path = "primary.obs";
    int param_R = 3 * AU;
    int param_deltaR = 0;
    std::string obscode_path = "data/ObsCodes.dat";
    std::string spice_path = "data/SPICE/";

    // handle command line arguments
    // there is a more compact version of doing this but this is easier for my brain
    int argtype = 0;
    for (int idx_cmd = 0; idx_cmd < argc; idx_cmd++)
    {
        if (strcmp(argv[idx_cmd], "-obs") || strcmp(argv[idx_cmd], "-observations")) // also the default argument
        {
            argtype = 0;
        }
        else if (strcmp(argv[idx_cmd], "-R"))
        {
            argtype = 1;
        }
        else if (strcmp(argv[idx_cmd], "-deltaR"))
        {
            argtype = 2;
        }
        else if (strcmp(argv[idx_cmd], "-obscode"))
        {
            argtype = 3;
        }
        else if (strcmp(argv[idx_cmd], "-spice"))
        {
            argtype = 4;
        }
        else if (strcmp(argv[idx_cmd], "-h") || strcmp(argv[idx_cmd], "--help"))
            printHelpMsg();
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
                param_deltaR = strtod(argv[idx_cmd], NULL);
                break;
            case 3:
                obscode_path = argv[idx_cmd];
                break;
            case 4:
                spice_path = argv[idx_cmd];
                break;
            }
        }
    }

    std::cout << "Loading SPICE kernels... ";
    std::string naif_path = spice_path + "/naif0012.tls";
    std::string de440_path = spice_path + "/de440.bsp";
    std::string pck_path = spice_path + "/pck00011.tpc";
    std::string earthfile_path = spice_path + "/earth_000101_250316_241218.bpc";

    // load SPICE kernels
    furnsh_c(naif_path.c_str());
    furnsh_c(de440_path.c_str());
    furnsh_c(pck_path.c_str());
    furnsh_c(earthfile_path.c_str());
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
    std::vector<Vec3> state_vec = determineOrbit(observations, obscode_map, param_R, param_deltaR);
    Vec3 p0 = state_vec[0];
    Vec3 v0 = state_vec[1];

    // post-processing equivalent of an orbit fit problem
    std::cout << "Calculating final orbital parameters...\n";
    std::vector<double> orbital_elements = stateVector2Kepler(p0, v0);

    std::cout << "Final orbital elements:\n";
    printOrbitalElements(orbital_elements);
    std::cout << "\nMPFT: Program end.\n";
}