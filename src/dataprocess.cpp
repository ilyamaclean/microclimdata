#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath> 
#include <limits>
#include <unordered_set>
using namespace Rcpp;
// Function to perform linear interpolation
// [[Rcpp::export]]
NumericVector na_approx(NumericVector v) {
    NumericVector vec = clone(v); // Copy the input vector
    if (v.size() > 0) {
        vec.push_back(v[0]); // Append the first element to the end
    }
    NumericVector result = clone(vec);
    int n = result.size();
    // Helper lambda to find the next non-NaN value
    auto find_next_valid = [&](int start) -> int {
        for (int i = start; i < n; ++i) {
            if (!R_IsNA(result[i]) && !R_IsNaN(result[i])) return i;
        }
        return -1;
    };
    int last_valid = -1;
    for (int i = 0; i < n; ++i) {
        if (!R_IsNA(result[i]) && !R_IsNaN(result[i])) {
            last_valid = i;
        }
        else {
            int next_valid = find_next_valid(i);
            if (next_valid == -1) {
                stop("Cannot interpolate: trailing NaNs at the end of the vector.");
            }
            if (last_valid == -1) {
                stop("Cannot interpolate: leading NaNs at the beginning of the vector.");
            }
            // Linear interpolation
            double step = (result[next_valid] - result[last_valid]) / (next_valid - last_valid);
            for (int j = last_valid + 1; j < next_valid; ++j) {
                result[j] = result[last_valid] + step * (j - last_valid);
            }
            i = next_valid - 1; // Move to the last valid position found
        }
    }
    NumericVector res2(n - 1);
    for (int i = 0; i < (n - 1); ++i) {
        res2[i] = result[i];
    }
    return res2;
}
// Function to compute the spline coefficients
void spline_coefficients(const NumericVector& x, const NumericVector& y,
    NumericVector& a, NumericVector& b, NumericVector& c, NumericVector& d) {
    int n = x.size() - 1;
    NumericVector h(n), alpha(n), l(n + 1), mu(n + 1), z(n + 1);
    // Step 1: Calculate h and alpha
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
        alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);
    }
    // Step 2: Solve the tridiagonal system
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    for (int i = 1; i < n; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n] = 1.0;
    z[n] = 0.0;
    c[n] = 0.0;
    // Step 3: Calculate coefficients
    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
        a[j] = y[j];
    }
}
// Spline interpolation function
double spline_interpolate_point(double xi, const NumericVector& x, const NumericVector& a,
    const NumericVector& b, const NumericVector& c, const NumericVector& d) {
    int n = x.size() - 1;
    int i = n - 1;
    while (i >= 0 && xi < x[i]) i--;
    double dx = xi - x[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}
// [[Rcpp::export]]
NumericVector splineCpp(NumericVector daily_values) {
    int days = daily_values.size();
    int hours_per_day = 24;
    // Extend the data by duplicating the first and last day
    NumericVector extended_values(days + 2);
    extended_values[0] = daily_values[0]; // Duplicate first day
    for (int i = 0; i < days; ++i) {
        extended_values[i + 1] = daily_values[i];
    }
    extended_values[days + 1] = daily_values[days - 1]; // Duplicate last day
    // Create x-axis (day indices) for the extended daily values (-1, 0, 1, ..., days)
    NumericVector x(days + 2);
    for (int i = 0; i < days + 2; ++i) {
        x[i] = i - 1; // Midday values are at -1, 0, ..., days
    }
    // Create x-axis for hourly intervals (fractional day values)
    NumericVector xout((days + 2) * hours_per_day);
    for (int h = 0; h < xout.size(); ++h) {
        xout[h] = h / (double)hours_per_day - 1; // Starts from day -1 to days
    }
    // Coefficients for cubic spline
    NumericVector a(days + 1), b(days + 1), c(days + 2), d(days + 1);
    // Compute spline coefficients
    spline_coefficients(x, extended_values, a, b, c, d);
    // Interpolate the values at the hourly intervals
    NumericVector interpolated_values(xout.size());
    for (int i = 0; i < xout.size(); ++i) {
        interpolated_values[i] = spline_interpolate_point(xout[i], x, a, b, c, d);
    }
    // Crop the first and last 12 hours (corresponding to the duplicated days)
    int num_hours = days * hours_per_day;
    NumericVector final_values = interpolated_values[Range(12, num_hours + 11)];
    return final_values;
}
NumericVector applymonth(NumericVector pai, NumericVector me) {
    // Check if input vectors have the same length
    if (pai.size() != me.size()) {
        stop("Vectors pai and me must have the same length.");
    }
    // Calculate mean pai and month effect
    double sumpai = 0;
    double summe = 0;
    int n = 0;
    for (int i = 0; i < pai.size(); ++i) {
        if (!R_IsNA(pai[i]) && !R_IsNaN(pai[i])) {
            sumpai += pai[i];
            summe += me[i];
            ++n;
        }
    }
    if (n == 0) {
        stop("No non-NA values in pai.");
    }
    double mepai = sumpai / static_cast<double>(n);
    double meme = summe / static_cast<double>(n);
    // Calculate mu: mean of non-NA pai / mean of non-NA month effects
    double mu = mepai / meme;
    // Apply month effects if pai is NA
    for (int i = 0; i < pai.size(); ++i) {
        if (R_IsNA(pai[i]) || R_IsNaN(pai[i])) {
            pai[i] = me[i] * mu;
        }
    }
    return pai;
}
int countMaxConsecutiveNAs(NumericVector vec) {
    int max_consecutive_nas = 0;
    int current_consecutive_nas = 0;
    for (double val : vec) {
        if (R_IsNA(val) || R_IsNaN(val)) {
            current_consecutive_nas++;
            if (current_consecutive_nas > max_consecutive_nas) {
                max_consecutive_nas = current_consecutive_nas;
            }
        }
        else {
            current_consecutive_nas = 0;
        }
    }
    return max_consecutive_nas;
}
NumericVector decidena(NumericVector pai, NumericVector me, int n = 3) {
    int cna = countMaxConsecutiveNAs(pai);
    if (cna > n) {
        pai = applymonth(pai, me);
    }
    else {
        pai = na_approx(pai);
    }
    return pai;
}
NumericVector aperm3D(NumericVector Tz, int rows, int cols, int tsteps) {
    // Initialize a new vector to store the permuted array
    NumericVector permuted(Tz.size());
    // Permute the dimensions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                // Calculate the index in the original array
                int index_original = k + tsteps * (j + cols * i);
                // Calculate the index in the permuted array
                int index_permuted = i + rows * (j + cols * k);
                // Copy the element from the original array to the permuted array
                permuted[index_permuted] = Tz[index_original];
            }
        }
    }
    // Set the dimensions of the permuted array
    permuted.attr("dim") = IntegerVector::create(rows, cols, tsteps);
    return permuted;
}
// Function to permute dimensions of a 3D array and convert to a NumericVector (the other way)
NumericVector aperm3D2(NumericVector Tz, int rows, int cols, int tsteps)
{
    // Initialize a new vector to store the permuted array
    NumericVector permuted(Tz.size());
    // Permute the dimensions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                // Calculate the index in the original array
                int index_original = i + rows * (j + cols * k);
                // Calculate the index in the permuted array
                int index_permuted = k + tsteps * (j + cols * i);
                // Copy the element from the original array to the permuted array
                permuted[index_permuted] = Tz[index_original];
            }
        }
    }
    // Set the dimensions of the permuted array
    return permuted;
}
// [[Rcpp::export]]
NumericVector napproxCpp(NumericVector pai, NumericMatrix landcover, NumericMatrix lpai)
{
    IntegerVector dims = pai.attr("dim");
    // Extract rows, columns, and layers
    int rows = dims[0];
    int cols = dims[1];
    int months = dims[2];
    pai = aperm3D2(pai, rows, cols, months);
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = pai[index];
            if (!NumericVector::is_na(val)) {
                int st = j * months + i * months * cols;
                int ed = st + months - 1;
                // subset to vector
                NumericVector v = pai[Range(st, ed)];
                // determine land cover class
                int s = landcover(i, j) - 1;
                if (s >= 0) {
                    NumericVector me = lpai(s, _);
                    NumericVector v2 = decidena(v, me);
                    for (int k = 0; k < months; ++k) {
                        pai[index] = v2[k];
                        ++index;
                    }
                }
                else {
                    for (int k = 0; k < months; ++k) ++index;
                }
            }
            else {
                for (int k = 0; k < months; ++k) ++index;
            }
        } // end j
    } // end i
    pai = aperm3D(pai, rows, cols, months);
    return pai;
}
// [[Rcpp::export]]
IntegerVector uniquecpp(const IntegerMatrix& matrix) {
    std::set<int> uniqueIntegers;
    // Get the dimensions of the matrix
    int nrow = matrix.nrow();
    int ncol = matrix.ncol();
    // Traverse the matrix and insert elements into the set
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int value = matrix(i, j);
            if (value != NA_INTEGER) {
                uniqueIntegers.insert(value);
            }
        }
    }
    // Convert the set to an IntegerVector
    IntegerVector result(uniqueIntegers.begin(), uniqueIntegers.end());
    return result;
}
// C++ code for calculating mean pai of a given land cover class
// [[Rcpp::export]]
NumericVector meanpai(NumericMatrix pai, IntegerMatrix landcover)
{
    int nrow = pai.nrow();
    int ncol = pai.ncol();
    IntegerVector lcclass = uniquecpp(landcover);
    NumericVector mpai(lcclass.size());
    for (R_xlen_t l = 0; l < lcclass.size(); ++l) {
        int lc = lcclass[l];
        double psum = 0;
        int pc = 0;
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < ncol; ++j) {
                double val = pai(i, j);
                int val2 = landcover(i, j);
                if (val2 == lc && !NumericMatrix::is_na(val)) {
                    psum = psum + val;
                    ++pc;
                }
            }
        }
        mpai[l] = psum / static_cast<double>(pc);
    }
    return mpai;
}
// Function to find the index of an element in a vector that is equal to a given value
// [[Rcpp::export]]
int whichcpp(IntegerVector v, int x) {
    for (R_xlen_t i = 0; i < v.size(); ++i) {
        if (v[i] == x) {
            return i; // Return the index if found
        }
    }
    return -1; // Return -1 if not found
}
// C++ code for filling pai by land cover
// [[Rcpp::export]]
NumericMatrix fillpai(NumericMatrix pai, IntegerMatrix landcover)
{
    // get mean pai for each lcover class
    IntegerVector lcclass = uniquecpp(landcover);
    NumericVector mpai = meanpai(pai, landcover);
    int nrow = pai.nrow();
    int ncol = pai.ncol();
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            double val = pai(i, j);
            if (Rcpp::NumericMatrix::is_na(val)) {
                int val2 = landcover(i, j);
                if (!Rcpp::IntegerMatrix::is_na(val2)) {
                    // get entry of land cover class
                    int s = whichcpp(lcclass, val2);
                    pai(i, j) = mpai[s];
                }
            }
        }
    }
    return pai;
}
// [[Rcpp::export]]
NumericMatrix expectedpai(IntegerMatrix landcover, IntegerMatrix landcoverm, NumericMatrix pai)
{
    // get vector of meanpais
    IntegerVector lcclass = uniquecpp(landcoverm);
    NumericVector mpai = meanpai(pai, landcoverm);
    int nrow = pai.nrow();
    int ncol = pai.ncol();
    NumericMatrix epai(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            int val = landcover(i, j);
            if (!Rcpp::IntegerMatrix::is_na(val)) {
                int s = whichcpp(lcclass, val);
                epai(i, j) = mpai[s];
            }
            else {
                epai(i, j) = NA_REAL;
            }
        }
    }
    return epai;
}
// Functions used for seasonal fill
// [[Rcpp::export]]
NumericVector seasoneffect(NumericVector lai) {
    // Get the dimensions of the 3D array
    IntegerVector dims = lai.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int mths = dims[2];
    // Output vector to store the means of each slice
    NumericVector means(mths);
    // Loop through each slice (3rd dimension)
    for (int k = 0; k < mths; ++k) {
        double sum = 0;
        int count = 0;
        // Loop through each element in the slice
        for (int i = 0; i < rows * cols; ++i) {
            double value = lai[k * rows * cols + i];
            if (!NumericVector::is_na(value)) {
                sum += value;
                count++;
            }
        }
        // If there are valid values, compute the mean, otherwise set NA
        if (count > 0) {
            means[k] = sum / count;
        }
        else {
            means[k] = NA_REAL;
        }
    }
    return means;
}
// Functions used for seasonal fill
// [[Rcpp::export]]
NumericVector seasonadjCpp(NumericVector lai, NumericVector seffect)
{
    IntegerVector dims = lai.attr("dim");
    // Extract rows, columns, and layers
    int rows = dims[0];
    int cols = dims[1];
    int months = dims[2];
    lai = aperm3D2(lai, rows, cols, months);
    int index1 = 0;
    int index2 = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Extract vector
            NumericVector v(months);
            for (int k = 0; k < months; ++k) {
                v[k] = lai[index1];
                ++index1;
            }
            // Check whther to interpolate or whether to just apply month effect
            int test = 1;
            if (NumericVector::is_na(v[0])) {
                test = 0;
            }
            else {
                int count = 0;
                for (int k = 0; k < months; ++k) {
                    if (!NumericVector::is_na(v[k])) ++count;
                }
                double frac = static_cast<double>(count) / months;
                if (frac < 0.5) test = 0;
            }
            // if test equals one interpolate
            if (test > 0) {
                NumericVector vo = na_approx(v);
                for (int k = 0; k < months; ++k) {
                    lai[index2] = vo[k];
                    ++index2;
                }
            }
            // if test equals zero apply month effect
            else {
                // ** Calculate mean of month effect
                double sum1 = 0.0;
                double sum2 = 0.0;
                int count = 0;
                for (int k = 0; k < months; ++k) {
                    sum1 = sum1 + seffect[k];
                    if (!NumericVector::is_na(v[k])) {
                        sum2 = sum2 + v[k];
                        ++count;
                    }
                }
                double mean1 = sum1 / months;
                double mean2 = sum2 / count;
                for (int k = 0; k < months; ++k) {
                    if (NumericVector::is_na(v[k])) {
                        lai[index2] = (seffect[k] / mean1) * mean2;
                    }
                    ++index2;
                } // end k
            } // end else
        } // end j
    } // end i
    lai = aperm3D(lai, rows, cols, months);
    return lai;
}
// ******************************************************************** //
// ~~~~~~~~~~~~~~~~~ Climate variable interpolations ~~~~~~~~~~~~~~~~~~ //
// ******************************************************************** //
// Calculates astronomical julian day
// [[Rcpp::export]]
int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// ** Calculates solar time ** //
double soltimeCpp(int jd, double lt, double lon)
{
    double m = 6.24004077 + 0.01720197 * (jd - 2451545);
    double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
    double st = lt + (4 * lon + eot) / 60;
    return st;
}
// ** Calculates solar zenith in radians ** //
double solzenCpp(double jd, double st, double lat, double lon)
{
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh);
    return z;
}
// ** Calculates day length ** //
double daylengthCpp(int jd, double lat)
{
    double declin = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double latr = lat * M_PI / 180;
    double hc = -0.01453808 / (cos(latr) * cos(declin)) - tan(latr) * tan(declin);
    double dl = 0;
    if (hc < -1) {
        dl = 24;
    }
    else if (hc < 1) {
        double ha = (acos(hc)) * 180 / M_PI;
        double m = 6.24004077 + 0.01720197 * (jd - 2451545);
        double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
        double sr = (720 - 4 * ha - eot) / 60;
        double ss = (720 + 4 * ha - eot) / 60;
        dl = ss - sr;
    }
    return(dl);
}
// ** Calculates 24 hourly temperatures for daily data
NumericVector tempintdayCpp(double tmn, double tmnn, double tmx, double dl, double stt, double lat, double lon, double srte = 0.09)
{
    // Calculate predicted night fraction
    double ngtp = 0.04187957 * ((tmx - tmn) * (1 - dl / 24)) + 0.4372056;
    if (ngtp < 0.01) ngtp = 0.01;
    if (ngtp > 0.99) ngtp = 0.99;
    // Calculate sunrise time
    double sr = 12 - 0.5 * dl;
    // Calculate solar time after sunrise
    NumericVector thour(24);
    for (int lt = 0; lt < 24; lt++) {
        double st = lt + stt - sr; // solar time after sunrise
        if (st < 0) st = st + 24;
        if (st > 24) st = st - 24;
        if (dl == 24) {
            thour[lt] = (tmx - tmn) * sin((M_PI * st) / 28) + tmn;
            double gr = (tmnn - tmn) * (st / 24);
            thour[lt] = thour[lt] + gr;
        }
        else if (dl == 0) {
            st = lt + stt; // solar time after sunrise
            if (st < 0) st = st + 24;
            if (st > 24) st = st - 24;
            thour[lt] = (tmx - tmn) * sin((M_PI * st) / 24) + tmn;
            double gr = (tmnn - tmn) * (st / 24);
            thour[lt] = thour[lt] + gr;
        }
        else {
            double k = -(24 - dl) / log(srte / ngtp);
            double ph = -0.5 * dl * ((M_PI / (asin(ngtp) - M_PI)) + 1);
            double rho = dl + 2 * ph;
            if (st > dl) {
                thour[lt] = ngtp * exp(-(st - dl) / k);
            }
            else {
                thour[lt] = sin((M_PI * st) / rho);
            }
            // Adjust by tmax and tmin
            thour[lt] = (tmx - tmn) * thour[lt] + tmn;
            // Apply gradient to times after sunset as tmn is next day
            if (lt + stt > dl) {
                double gr = (tmnn - tmn) * (st - dl) / (24 - dl);
                thour[lt] = thour[lt] + gr;
            }
        }
    }
    // Adjust to ensure tmax and tmin match
    double ptmx = thour[0];
    double ptmn = thour[0];
    for (int lt = 1; lt < 24; lt++) {
        if (thour[lt] > ptmx) ptmx = thour[lt];
        if (thour[lt] < ptmn) ptmn = thour[lt];
    }
    double b = (ptmx - ptmn) / (tmx - tmn);
    double a = b * tmn - ptmn;
    for (int lt = 0; lt < 24; lt++) thour[lt] = (thour[lt] + a) / b;
    return(thour);
}
// ** Interpolates hourly temperature (vector) ** //
// [[Rcpp::export]]
NumericVector hourlytempv(NumericVector tmn, NumericVector tmx,
    IntegerVector year, IntegerVector month, IntegerVector day,
    double lat, double lon, double srte = 0.09)
{
    int nrow = tmn.size();
    NumericMatrix thour(nrow, 24);
    for (int i = 0; i < nrow; ++i) {
        // Calculate jd
        int jd = juldayCpp(year[i], month[i], day[i]);
        // Calculate sunrise time
        double dl = daylengthCpp(jd, lat);
        double stt = soltimeCpp(jd, 0, lon);
        double sr = stt + 12 - 0.5 * dl;
        if (sr > 24) sr = sr - 24;
        if (sr > 0) { // Sunrise in current day
            double tmnn = tmn[i + 1];
            if (i == (nrow - 1)) tmnn = tmn[nrow - 1];
            thour(i, _) = tempintdayCpp(tmn[i], tmnn, tmx[i], dl, stt, lat, lon, srte);
        }
        else { // sunrise in previous day
            double tmnp = tmn[0];
            if (i > 0) tmnp = tmn[i - 1];
            thour(i, _) = tempintdayCpp(tmnp, tmn[i], tmx[i], dl, stt, lat, lon, srte);
        }
    }
    NumericVector thourv(nrow*24);
    // Flatten the matrix row-wise
    int index = 0;
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < 24; ++j) {
            thourv[index] = thour(i, j);
            ++index;
        }
    }
    return thourv;
}
// ** converts humidity
// [[Rcpp::export]]
NumericVector converthumidityCpp(NumericVector sph, NumericVector tc, NumericVector pk)
{
    NumericVector rh(sph.size());
    for (R_xlen_t i = 0; i < sph.size(); ++i) {
        double es = 0.6108 * exp(17.27 * tc[i] / (tc[i] + 237.3));
        double ea = (sph[i] * pk[i]) / (0.622 + (0.378 * sph[i]));
        rh[i] = (ea / es) * 100.0;
        if (rh[i] > 100.0) rh[i] = 100.0;
        if (rh[i] < 15.0) rh[i] = 15.0;
    }
    return rh;
}
// ************************************************************************* //
// ~~~~~~~~~~~~~ Calculate clearsky radiation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************************************* //
// Calculate solar altitude corrected for solar refraction
double solarref(double Z) {
    // Convert to solar altitude etc
    double Zd = (180.0 / M_PI) * Z;
    double Ad = 90.0 - Zd;
    double A = (M_PI / 180.0) * Ad;
    double top = 0.1594 + 1.123 * A + 0.065656 * pow(A, 2);
    double btm = 1.0 + 28.9344 * A + 277.3971 * pow(A, 2);
    double dAR = 0.061359 * (top / btm);
    double At = A + dAR;
    return At;
}
// Calculate Kasten and Youing airmass coefficient
// [[Rcpp::export]]
double airMass(double Z, double z = 0.0) {
    // Derive refraction corrected solar altitude
    double At = solarref(Z); 
    double At_deg = At * 180.0 / M_PI;  // Convert radians to degrees
    double altc = exp(-z / 8434.5);
    double btm = sin(At) + 0.50572 * pow(At_deg + 6.07995, -1.6364);
    double m = altc / btm;
    return m;
}
// Calculate Rayleigh optical thickness
double Rayleigh(double m) {
    double ir = 10.4 + 0.718 * m;
    if (m <= 20) {
        ir = 6.62960 + 1.7513 * m - 0.1202 * pow(m, 2) + 0.0065 * pow(m, 3) - 0.00013 * pow(m, 4);
    }
    return 1 / ir;
}
// Calculate Linke turbidity (sh = specific humidity in kg/kg)
double Linke(double sh, double pk = 101.3, double aod = 0.1)
{
    // Integrated precipitable water vapor content of the atmosphere
    double w = sh * 101.3 * 1000.0 / 9.81;
    // Optical water depth
    double owd = 0.07649825 * pow(w, 0.34);
    double TL = 11.2 * (0.1093309 + owd + 0.1);
    return TL;
}
// Calculate direct (beam) radiation tranmission
// [[Rcpp::export]]
double dirtran(double Z, double TL, double z = 0)
{
    double sa_t = solarref(Z);
    double Tr = 0.0;
    if (sa_t > 0.0) {
        double m = airMass(Z, z);
        double r = Rayleigh(m);
        Tr = exp(-0.8662 * TL * m * r);
    }
    return Tr;
}
// Calculate diffuse radiation tranmission
// [[Rcpp::export]]
double diftran(double Z, double TL)
{
    double Td = -0.015843 + 0.030543 * TL + 0.0003797 * TL * TL;
    double a0 = 0.26463 - 0.061581 * TL + 0.0031408 * TL * TL;
    if (a0 < 0.02) a0 = 0.02 / Td;
    double a1 = 2.0402 + 0.018945 * TL - 0.011161 * TL * TL;
    double a2 = -1.3025 + 0.039231 * TL + 0.0085079 * TL * TL;
    double Fd = a0 + a1 * cos(Z) + a2 * cos(Z) * cos(Z);
    double Tr = Td * Fd;
    if (Tr < 0.0) Tr = 0.0;
    return Tr;
}
// ** Calculate extra terrestrial radiation
// [[Rcpp::export]]
double extrarad(int jd) {
    // Calculate the day of the year
    double dayOfYear = jd - 2451545.0 + 0.0008; // Julian Day for Jan 1, 2000
    // Compute the number of days since the start of the year
    double N = dayOfYear;
    // Calculate the Earth's eccentricity factor
    double eccentricity = 1.00011 + 0.034221 * cos(2 * M_PI * N / 365.0)
        + 0.00128 * sin(2 * M_PI * N / 365.0)
        + 0.000719 * cos(4 * M_PI * N / 365.0)
        + 0.000077 * sin(4 * M_PI * N / 365.0);
    // Compute the extraterrestrial radiation
    double extraterrestrialRad = 1361.0 * eccentricity;
    return extraterrestrialRad;
}
// ** Calculates clear sky radiation ** //
double clearskyradCpp(int jd, double lt, double lat, double lon,
    double z = 0, double TL = 1.0)
{
    // compute extra terrestrial rad
    double Io = extrarad(jd);
    // solar zenith
    double st = soltimeCpp(jd, lt, lon);
    double Z = solzenCpp(jd, st, lat, lon);
    // compute Linke turbidity
    double Icb = 0.0;
    if (Z <= M_PI / 2) {
        double Trb = dirtran(Z, TL, z);
        Icb = Io * cos(Z) * Trb;
    }
    double Trd = diftran(Z, TL);
    double Icd = Io * Trd;
    double Ic = Icb + Icd;
    return Ic;
}
// ** Calculates clear sky radiation - 6 hour averages starting at the reference time //
double clearskyradsixhr(int jd, double hr, double lat, double lon,
    double z = 0, double TL = 1.0)
{
    double csr = 0.0;
    for (int i = 0; i < 360; ++i) {
        double lt = hr + i / 60;
        csr = csr + clearskyradCpp(jd, lt, lat, lon, z, TL);
    }
    csr = csr / 360;
    return csr;
}
// Calculate Linke turbidity as a vector
// [[Rcpp::export]]
NumericVector Linkev(NumericVector sh, NumericVector pk, NumericVector aod)
{
    NumericVector TL(sh.size());
    for (R_xlen_t i = 0; i < sh.size(); ++i) {
        TL[i] = Linke(sh[i], pk[i], aod[i]);
    }
    return TL;
}
// ** Above for vector
// [[Rcpp::export]]
NumericVector clearskyradv(IntegerVector year, IntegerVector month, 
    IntegerVector day, NumericVector hr, double lat, double lon, double z,
    NumericVector TL)
{
    NumericVector csr(hr.size());
    for (R_xlen_t i = 0; i < hr.size(); ++i) {
        int jd = juldayCpp(year[i], month[i], day[i]);
        csr[i] = clearskyradsixhr(jd, hr[i], lat, lon, z, TL[i]);
    }
    return csr;
}
// [[Rcpp::export]]
NumericVector clearskyradhourly(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector hr, double lat, double lon, double z,
    NumericVector TL)
{
    NumericVector csr(hr.size());
    for (R_xlen_t i = 0; i < hr.size(); ++i) {
        int jd = juldayCpp(year[i], month[i], day[i]);
        csr[i] = clearskyradCpp(jd, hr[i], lat, lon, z, TL[i]);
    }
    return csr;
}
// ** Calculates diffuse fraction ** //
double difpropCpp(double swrad, int jd, double lt, double lat, double lon)
{
    double d = 1.0;
    if (swrad > 0) {
        double st = soltimeCpp(jd, lt, lon);
        double z = solzenCpp(jd, st, lat, lon);
        if (z < M_PI / 2) {
            double zd = z * 180 / M_PI;
            double k1 = 0.83 - 0.56 * exp(-0.06 * (90 - zd));
            double si = 0.0;
            if (z <= M_PI / 2) si = cos(z);
            double k = 0.0;
            if (si > 0) k = swrad / (1352 * si);
            if (k > k1) k = k1;
            if (k < 0) k = 0;
            double rho = k / k1;
            double sigma3 = 0;
            if (rho > 1.04) {
                sigma3 = 0.12 + 0.65 * (rho - 1.04);
            }
            else {
                sigma3 = 0.021 + 0.397 * rho - 0.231 * pow(rho, 2) - 0.13 * exp(-1 * pow((rho - 0.931) / 0.134, 2) * 0.834);
            }
            double k2 = 0.95 * k1;
            double d1 = 1.0;
            if (zd < 88.6) d1 = 0.07 + 0.046 * zd / (93 - zd);
            double K = 0.5 * (1 + sin(M_PI * (k - 0.22) / (k1 - 0.22) - M_PI / 2));
            double d2 = 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K * K));
            double d3 = (d2 * k2) * (1 - k) / (k * (1 - k2));
            double alpha = pow(1 / cos(z), 0.6);
            double kbmax = pow(0.81, alpha);
            double kmax = (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2));
            double dmax = (d2 * k2) * (1 - kmax) / (kmax * (1 - k2));
            d = 1 - kmax * (1 - dmax) / k;
            if (k <= kmax) d = d3;
            if (k <= k2) d = d2;
            if (k <= 0.22) d = 1;
            double kX = 0.56 - 0.32 * exp(-0.06 * (90 - zd));
            double kL = (k - 0.14) / (kX - 0.14);
            double kR = (k - kX) / 0.71;
            double delta = (k >= 0.14 && k < kX) ? (-3 * pow(kL, 2) * (1 - kL) * pow(sigma3, 1.3)) : 0;
            if (k >= kX && k < (kX + 0.71)) delta = 3 * kR * pow((1 - kR), 2) * pow(sigma3, 0.6);
            if (sigma3 > 0.01) d = d + delta;
        }
    }
    return d;
}
// ** Calculates diffuse fraction on a vector of data ** //
// [[Rcpp::export]]
NumericVector difpropvCpp(NumericVector swrad, IntegerVector year, 
    IntegerVector month, IntegerVector day, NumericVector hr, double lat, 
    double lon)
{
    NumericVector d(hr.size());
    for (R_xlen_t i = 0; i < hr.size(); ++i) {
        int jd = juldayCpp(year[i], month[i], day[i]);
        d[i] = difpropCpp(swrad[i], jd, hr[i], lat, lon);
    }
    return d;
}
// ** Compute six hour average temperature starting at reference time
// [[Rcpp::export]]
NumericVector tempsix(NumericVector tc) {
    int sh = tc.size() / 6;
    NumericVector tc6(sh);
    int index = 0;
    for (int i = 0; i < sh; ++i) {
        double st = 0.0;
        for (int j = 0; j < 6; ++j) {
            st = st + tc[index];
            ++index;
        }
        tc6[i] = st / 6.0;
    }
    return tc6;
}
// [[Rcpp::export]]
NumericVector prectohour(NumericVector prec) {
    int n = prec.size() * 6;
    NumericVector prech(n);
    int index = 0;
    for (R_xlen_t i = 0; i < prec.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            prech[index] = prec[i] * 3600.0;
            ++index;
        }
    }
    return prech;
}
// ** Blend met office and e.g. era5 climate data
// [[Rcpp::export]]
NumericVector blendtempCpp(NumericVector tasmin, NumericVector tasmax,
    NumericVector temp)
{
    // Do dimension and length checks
    IntegerVector dims1 = tasmin.attr("dim");
    IntegerVector dims2 = tasmax.attr("dim");
    IntegerVector dims3 = temp.attr("dim");
    for (int i = 0; i < 3; ++i) {
        if (dims1[i] != dims2[i]) stop("Dimensions of tasmin and tasmax don't match");
    }
    for (int i = 0; i < 2; ++i) {
        if (dims1[i] != dims3[i]) stop("xy dimensions of tasmin and temp don't match");
    }
    if (dims1[2] * 24 != dims3[2]) stop("temp must have 24 times more entries than tasmin");
    // Get dimensions
    int rows = dims1[0];
    int cols = dims1[1];
    int days = dims1[2];
    int hours = days * 24;
    NumericVector tcout(rows * cols * hours);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!NumericVector::is_na(temp[i + rows * (j + cols * 0)])) {
                // Extract vector of hourly temperatures
                NumericVector tempv(hours);
                for (int k = 0; k < hours; ++k) {
                    tempv[k] = temp[i + rows * (j + cols * k)];
                }
                // Extract vector of daily mins and maxes temperatures
                NumericVector tmxd(days);
                NumericVector tmnd(days);
                for (int k = 0; k < days; ++k) {
                    tmnd[k] = tasmin[i + rows * (j + cols * k)];
                    tmxd[k] = tasmax[i + rows * (j + cols * k)];
                }
                int index1 = 0;
                int index2 = 0;
                for (int day = 0; day < days; ++day) {
                    // For each day, calculate fraction of dtr in temp
                    double tmx = -273.15;
                    double tmn = 273.15;
                    for (int hr = 0; hr < 24; ++hr) {
                        if (tempv[index1] < tmn) tmn = tempv[index1];
                        if (tempv[index1] > tmx) tmx = tempv[index1];
                        ++index1;
                    }
                    double dtr1 = tmxd[day] - tmnd[day]; // dtr of met office
                    double dtr2 = tmx - tmn; // dtr of era5
                    for (int hr = 0; hr < 24; ++hr) {
                        double tfrac = (tempv[index2] - tmn) / dtr2;
                        double tc = tfrac * dtr1 + tmnd[day];
                        tcout[i + rows * (j + cols * index2)] = tc;
                        ++index2;
                    } // end hours
                } // end days
            } // NA check
            else {
                for (int k = 0; k < hours; ++k) {
                    tcout[i + rows * (j + cols * k)] = NA_REAL;
                } // end hours
            } // end NA
        } // end col
    } // end row
    tcout.attr("dim") = IntegerVector::create(rows, cols, hours);
    return tcout;
}
// Generic function to convert to daily
// [[Rcpp::export]]
NumericVector applytodaily(NumericVector hourly, std::string func)
{
    // Get dimensions of hourly
    IntegerVector dims = hourly.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hrs = dims[2];
    int days = hrs / 24.0;
    NumericVector daily(rows * cols * days);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hourly[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int day = 0; day < days; ++day) {
                    double x = 0.0;
                    if (func == "min") x = 99999999.0;
                    if (func == "max") x = -99999999.0;
                    for (int h = 0; h < 24; ++h) {
                        int idx = i + j * rows + (day * 24 + h) * rows * cols;
                        if (func == "mean") x += hourly[idx];
                        if (func == "sum") x += hourly[idx];
                        if (func == "min") {
                            if (hourly[idx] < x) x = hourly[idx];
                        }
                        if (func == "max") {
                            if (hourly[idx] > x) x = hourly[idx];
                        }
                    }
                    if (func == "mean") x = x / 24.0;
                    // Store the daily mean in the result vector
                    daily[i + j * rows + day * rows * cols] = x;
                } // end day
            } // end NA check
            else {
                for (int day = 0; day < days; ++day) {
                    daily[i + j * rows + day * rows * cols] = NA_REAL;
                } // end day
            } // end NA check
        } // end j
    } // end i
    daily.attr("dim") = IntegerVector::create(rows, cols, days);
    return daily;
}
// [[Rcpp::export]]
NumericVector blendprecCpp(NumericVector precd, NumericVector prech)
{
    // Do dimension and length checks
    IntegerVector dims1 = precd.attr("dim");
    // Get dimensions
    int rows = dims1[0];
    int cols = dims1[1];
    int days = dims1[2];
    int hrs = days * 24;
    // Calculate daily from hourly
    NumericVector precdh = applytodaily(prech, "sum");
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = precd[i + rows * j];
            if (!NumericVector::is_na(val)) {
                for (int h = 0; h < hrs; ++h) {
                    int d = floor(h / 24.0);
                    int idxd = i + rows * (j + cols * d);
                    int idxh = i + rows * (j + cols * h);
                    double mu = 1.0;
                    if (precdh[idxd] > 0.0) {
                        mu = precd[idxd] / precdh[idxd];
                    }
                    prech[idxh] = prech[idxh] * mu;
                    if (prech[idxh] < 0.0) prech[idxh] = 0.0;
                }
            }
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idxh = i + rows * (j + cols * h);
                    prech[idxh] = NA_REAL;
                } // end hours
            } // end NA check
        } // end j
    }  // end i
    return prech;
}


// [[Rcpp::export]]
NumericVector satvapCpp(NumericVector tc) {
    IntegerVector dims = tc.attr("dim");
    // Get dimensions
    int rows = dims[0];
    int cols = dims[1];
    int hrs = dims[2];
    int n = rows * cols * hrs;
    NumericVector es(n);
    for (int i = 0; i < n; ++i) {
        if (!NumericVector::is_na(tc[i])) {
            if (tc[i] > 0) {
                es[i] = 0.61078 * exp(17.27 * tc[i] / (tc[i] + 237.3));
            }
            else {
                es[i] = 0.61078 * exp(21.875 * tc[i] / (tc[i] + 265.5));
            }
        }
        else {
            es[i] = NA_REAL;
        }
    }
    es.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return es;
}
// [[Rcpp::export]]
NumericVector clearskya(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector hr, NumericMatrix lats, NumericMatrix lons)
{
    // Get dimensions
    IntegerVector dims = lats.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hrs = year.size();
    int n = rows * cols * hrs;
    NumericVector csr(n);
    IntegerVector jd(hrs);
    for (int k = 0; k < hrs; ++k) {
        jd[k] = juldayCpp(year[k], month[k], day[k]);
    }
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < hrs; ++k) {
                csr[index] = clearskyradCpp(jd[k], hr[k], lats(i, j),
                    lons(i, j), 0.0, 1.0);
                ++index;
            }
        }
    }
    csr = aperm3D(csr, rows, cols, hrs);
    csr.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return csr;
}
// ~~~~~~~~~~~~~ Coastal fill functions ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// get distance weighted value of all non NA cells in a matrix
double distwgt(NumericMatrix m, int ii, int jj, double pw = 2.0) 
{
    int rows = m.nrow();
    int cols = m.ncol();
    double swgt = 0.0;
    double val = 0.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = m(i, j);
            if (!R_IsNA(value) && !R_IsNaN(value)) {
                int r_dist = ii - i;
                int c_dist = jj - j;
                double dist = static_cast<double>(sqrt(pow(r_dist, 2.0) + pow(c_dist, 2.0)));
                double wgt = 1 / pow(dist, pw);
                swgt = swgt + wgt;
                val = val + wgt * value;
            } // end NA check
        } // end j
    } // end i
    val = val / swgt;
    return val;
}
// NA fill matrix by inverse distance
// [[Rcpp::export]]
NumericMatrix na_fill_matrix(NumericMatrix m, double pw = 2.0) 
{
    int rows = m.nrow();
    int cols = m.ncol();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = m(i, j);
            if (R_IsNA(value) || R_IsNaN(value)) {
                m(i, j) = distwgt(m, i, j, pw);
            }
        }
    }
    return m;
}
// Define a struct to hold cell information
struct Cell {
    int row;
    int col;
    double dist;
};
// Compare function for sorting cells by distance
bool compareDistance(const Cell& a, const Cell& b) {
    return a.dist < b.dist;
}
// [[Rcpp::export]]
List ordered_cells(NumericMatrix mat, int i, int j, double pw) {
    // Create a vector to store valid cells (non-NA, non-NaN)
    std::vector<Cell> valid_cells;
    // Get matrix dimensions
    int nrow = mat.nrow();
    int ncol = mat.ncol();
    // Iterate over all cells in the matrix
    for (int r = 0; r < nrow; r++) {
        for (int c = 0; c < ncol; c++) {
            // Check if the cell is not NA and not NaN
            double value = mat(r, c);
            if (!R_IsNA(value) && !R_IsNaN(value)) {
                // Calculate the Euclidean distance from the focal cell (i, j)
                double dist = std::sqrt(std::pow(r - i, 2) + std::pow(c - j, 2));
                // Store the cell coordinates and distance
                valid_cells.push_back({ r, c, dist });
            }
        }
    }
    // Sort valid cells by distance
    std::sort(valid_cells.begin(), valid_cells.end(), compareDistance);
    // Limit to the first 10 cells if there are more
    size_t n = valid_cells.size();
    if (n > 10) n = 10;
    // Create vectors to store the sorted row and column indices
    IntegerVector sorted_rows(n);
    IntegerVector sorted_cols(n);
    // Extract the sorted row and column indices
    for (size_t k = 0; k < n; k++) {
        sorted_rows[k] = valid_cells[k].row + 1; // R is 1-based
        sorted_cols[k] = valid_cells[k].col + 1; // R is 1-based
    }
    // Get the weights based on the distances
    NumericVector wgts(n);
    double swgt = 0.0;
    for (size_t k = 0; k < n; ++k) {
        double r_dist = sorted_rows[k] - i - 1; // Convert back to 0-based
        double c_dist = sorted_cols[k] - j - 1; // Convert back to 0-based
        double dist = std::sqrt(std::pow(r_dist, 2) + std::pow(c_dist, 2));
        if (dist > 0) { // Avoid division by zero
            wgts[k] = 1.0 / std::pow(dist, pw);
            swgt += wgts[k];
        }
        else {
            wgts[k] = 0.0;
        }
    }
    // Normalize the weights
    if (swgt > 0) {
        for (size_t k = 0; k < n; ++k) {
            wgts[k] /= swgt;
        }
    }
    // Substract one from sorted_rows and cols
    for (size_t k = 0; k < n; ++k) {
        sorted_rows[k] = sorted_rows[k] - 1;
        sorted_cols[k] = sorted_cols[k] - 1;
    }
    // Return a list of sorted row and column indices and the weights
    return List::create(Named("rows") = sorted_rows,
        Named("cols") = sorted_cols,
        Named("wgts") = wgts);
}
// [[Rcpp::export]]
NumericMatrix slicearray(NumericVector a, int d1, int d2, int d3, int k) {
    if (k < 0 || k >= d3) {
        stop("Invalid slice index.");
    }
    // Create a 2D matrix to store the slice
    NumericMatrix matrix_slice(d1, d2);
    // Loop over the rows and columns to extract the slice at index k
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            // Calculate the index in the 1D NumericVector
            int idx = i + d1 * (j + d2 * k);
            matrix_slice(i, j) = a[idx];
        }
    }
    return matrix_slice;
}
// Fill NAs in an array with inverse distance-weighted value of closests cells
// [[Rcpp::export]]
NumericVector na_fill_array(NumericVector clima, NumericMatrix lsea,
    double pw = 2.0)
{
    // get dims of clima
    IntegerVector dims = clima.attr("dim");
    int rws = dims[0];
    int cls = dims[1];
    int tsteps = dims[2];
    // get dims of lsea
    int nrow = lsea.nrow();
    int ncol = lsea.ncol();
    if (rws != nrow) stop("dims of clima and lsea must match");
    if (cls != ncol) stop("dims of clima and lsea must match");
    // slice out first slice of 3d array
    NumericMatrix m = slicearray(clima, rws, cls, tsteps, 0);
    for (int i = 0; i < rws; ++i) {
        for (int j = 0; j < cls; ++j) {
            // work out whether to bother doing anything
            double v1 = m(i, j);
            double v2 = lsea(i, j);
            if (R_IsNaN(v1)) v1 = NA_REAL;
            if (R_IsNaN(v2)) v2 = NA_REAL;
            if (!R_IsNA(v2) && R_IsNA(v1)) {
                // get the 10 closest cells and their weights
                List oc = ordered_cells(m, i, j, pw);
                NumericVector iiv = as<NumericVector>(oc["rows"]);
                NumericVector jjv = as<NumericVector>(oc["cols"]);
                NumericVector wgts = as<NumericVector>(oc["wgts"]);
                size_t n = iiv.size();
                // Loop through each k
                for (int k = 0; k < tsteps; ++k) {
                    // Get weighted vals of 10 closest non NA cells
                    double val = 0.0;
                    for (size_t h = 0; h < n; ++h) {
                        int idx = iiv[h] + rws * (jjv[h] + cls * k);
                        val = val + wgts[h] * clima[idx];
                    }
                    int idx2 = i + rws * (j + cls * k);
                    clima[idx2] = val;
                } // end k
            } // end check of lsea is not null
        } // end j
    } // end i
    return clima;
}
// [[Rcpp::export]]
NumericVector expandssttohour(NumericVector sst)
{
    // get dims of sst
    IntegerVector dims = sst.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // create sst hourly output
    NumericVector ssth(rows * cols * tsteps * 24);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                int idx = i + rows * (j + cols * k);
                double sstv = sst[idx];
                for (int h = 0; h < 24; ++h) {
                    int idx2 = i + rows * (j + cols * (k + 24 * h));
                    ssth[idx2] = sstv;
                }
            }
        }
    }
    ssth.attr("dim") = IntegerVector::create(rows, cols, tsteps * 24);
    return ssth;
}
// [[Rcpp::export]]
NumericVector remove_nas(NumericVector v) {
    // Create a NumericVector to store non-NA and non-NaN values
    NumericVector non_na_nan;
    // Iterate through the input vector
    for (R_xlen_t i = 0; i < v.size(); i++) {
        // Check if the element is neither NA nor NaN
        if (!R_IsNA(v[i]) && !R_IsNaN(v[i])) {
            non_na_nan.push_back(v[i]); // Add the valid value to the new vector
        }
    }
    return non_na_nan; // Return the vector without NAs and NaNs
}
// [[Rcpp::export]]
double calculate_mode(NumericVector x) {
    std::map<double, int> frequency; // Map to store the frequency of each value
    // Calculate the frequency of each value in the NumericVector
    for (R_xlen_t i = 0; i < x.size(); i++) {
        frequency[x[i]]++;
    }
    // Find the mode
    double mode_value = R_NaReal; // Initialize to NA
    int max_count = 0;
    for (const auto& pair : frequency) {
        if (pair.second > max_count) {
            max_count = pair.second;
            mode_value = pair.first;
        }
    }
    // Handle the case where no mode exists (all values are unique)
    if (max_count <= 1) {
        return R_NaReal; // If no mode, return NA
    }
    return mode_value; // Return the mode value as a double
}
// [[Rcpp::export]]
NumericVector modalwinddir(NumericVector wd)
{
    // get dims of wind direction
    IntegerVector dims = wd.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // Round wind direction
    int n = wd.size();
    NumericVector rwd(n);
    for (int i = 0; i < n; i++) rwd[i] = round(wd[i] / 5.0) * 5.0;
    rwd.attr("dim") = IntegerVector::create(rows, cols, tsteps);
    NumericVector wmode(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        NumericVector wdv(rows * cols);
        int index = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx = i + rows * (j + cols * k);
                wdv[index] = rwd[idx];
                ++index;
            } // end j
        } // end i
        wdv = remove_nas(wdv);
        wmode[k] = calculate_mode(wdv);
    } // end k
    return wmode;
}
// [[Rcpp::export]]
NumericVector mattoarray(NumericMatrix m, int n)
{
    // get dims of lsea
    int nrow = m.nrow();
    int ncol = m.ncol();
    NumericVector a(nrow * ncol * n);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            for (int k = 0; k < n; ++k) {
                int idx = i + nrow * (j + ncol * k);
                a[idx] = m(i, j);
            }
        }
    }
    a.attr("dim") = IntegerVector::create(nrow, ncol, n);
    return a;
}
// [[Rcpp::export]]
NumericVector solarindexarray(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector hour, NumericMatrix lats, NumericMatrix lons)
{
    // Get dims
    int nrow = lats.nrow();
    int ncol = lats.ncol();
    int n = year.size();
    // Calculate julian day vector
    NumericVector jd(n);
    for (int k = 0; k < n; ++k) jd[k] = juldayCpp(year[k], month[k], day[k]);
    // Get solar index
    NumericVector si(nrow * ncol * n);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            for (int k = 0; k < n; ++k) {
                int idx = i + nrow * (j + ncol * k);
                double st = soltimeCpp(jd[k], hour[k], lons(i, j));
                double z = solzenCpp(jd[k], st, lats(i, j), lons(i, j));
                si[idx] = cos(z);
                if (si[idx] < 0.0) si[idx] = 0.0;
            }
        }
    }
    si.attr("dim") = IntegerVector::create(nrow, ncol, n);
    return si;
}
// [[Rcpp::export]]
NumericMatrix applymeanCpp(NumericVector a)
{
    // get dims of wind direction
    IntegerVector dims = a.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int n = dims[2];
    NumericMatrix m(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = a[i + rows * j];
            if (!Rcpp::NumericVector::is_na(val)) {
                double sm = 0.0;
                for (int k = 0; k < n; ++k) {
                    int idx = i + rows * (j + cols * k);
                    sm = sm + a[idx];
                }
                m(i, j) = sm / n;
            }
            else {
                m(i, j) = NA_REAL;
            }
        }
    }
    return m;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********** Code for doing bias correct and temporal downscaling ********* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// *********************** Convert hourly to daily ************************* //
// hourly temperature to dtr and mean
// [[Rcpp::export]]
List temptoday(NumericVector tc) {
    // Calculate mean temperature
    NumericVector tmean = applytodaily(tc, "mean");
    NumericVector tmx = applytodaily(tc, "max");
    NumericVector tmn = applytodaily(tc, "min");
    NumericVector dtr(tmn.size());
    for (R_xlen_t i = 0; i < tmx.size(); ++i) {
        if (!NumericVector::is_na(tmx[i])) {
            dtr[i] = tmx[i] - tmn[i];
        }
        else {
            dtr[i] = NA_REAL;
        }
    }
    IntegerVector dims = tmn.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    dtr.attr("dim") = IntegerVector::create(rows, cols, days);
    List result = List::create(Named("tmean") = tmean,
        Named("dtr") = dtr);
    return result;
}
// hourly relative humidity to daily vapour pressure
// [[Rcpp::export]]
NumericVector rhtoday(NumericVector rh, NumericVector tc)
{
    NumericVector es = satvapCpp(tc);
    NumericVector ea(es.size());
    for (R_xlen_t i = 0; i < es.size(); ++i) {
        if (!NumericVector::is_na(es[i])) {
            ea[i] = es[i] * rh[i] / 100.0;
        }
        else {
            ea[i] = NA_REAL;
        }
    }
    IntegerVector dims = rh.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hours = dims[2];
    ea.attr("dim") = IntegerVector::create(rows, cols, hours);
    NumericVector ead = applytodaily(ea, "mean");
    return ead;
}
// hourly longwave radiation to daily sky emissivity
// [[Rcpp::export]]
NumericVector lwtoday(NumericVector lwdown, NumericVector tc)
{
    // Hourly sky emissivity
    NumericVector skyem(tc.size());
    for (R_xlen_t i = 0; i < tc.size(); ++i) {
        if (!NumericVector::is_na(tc[i])) {
            double lwup = 5.67 * pow(10.0, -8.0) * pow(tc[i] + 273.15, 4.0);
            skyem[i] = lwdown[i] / lwup;
            if (skyem[i] > 1.0) skyem[i] = 1.0;
            if (skyem[i] < 0.0) skyem[i] = 0.0;
        }
        else {
            skyem[i] = NA_REAL;
        }
    }
    IntegerVector dims = tc.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hours = dims[2];
    skyem.attr("dim") = IntegerVector::create(rows, cols, hours);
    NumericVector skyemd = applytodaily(skyem, "mean");
    return skyemd;
}
// *********************** Convert tasmin and tasmax ********************** //
// [[Rcpp::export]]
List tmnmxtodtr(NumericVector tmn, NumericVector tmx)
{
    // Get dimension
    IntegerVector dims = tmn.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create variables
    NumericVector tmean(rows * cols * days);
    NumericVector dtr(rows * cols * days);
    // Loop through pixels
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = tmn[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    tmean[idx] = (tmn[idx] + tmx[idx]) / 2.0;
                    dtr[idx] = tmx[idx] - tmn[idx];
                }
            }
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    tmean[idx] = NA_REAL;
                    dtr[idx] = NA_REAL;
                } // end day
            } // end NA check
        } // end j
    } // end i
    tmean.attr("dim") = IntegerVector::create(rows, cols, days);
    dtr.attr("dim") = IntegerVector::create(rows, cols, days);
    List result = List::create(Named("tmean") = tmean,
        Named("dtr") = dtr);
    return result;
}
// [[Rcpp::export]]
NumericVector sphtoea(NumericVector sph, NumericVector pk)
{
    // Get dimension
    IntegerVector dims = sph.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create ea vector
    NumericVector ea(rows * cols * days);
    // Loop through pixels
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = sph[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    ea[idx] = (sph[idx] * pk[idx]) /
                        (0.622 + (0.378 * sph[idx]));
                } // end d
            } // end NA check
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    ea[idx] = NA_REAL;
                } // end d
            } // end NA check
        } // end j
    } // end i
    ea.attr("dim") = IntegerVector::create(rows, cols, days);
    return ea;
}
// [[Rcpp::export]]
NumericVector netshorttodownshort(NumericVector rss, NumericVector alb)
{
    // Get dimension
    IntegerVector dims = rss.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create ea vector
    NumericVector rsw(rows * cols * days);
    // Loop through pixels
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = alb[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    rsw[idx] = rss[idx] / (1 - alb[idx]);
                } // end d
            } // end NA check
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    rsw[idx] = NA_REAL;
                } // end d
            } // end NA check
        } // end j
    } // end i
    rsw.attr("dim") = IntegerVector::create(rows, cols, days);
    return rsw;
}
// [[Rcpp::export]]
NumericVector netlwtoskyem(NumericVector rsl, NumericVector tmx, NumericVector tmn)
{
    // Get dimension
    IntegerVector dims = rsl.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create ea vector
    NumericVector skyem(rows * cols * days);
    // Loop through pixels
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = rsl[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    double tc = (tmx[idx] + tmn[idx]) / 2;
                    double lwup = 5.67 * pow(10, -8) * 
                        pow(tc + 273.15, 4);
                    double lwd = rsl[idx] + lwup;
                    skyem[idx] = lwd / lwup;
                } // end d
            } // end NA check
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    skyem[idx] = NA_REAL;
                } // end d
            } // end NA check
        } // end j
    } // end i
    skyem.attr("dim") = IntegerVector::create(rows, cols, days);
    return skyem;
}
// [[Rcpp::export]]
NumericVector uvtows(NumericVector u, NumericVector v)
{
    // Get dimension
    IntegerVector dims = u.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create ea vector
    NumericVector ws(rows * cols * days);
    // Loop through pixels
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = u[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    ws[idx] = sqrt(u[idx] * u[idx] + v[idx] * v[idx]);
                } // end d
            } // end NA check
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    ws[idx] = NA_REAL;
                } // end d
            } // end NA check
        } // end j
    } // end i
    ws.attr("dim") = IntegerVector::create(rows, cols, days);
    return ws;
}
// *********************** Correct rainfall ********************** //
std::vector<std::vector<double>> convertoCppmatrix(NumericMatrix mat) {
    std::vector<std::vector<double>> result(mat.nrow(), std::vector<double>(mat.ncol()));
    for (int i = 0; i < mat.nrow(); ++i) {
        for (int j = 0; j < mat.ncol(); ++j) {
            result[i][j] = mat(i, j);
        }
    }
    return result;
}
// ** Convert C++ matrix to R matrix ** //
NumericMatrix convertoRmatrix(std::vector<std::vector<double>>& mat) {
    int nrow = mat.size();
    int ncol = mat[0].size(); // Assuming all inner vectors have the same size
    NumericMatrix mat2(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            mat2(i, j) = mat[i][j];
        }
    }
    return mat2;
}
std::vector<double> rainadjustv(std::vector<double> rain, std::vector<double> rrain, double rfrac, double rtot)
{
    // Calculate expected rianfall days/hours
    if (!std::isnan(rain[0])) {
        double pr = rfrac * rain.size();
        int prd = std::round(pr);
        int ard = 0;
        // Calculate actual rainfall days/hours
        for (size_t i = 0; i < rain.size(); ++i) if (rain[i] > 0) ard++;
        // if actual rain days/hours < predicted rain days/hours give some random no rain days a small amount of rain
        if (ard < prd) {
            // Calculate number of zero days/hours to assign rain to
            int rta = prd - ard;
            // Calculate minimum none-zero in rain
            double nz = 5001;
            for (size_t i = 0; i < rain.size(); ++i) if (rain[i] > 0 && rain[i] < nz) nz = rain[i];
            // Create a vector of indices corresponding to zeros in rain
            std::vector<int> zeroIndices;
            for (size_t i = 0; i < rain.size(); ++i) {
                if (rain[i] == 0) {
                    zeroIndices.push_back(i);
                }
            }
            // Sort zeroIndices based on corresponding values in rrain
            std::sort(zeroIndices.begin(), zeroIndices.end(), [&rrain](int i, int j) { return rrain[i] < rrain[j]; });
            // Assign the lowest rta indices a value nz
            int limit = std::min(rta, static_cast<int>(zeroIndices.size())); // Cast size to int
            for (int i = 0; i < limit; ++i) {
                rain[zeroIndices[i]] = nz; // or any desired non-zero value
            }
        }
        // if actual rain days > predicted rain days give lowest rainfall days zero rain
        if (ard > prd) {
            // Find number of days/hours that should be zero
            int nrd = rain.size() - prd;
            // Create a vector of indices ordered by the values in rain
            std::vector<size_t> indices(rain.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&rain](size_t i, size_t j) { return rain[i] < rain[j]; });
            // Assign 0 to the lowest nrd elements
            int limit = std::min(nrd, static_cast<int>(rain.size())); // Cast size to int
            for (int i = 0; i < limit; ++i) {
                rain[indices[i]] = 0;
            }
        }
        // adjust rainfall to match total
        double rsum = rain[0];
        for (size_t i = 1; i < rain.size(); ++i) rsum = rsum + rain[i];
        double mu = rtot / rsum;
        for (size_t i = 0; i < rain.size(); ++i) rain[i] = rain[i] * mu;
    }
    return rain;
}
// [[Rcpp::export]]
NumericMatrix rainadjustm(NumericMatrix rainm, std::vector<double> rrain, std::vector<double> rfrac, std::vector<double> rtot)
{
    std::vector<std::vector<double>> rainc = convertoCppmatrix(rainm);
    for (size_t i = 0; i < rainc.size(); ++i) {
        rainc[i] = rainadjustv(rainc[i], rrain, rfrac[i], rtot[i]);
    }
    rainm = convertoRmatrix(rainc);
    return rainm;
}
// *********************** Apply range lims to gam ********************** //
// [[Rcpp::export]]
NumericVector rangelimapply(NumericVector v1, NumericVector v2,
    NumericVector v3, NumericVector v3c, double rangelim = 1.05)
{
    // Calculate means
    double v1m = 0.0;
    double v2m = 0.0;
    double v3m = 0.0;
    for (R_xlen_t i = 0; i < v1.size(); ++i) v1m += v1[i];
    for (R_xlen_t i = 0; i < v2.size(); ++i) v2m += v2[i];
    for (R_xlen_t i = 0; i < v3.size(); ++i) v3m += v3[i];
    v1m = v1m / v1.size();
    v2m = v2m / v2.size();
    v3m = v3m / v3.size();
    // Calculate expected mean
    double exp_mean = (v1m - v2m) + v3m;
    // Calculate max differences from means
    double v1d = 0.0;
    double v2d = 0.0;
    double v3d = 0.0;
    for (R_xlen_t i = 0; i < v1.size(); ++i) {
        double dif = abs(v1[i] - v1m);
        if (dif > v1d) v1d = dif;
    }
    for (R_xlen_t i = 0; i < v2.size(); ++i) {
        double dif = abs(v2[i] - v2m);
        if (dif > v2d) v2d = dif;
    }
    for (R_xlen_t i = 0; i < v3.size(); ++i) {
        double dif = abs(v3[i] - v3m);
        if (dif > v3d) v3d = dif;
    }
    // Calculate expected max difference
    double exp_maxdif = (v1d / v2d) * v3d * rangelim;
    // Calculate min and max
    double mn = exp_mean - exp_maxdif;
    double mx = exp_mean + exp_maxdif;
    // Apply limits
    for (R_xlen_t i = 0; i < v3.size(); ++i) {
        if (v3c[i] < mn) v3c[i] = mn;
        if (v3c[i] > mx) v3c[i] = mx;
    }
    return v3c;
}
// Apply range limit - precipitation
// [[Rcpp::export]]
NumericVector prangelimapply(NumericVector v1, NumericVector v2,
    NumericVector v3, NumericVector v3c, double rangelim = 1.05)
{
    // Calculate maximums
    double v1m = 0.0;
    double v2m = 0.0;
    double v3m = 0.0;
    for (R_xlen_t i = 0; i < v1.size(); ++i) if (v1[i] > v1m) v1m = v1[i];
    for (R_xlen_t i = 0; i < v2.size(); ++i) if (v2[i] > v2m) v2m = v1[i];
    for (R_xlen_t i = 0; i < v3.size(); ++i) if (v3[i] > v3m) v3m = v1[i];
    // Calculate maximum
    double mx = (v1m / v2m) * v3m;
    for (R_xlen_t i = 0; i < v3c.size(); ++i) if (v3c[i] > mx) v3c[i] = mx;
    return v3c;
}
// ************ Expand corrected daily data to hourly ********************** //
// Extract vector
NumericVector extractv(NumericVector a, int i, int j, int rows, int cols, int days) {
    NumericVector v(days);
    for (int t = 0; t < days; ++t) {
        // Calculate the index for (i, j, t) in the flattened array
        int idx = i + j * rows + t * (rows * cols);
        v[t] = a[idx];
    }
    return v;
}
// Temperature
// [[Rcpp::export]]
NumericVector tempha(NumericVector tmean, NumericVector dtr,
    IntegerVector year, IntegerVector month, IntegerVector day,
    NumericMatrix lats, NumericMatrix lons, double srte = 0.09)
{
    // Get dimensions of hourly
    IntegerVector dims = tmean.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    int hrs = days * 24.0;
    NumericVector temp(rows * cols * hrs);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = tmean[i + j * rows];
            if (!NumericVector::is_na(val)) {
                // Extract vectors and convert to tmn and max
                NumericVector tmev = extractv(tmean, i, j, rows, cols, days);
                NumericVector dtrv = extractv(dtr, i, j, rows, cols, days);
                NumericVector tmnv(tmev.size());
                NumericVector tmxv(tmev.size());
                for (R_xlen_t k = 0; k < tmev.size(); ++k) {
                    tmnv[k] = tmev[k] - 0.5 * dtrv[k];
                    tmxv[k] = tmev[k] + 0.5 * dtrv[k];
                }
                NumericVector tcv = hourlytempv(tmnv, tmxv, year, month,
                    day, lats(i, j), lons(i, j), srte);
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    temp[idx] = tcv[h];
                } // end h
            } // end NA check
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    temp[idx] = NA_REAL;
                } // end h
            } // end NA check
        } // end j
    } // end i
    temp.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return temp;
}
// Atmospheric Pressure
// [[Rcpp::export]]
NumericVector splina(NumericVector a)
{
    // Get dimensions of hourly
    IntegerVector dims = a.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    int hrs = days * 24.0;
    NumericVector ah(rows * cols * hrs);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = a[i + j * rows];
            if (!NumericVector::is_na(val)) {
                NumericVector av = extractv(a, i, j, rows, cols, days);
                NumericVector ahv = splineCpp(av);
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    ah[idx] = ahv[h];
                } // end h
            } // end NA check
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    ah[idx] = NA_REAL;
                } // end h
            } // end NA check
        } // end j
    } // end i
    ah.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return ah;
}
// Relative humidity
// [[Rcpp::export]]
NumericVector relhuma(NumericVector ea, NumericVector tch)
{
    NumericVector eah = splina(ea);
    // Get dimensions of hourly
    IntegerVector dims = eah.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hrs = dims[2];
    NumericVector rh(rows * cols * hrs);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = eah[i + j * rows];
            if (!NumericVector::is_na(val)) {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    double es = 0.61078 * exp(17.27 * tch[idx] / (tch[idx] + 237.3));
                    if (tch[idx] <= 0.0) {
                        es = 0.61078 * exp(21.875 * tch[idx] / (tch[idx] + 265.5));
                    }
                    rh[idx] = (eah[idx] / es) * 100.0;
                    if (rh[idx] > 100.0) rh[idx] = 100.0;
                    if (rh[idx] < 10.0) rh[idx] = 10.0;
                }
            }
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    rh[idx] = NA_REAL;
                } // end h
            } // end NA check
        } // end j
    }// end i
    rh.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return rh;
} 
// Shortwave Radiation
// ~~ Calculate daily clear-sky
// [[Rcpp::export]]
NumericMatrix dailyclm(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector lats)
{
    // Get dimensions
    int n = lats.size();
    int days = year.size();
    NumericMatrix csr(n, days);
    // Calculate julian day
    IntegerVector jd(days);
    for (int k = 0; k < days; ++k) {
        jd[k] = juldayCpp(year[k], month[k], day[k]);
    }
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < days; ++d) {
            double csrv = 0.0;
            for (int h = 0; h < 240; ++h) {
                csrv += clearskyradCpp(jd[d], h / 10, lats[i],
                    0.0, 0.0, 1.0);
            } // end 6 minute slice
            csr(i, d) = csrv / 240;
        } // end d
    } // end j
    return csr;
}
// [[Rcpp::export]]
NumericVector dailycla(NumericMatrix csr, int cols) {
    int rows = csr.nrow();  
    int n = csr.ncol(); 
    NumericVector csra(rows * cols * n);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < n; ++k) {
                int idx = i + j * rows + k * rows * cols;
                csra[idx] = csr(i, k);
            } // end k
        } // end j
    } // end i
    csra.attr("dim") = IntegerVector::create(rows, cols, n);
    return csra;
}
// ~~ Calculate daily clear-sky
// [[Rcpp::export]]
NumericMatrix hourlyclm(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericVector lats)
{
    // Get dimensions
    int n = lats.size();
    int days = year.size();
    int hrs = days * 24;
    NumericMatrix csr(n, hrs);
    // Calculate julian day
    IntegerVector jd(days);
    for (int k = 0; k < days; ++k) {
        jd[k] = juldayCpp(year[k], month[k], day[k]);
    }
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < days; ++d) {
            for (int h = 0; h < 24; ++h) {
                int j = d * 24 + h;
                double csrv = 0.0;
                for (int m = 0; m < 10; ++m) {
                    double hh = h + m / 10;
                    csrv += clearskyradCpp(jd[d], hh, lats[i],
                        0.0, 0.0, 1.0);
                } // end 6 minute slice
                csr(i, j) = csrv / 10;
            }
            
        } // end d
    } // end j
    return csr;
}
// [[Rcpp::export]]
NumericVector hourlycla(IntegerVector year, IntegerVector month,
    IntegerVector day, NumericMatrix lats, NumericMatrix lons)
{
    // Get dimensions
    int rows = lats.nrow();
    int cols = lats.ncol();
    int days = year.size();
    int hrs = days * 24;
    NumericVector csrh(rows * cols * hrs);
    // Calculate julian day
    IntegerVector jd(days);
    for (int k = 0; k < days; ++k) {
        jd[k] = juldayCpp(year[k], month[k], day[k]);
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int h = 0; h < hrs; ++h) {
                int d = static_cast<int>(std::floor(h));
                double csrv = 0.0;
                for (int m = 0; m < 10; ++m) {
                    double hh = h % 24  + m / 10;
                    csrv += clearskyradCpp(jd[d], hh, lats(i, j),
                        lons(i, j), 0.0, 1.0);
                }
                int idx = i + j * rows + h * rows * cols;
                csrh[idx] = csrv / 10;
            } // end h
        } // end j
    } // end i
    csrh.attr("dim") = NumericVector::create(rows, cols, hrs);
    return csrh;
}
// Diffuse Radiation
// [[Rcpp::export]]
NumericVector hourlydifa(NumericVector swrad, IntegerVector year, IntegerVector month,
    IntegerVector day, NumericMatrix lats, NumericMatrix lons)
{
    // Get dimensions
    int rows = lats.nrow();
    int cols = lats.ncol();
    int days = year.size();
    int hrs = days * 24;
    NumericVector difrh(rows * cols * hrs);
    // Calculate julian day
    IntegerVector jd(days);
    for (int k = 0; k < days; ++k) {
        jd[k] = juldayCpp(year[k], month[k], day[k]);
    }
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int val = swrad[i + j * rows];
            if (!Rcpp::NumericVector::is_na(val)) {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    if (swrad[idx] > 0.0) {
                        int d = static_cast<int>(std::floor(h));
                        double hh = h % 24;
                        double dp = difpropCpp(swrad[idx], jd[d], hh,
                            lats(i, j), lons(i, j));
                        difrh[idx] = swrad[idx] * dp;
                    }
                    else difrh[idx] = 0.0;
                } // end h
            } // end NA check
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    difrh[idx] = NA_REAL;
                } // end h
            } // end NA check
        } // end j
    } // end i
    difrh.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return difrh;
}
// Longwave Radiation
// [[Rcpp::export]]
NumericVector lwrada(NumericVector skyem, NumericVector tch)
{
    NumericVector skyemh = splina(skyem);
    // Get dimensions of hourly
    IntegerVector dims = skyemh.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int hrs = dims[2];
    NumericVector lwdown(rows * cols * hrs);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int val = lwdown[i + j * rows];
            if (!Rcpp::NumericVector::is_na(val)) {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    double lwup = 5.67 * pow(10, -8) * pow(tch[idx] + 273.15, 4);
                    lwdown[idx] = skyemh[idx] * lwup;
                } // end h
            } // end NA check
            else {
                for (int h = 0; h < hrs; ++h) {
                    int idx = i + j * rows + h * rows * cols;
                    lwdown[idx] = NA_REAL;
                } // end h
            } // end NA check
        } // end j
    } // end i
    // interpolate sky em to hourly
    lwdown.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return lwdown;
}
// Precipitation
// [[Rcpp::export]]
NumericVector precha(NumericVector precd)
{
    IntegerVector dims = precd.attr("dim");
    // Extract rows, columns, and layers
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    int hrs = days * 24;
    // spline interpolate rainfall
    NumericVector prech = splina(precd);
    prech.attr("dim") = IntegerVector::create(rows, cols, hrs);
    // Calculate and adjust daily total of hourly
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int val = precd[i + j * rows];
            if (!Rcpp::NumericVector::is_na(val)) {
                int hh = 0;
                int h2 = 0;
                for (int d = 0; d < days; ++d) {
                    // Calculate daily total of splined hourly rainfall
                    double smr = 0.0;
                    for (int h = 0; h < 24; ++h) {
                        int idxh = i + j * rows + hh * rows * cols;
                        smr += (prech[idxh] / 24.0);
                        ++hh;
                    }
                    // Calculate and apply adjustment to ensure consistency 
                    int idxd = i + j * rows + d * rows * cols;
                    double mu = precd[idxd] / smr;
                    for (int h = 0; h < 24; ++h) {
                        int idxh = i + j * rows + h2 * rows * cols;
                        prech[idxh] = mu * prech[idxh];
                        // Do NA and negative check
                        if (prech[idxh] <= 0.0) prech[idxh] = precd[idxd] / 24.0;
                        ++h2;
                    } // end hour
                } // end day
            } // end NA
        } // end j
    } // end i
    prech.attr("dim") = IntegerVector::create(rows, cols, hrs);
    return prech;
}
// ~~~~~~~~~~~~~  New wind direction coastal exposure function ~~~~~~~~~~~~ //
// [[Rcpp::export]]
NumericVector coastexpm(NumericMatrix coastm, int n)
{
    int rows = coastm.nrow(); 
    int cols = coastm.ncol(); 
    NumericVector coastma(rows * cols * n);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int val = coastm(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                for (int k = 0; k < n; ++k) {
                    int idx = i + j * rows + k * rows * cols;
                    coastma[idx] = coastm(i, j);
                } // end k
            } // end NA check
            else {
                for (int k = 0; k < n; ++k) {
                    int idx = i + j * rows + k * rows * cols;
                    coastma[idx] = NA_REAL;
                } // end k
            } // end NA check
        } // end j
    } // end i
    coastma.attr("dim") = IntegerVector::create(rows, cols, n);
    return coastma;
}
// [[Rcpp::export]]
NumericVector coastexpa(NumericVector winddir, NumericVector coastexp, int n = 8)
{
    // get dims
    IntegerVector dims = winddir.attr("dim");
    // Extract rows, columns, and layers
    int rows = dims[0];
    int cols = dims[1];
    int days = dims[2];
    // Create output vector
    NumericVector coasta(rows * cols * days);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int val = winddir[i + j * rows];
            if (!Rcpp::NumericVector::is_na(val)) {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    int i8 = winddir[idx] / (360 / n);
                    int idx8 = i + j * rows + i8 * rows * cols;
                    coasta[idx] = coastexp[idx8];
                } // end d
            } // end NA check
            else {
                for (int d = 0; d < days; ++d) {
                    int idx = i + j * rows + d * rows * cols;
                    coasta[idx] = NA_REAL;
                } // end d
            } // end NA check
        } // end j
    } // end i
    coasta.attr("dim") = IntegerVector::create(rows, cols, days);
    return coasta;
}
// [[Rcpp::export]]
NumericMatrix apply3D(NumericVector a) {
    // Get the dimensions of the 3D array
    IntegerVector dims = a.attr("dim");
    int nrows = dims[0];
    int ncols = dims[1];
    int ntime = dims[2];
    // Initialize the result matrix
    NumericMatrix result(nrows, ncols);
    // Iterate through each pixel
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            double sum = 0.0;
            int count = 0;
            // Iterate through the time dimension
            for (int t = 0; t < ntime; ++t) {
                double value = a[i + nrows * (j + ncols * t)];
                if (!NumericVector::is_na(value)) { // Skip NA values
                    sum += value;
                    count++;
                }
            }
            // Compute the mean, assign NA if no valid values
            if (count > 0) {
                result(i, j) = sum / count;
            }
            else {
                result(i, j) = NA_REAL;
            }
        }
    }
    return result;
}
// ====================================================================== //
// ~ Apply two-stream model: used for relating NDVI to reflectance etc. ~ //
// ====================================================================== //
// Function used to solve ndvi for a given pai
double ndvicpp(double pai, double ndviin)
{
    // Red
    double S1 = exp(-0.956782804 * pai);
    double D1 = -0.682265131 / S1 + 0.00796106 * S1;
    double p1 = (0.061666667 / (D1 * S1)) * -0.357497089;
    double p2 = (-0.061666667 * S1 / D1) * 1.556068518;
    double alb_red = p1 + p2;
    // NIR
    double alb_nir = 0.212278228;
    // calculate NDVI
    double ndvi = (alb_nir - alb_red) / (alb_nir + alb_red);
    // Output
    double out = ndvi - ndviin;
    return out;
}
// Root-finding function using the bisection method: ndvi
double solve_ndvi(double ndviin, double tol = 1e-6, int max_iter = 100) {
    double lower = 0.0;
    double upper = 50.0;
    double mid = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        mid = (lower + upper) / 2.0;
        double f_lower = ndvicpp(lower, ndviin);
        double f_mid = ndvicpp(mid, ndviin);
        if (std::abs(f_mid) < tol) {
            return mid; // Root found
        }
        if (f_lower * f_mid < 0) {
            upper = mid; // Root lies in the lower half
        }
        else {
            lower = mid; // Root lies in the upper half
        }
    }
    mid = NA_REAL;
    return mid; // Should never reach here
}
// [[Rcpp::export]]
NumericMatrix find_pai(NumericMatrix ndvi)
{
    // Get dimensions
    int rows = ndvi.nrow();
    int cols = ndvi.ncol();
    NumericMatrix pai(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = ndvi(i, j);
            if (!NumericMatrix::is_na(val)) {
                pai(i, j) = solve_ndvi(val);
            }
            else {
                pai(i, j) = NA_REAL;
            }
        }
    }
    return pai;
}
// Function used to calculate leaf or ground reflectance solve
// [[Rcpp::export]]
double leafrcpp(double lref, double pai, double gref, double x, double albin)
{
    // Base parameters
    double ltra = 0.33 * lref;
    double om = lref + ltra;
    double a = 1 - om;
    double del = lref - ltra;
    double J = 1.0 / 3.0;
    if (x != 1.0) {
        double mla = 9.65 * pow((3 + x), -1.65);
        if (mla > M_PI / 2) mla = M_PI / 2;
        J = cos(mla) * cos(mla);
    }
    double gma = 0.5 * (om + J * del);
    double h = sqrt(a * a + 2 * a * gma);
    // Calculate base parameters: diffuse
    double S1 = exp(-h * pai);
    double u1 = a + gma * (1 - 1 / gref);
    double D1 = (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h)
        * (u1 + h) * S1;
    // Calculate parameters: diffuse
    double p1 = (gma / (D1 * S1)) * (u1 - h);
    double p2 = (-gma * S1 / D1) * (u1 + h);
    double albd = p1 + p2;
    // Output
    double out = albd - albin;
    return out;
}
// Root-finding function using the bisection method: leafr
double solve_lref(double pai, double gref, double x, double albin, double tol = 1e-6, int max_iter = 100) {
    double lower = 0.0001;
    double upper = 0.6665;
    // check whether between lower and mid or mid and upper
    double mid = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        mid = (lower + upper) / 2.0;
        double f_lower = leafrcpp(lower, pai, gref, x, albin);
        double f_mid = leafrcpp(mid, pai, gref, x, albin);
        if (std::abs(f_mid) < tol) {
            return mid; // Root found
        }
        if (f_lower * f_mid < 0) {
            upper = mid; // Root lies in the lower half
        }
        else {
            lower = mid; // Root lies in the upper half
        }
    }
    return mid; // Should never reach here
}
// Root-finding function using the bisection method: gref
// [[Rcpp::export]]
double solve_gref(double lref, double pai, double x, double albin, double tol = 1e-6, int max_iter = 100) {
    double lower = 0.0001;
    double upper = 0.9999;
    // check whether between lower and mid or mid and upper
    double mid = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        mid = (lower + upper) / 2.0;
        double f_lower = leafrcpp(lref, pai, lower, x, albin);
        double f_mid = leafrcpp(lref, pai, mid, x, albin);
        if (std::abs(f_mid) < tol) {
            return mid; // Root found
        }
        if (f_lower * f_mid < 0) {
            upper = mid; // Root lies in the lower half
        }
        else {
            lower = mid; // Root lies in the upper half
        }
    }
    return mid; // Should never reach here
}
// [[Rcpp::export]]
NumericMatrix find_gref(NumericMatrix lref, NumericMatrix pai, 
    NumericMatrix x, NumericMatrix albin) 
{
    // Get dimensions
    int rows = pai.nrow();
    int cols = pai.ncol();
    NumericMatrix gref(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = pai(i, j);
            if (NumericMatrix::is_na(lref(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(x(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(albin(i, j))) val = NA_REAL;
            if (!NumericMatrix::is_na(val)) {
                gref(i, j) = solve_gref(lref(i, j), val, x(i, j), albin(i, j));
            }
            else {
                gref(i, j) = NA_REAL;
            }
        }
    }
    return gref;
}
// [[Rcpp::export]]
NumericMatrix find_lref(NumericMatrix pai, NumericMatrix gref,
    NumericMatrix x, NumericMatrix albin)
{
    // Get dimensions
    int rows = pai.nrow();
    int cols = pai.ncol();
    NumericMatrix lref(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = pai(i, j);
            if (NumericMatrix::is_na(gref(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(x(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(albin(i, j))) val = NA_REAL;
            if (!NumericMatrix::is_na(val)) {
                lref(i, j) = solve_lref(val, gref(i, j), x(i, j),
                    albin(i, j));
            }
            else {
                lref(i, j) = NA_REAL;
            }
        }
    }
    return lref;
}
// [[Rcpp::export]]
IntegerMatrix getsoiltypecpp(NumericMatrix bulkden, NumericMatrix clay,
    NumericMatrix sand, NumericMatrix silt)
{
    double bdensd = 0.08835813;
    double claysd = 0.1378817;
    double sandsd = 0.237984;
    double siltsd = 0.192472;
    NumericVector bulkv = { 1.597779, 1.587082, 1.578984, 1.513506,
        1.358636, 1.617506, 1.529643, 1.472509, 1.642237, 1.670682,
        1.604273 };
    NumericVector clayv = { 0, 0.075, 0.15, 0.25, 0.175, 0.275,
        0.325, 0.275, 0.3, 0.35, 0.5 };
    NumericVector sandv = { 0.5, 0.8, 0.7, 0.5, 0.35, 0.5, 0.75, 0.3, 
        0.55, 0.1, 0.1 };
    NumericVector siltv = { 0, 0.15, 0.2, 0.3, 0.65, 0.2, 0.3, 0.5, 
        0.1, 0.1, 0.1 };
    // Loop through and calculate most likely soil type
     // Get dimensions
    int rows = bulkden.nrow();
    int cols = bulkden.ncol();
    IntegerMatrix soiltype(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = bulkden(i, j);
            if (NumericMatrix::is_na(clay(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(sand(i, j))) val = NA_REAL;
            if (NumericMatrix::is_na(silt(i, j))) val = NA_REAL;
            if (!NumericMatrix::is_na(val)) {
                double difs = 9999.99;
                for (int k = 0; k < 11; ++k) {
                    double difn = abs(val - bulkv[k]) / bdensd +
                        abs(clay(i, j) - clayv[k]) / claysd +
                        abs(sand(i, j) - sandv[k]) / sandsd +
                        abs(silt(i, j) - siltv[k]) / sandsd;
                    if (difn < difs) soiltype(i, j) = k + 1;
                    difs = difn;

                } // end k
            } // end NA check
            else {
                soiltype(i, j) = NA_REAL;
            } // end NA check
        } // end j
    } // end i
    return soiltype;
}
// [[Rcpp::export]]
NumericMatrix calcclumpcppone(NumericMatrix leafd, NumericMatrix hgt, 
    NumericMatrix lai, double n)
{
    int rows = leafd.nrow();
    int cols = leafd.ncol();
    // leaf size x 1
    NumericMatrix clump(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = leafd(i, j);
            if (!NumericMatrix::is_na(val)) {
                double lyrs = hgt(i, j) / (leafd(i, j) * n);
                double lvs = lai(i, j) / (leafd(i, j) * n);
                double pr = lvs / lyrs;
                if (pr > 1.0) pr = 1.0;
                if (pr < 0.0) pr = 0.0;
                clump(i, j) = pow((1 - pr), lvs);
                //clump(i, j) = pr;
            }
            else {
                clump(i, j) = NA_REAL;
            }
        }
    }
    return clump;
}
// [[Rcpp::export]]
NumericMatrix calcclumpcpp(NumericMatrix leafd, NumericMatrix hgt,
    NumericMatrix lai, double mx = 0.48)
{
    NumericMatrix cl1 = calcclumpcppone(leafd, hgt, lai, 1.0);
    NumericMatrix cl2 = calcclumpcppone(leafd, hgt, lai, 3.0);
    NumericMatrix cl3 = calcclumpcppone(leafd, hgt, lai, 10.0);
    int rows = cl1.nrow();
    int cols = cl2.ncol();
    NumericMatrix clump(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = cl1(i, j);
            if (!NumericMatrix::is_na(val)) {
                clump(i, j) = 0.6 * cl1(i, j) + 0.3 * cl2(i, j)
                    + 0.1 * cl3(i, j);
                if (clump(i, j) > mx) clump(i, j) = mx;
            }
            else {
                clump(i, j) = NA_REAL;
            }
        }
    }
    return clump;
}