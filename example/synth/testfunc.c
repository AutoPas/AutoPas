/*
 * Copyright 2003-2016 Jeffrey K. Hollingsworth
 *
 * This file is part of Active Harmony.
 *
 * Active Harmony is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Active Harmony is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
 */
#define _XOPEN_SOURCE 500 // Needed for M_PI and M_E

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "testfunc.h"

benchfunc_t f_dejong;
benchfunc_t f_axispar;
benchfunc_t f_axisrot;
benchfunc_t f_rosenbrock;
benchfunc_t f_ackley;
benchfunc_t f_michalewicz;
benchfunc_t f_cos;
benchfunc_t f_flat;
benchfunc_t f_sum;
benchfunc_t f_oka1_1;
benchfunc_t f_oka1_2;
benchfunc_t f_oka2_1;
benchfunc_t f_oka2_2;
benchfunc_t f_vlmop2_1;
benchfunc_t f_vlmop2_2;
benchfunc_t f_vlmop3_1;
benchfunc_t f_vlmop3_2;
benchfunc_t f_vlmop3_3;
benchfunc_t f_kno1_1;
benchfunc_t f_kno1_2;
benchfunc_t f_dtlz1a_1;
benchfunc_t f_dtlz1a_2;
benchfunc_t f_dtlz2a_1;
benchfunc_t f_dtlz2a_2;
benchfunc_t f_dtlz2a_3;
benchfunc_t f_dtlz4a_1;
benchfunc_t f_dtlz4a_2;
benchfunc_t f_dtlz4a_3;
benchfunc_t f_dtlz7a_1;
benchfunc_t f_dtlz7a_2;
benchfunc_t f_dtlz7a_3;

finfo_t flist[] = {
    {"dejong", "De Jong’s first function",
     0, FTYPE_REAL, -64.0, 64.0, 1, 0.0, f_dejong,
     "    De Jong's first function is continuous, convex, and unimodal.\n"},

    {"axispar", "Axis parallel hyper-ellipsoid function",
     0, FTYPE_REAL, -64.0, 64.0, 1, 0.0, f_axispar,
     "    The axis parallel hyper-ellipsoid is similar the first De Jong\n"
     "    function.  It is also known as the weighted sphere model.  This\n"
     "    function is continuous, convex and unimodal.\n"},

    {"axisrot", "Rotated hyper-ellipsoid function",
     0, FTYPE_REAL, -64.0, 64.0, 1, 0.0, f_axisrot,
     "    The rotated hyper-ellipsoid function.\n"},

    {"rosenbrock", "Rosenbrock's Valley",
     0, FTYPE_REAL, -2, 2, 1, 0.0, f_rosenbrock,
     "    De Jong's second function, Rosenbrock's valley, is a classic\n"
     "    optimization problem, also known as Rosenbrock's banana function.\n"
     "    The global optimum is inside a long, narrow, parabolic shaped\n"
     "    flat valley.  To find the valley is trivial, however convergence\n"
     "    to the global optimum is difficult.\n"},

    {"ackley", "Ackley's Function",
     0, FTYPE_REAL, -32.0, 32.0, 1, 0.0, f_ackley,
     "    Ackley's Function is a continuous, multimodal function obtained\n"
     "    by modulating an exponential function with a cosine wave of\n"
     "    moderate amplitude.  Its topology is characterized by an almost\n"
     "    flat (due to the dominating exponential) outer region and a\n"
     "    central hole or peak where the modulations by the cosine wave\n"
     "    become more and more influential.\n"},

    {"michalewicz", "Michalewicz's Function",
     0, FTYPE_REAL, 0, 4.0, 0, -0.0, f_michalewicz,
     "    Michalewicz's Function is a multimodal function with N! local\n"
     "    optima.  It takes one parameter that controls the slope of its\n"
     "    valleys and edges.  As the parameter grows larger, the function\n"
     "    becomes more difficult and effectively becomes a needle in the\n"
     "    haystack problem.  10 is used by default.\n"},

    {"cos", "Cosine Sum",
     0, FTYPE_REAL, -10.0, 10.0, 0, -0.0, f_cos,
     "    Returns the running sum of the cosine of each parameter.\n"},

    {"flat", "Flat Function",
     0, FTYPE_REAL, -10.0, 10.0, 0, -0.0, f_flat,
     "    Always returns 1.0 (or the value of the first option).\n"},

    {"sum", "Sum",
     0, FTYPE_REAL, -10.0, 10.0, 0, -0.0, f_sum,
     "    Returns the simple sum of all parameters.\n"},

    {"oka1_1", "oka1 Objective #1",
     2, FTYPE_REAL, 0.027415472, 6.310535189, 0, -0.0, f_oka1_1,
     "    oka1 from J.Knowles.\n"},

    {"oka1_2", "oka1 Objective #2",
     2, FTYPE_REAL, -0.028709416, 5.999937366, 0, -0.0, f_oka1_2,
     "    oka1 from J.Knowles.\n"},

    {"oka2_1", "oka2 Objective #1",
     3, FTYPE_REAL, -5.0, 5.0, 0, -0.0, f_oka2_1,
     "    oka2 from J.Knowles.\n"},

    {"oka2_2", "oka2 Objective #2",
     3, FTYPE_REAL, -5.0, 5.0, 0, -0.0, f_oka2_2,
     "    oka2 from J.Knowles.\n"},

    {"vlmop2_1", "vlmop2 Objective #1",
     2, FTYPE_REAL, -2.0, 2.0, 0, -0.0, f_vlmop2_1,
     "    vlmop2 from J.Knowles.\n"},

    {"vlmop2_2", "vlmop2 Objective #2",
     2, FTYPE_REAL, -2.0, 2.0, 0, -0.0, f_vlmop2_2,
     "    vlmop2 from J.Knowles.\n"},

    {"vlmop3_1", "vlmop3 Objective #1",
     2, FTYPE_REAL, -3.0, 3.0, 0, -0.0, f_vlmop3_1,
     "    vlmop2 from J.Knowles.\n"},

    {"vlmop3_2", "vlmop3 Objective #2",
     2, FTYPE_REAL, -3.0, 3.0, 0, -0.0, f_vlmop3_2,
     "    vlmop2 from J.Knowles.\n"},

    {"vlmop3_3", "vlmop3 Objective #3",
     2, FTYPE_REAL, -3.0, 3.0, 0, -0.0, f_vlmop3_3,
     "    vlmop2 from J.Knowles.\n"},

    {"kno1_1", "kno1 Objective #1",
     2, FTYPE_REAL, 0.0, 3.0, 0, -0.0, f_kno1_1,
     "    kno1 from J.Knowles.\n"},

    {"kno1_2", "kno1 Objective #2",
     2, FTYPE_REAL, 0.0, 3.0, 0, -0.0, f_kno1_2,
     "    kno1 from J.Knowles.\n"},

    {"dtlz1a_1", "dtlz1a Objective #1",
     6, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz1a_1,
     "    dtlz1a from J.Knowles.  Uses 3 decision variables.\n"},

    {"dtlz1a_2", "dtlz1a Objective #2",
     6, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz1a_2,
     "    dtlz1a from J.Knowles.  Uses 3 decision variables.\n"},

    {"dtlz2a_1", "dtlz2a Objective #1",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz2a_1,
     "    dtlz2a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz2a_2", "dtlz2a Objective #2",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz2a_2,
     "    dtlz2a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz2a_3", "dtlz2a Objective #3",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz2a_3,
     "    dtlz2a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz4a_1", "dtlz4a Objective #1",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz4a_1,
     "    dtlz4a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz4a_2", "dtlz4a Objective #2",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz4a_2,
     "    dtlz4a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz4a_3", "dtlz4a Objective #3",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz4a_3,
     "    dtlz4a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz7a_1", "dtlz7a Objective #1",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz7a_1,
     "    dtlz7a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz7a_2", "dtlz7a Objective #2",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz7a_2,
     "    dtlz7a from J.Knowles.  Uses 8 decision variables.\n"},

    {"dtlz7a_3", "dtlz7a Objective #3",
     8, FTYPE_REAL, 0.0, 1.0, 0, -0.0, f_dtlz7a_3,
     "    dtlz7a from J.Knowles.  Uses 8 decision variables.\n"},

    {NULL, NULL, 0, FTYPE_REAL, 0.0, 0.0, 0, -0.0, NULL, NULL}
};

void flist_print(FILE* fd, int verbose)
{
    finfo_t* ptr;
    int len;

    len = 0;
    if (!verbose) {
        for (ptr = flist; ptr->name; ++ptr) {
            if (len < strlen(ptr->name))
                len = strlen(ptr->name);
        }
    }

    for (ptr = flist; ptr->name; ++ptr) {
        fprintf(fd, "%-*s - %s - [%.3lf, %.3lf]",
                len, ptr->name, ptr->title,
                ptr->b_min, ptr->b_max);
        if (ptr->n_max)
            fprintf(fd, " (N<=%d)", ptr->n_max);
        fprintf(fd, "\n");

        if (verbose)
            fprintf(fd, "%s\n", ptr->description);
    }
}

finfo_t* flist_find(const char* name)
{
    finfo_t* ptr;

    for (ptr = flist; ptr->name; ++ptr) {
        if (strcmp(ptr->name, name) == 0)
            break;
    }
    if (ptr->name)
        return ptr;
    return NULL;
}

double f_dejong(int n, double x[], double option[])
{
    int i;
    double d;

    d = 0.0;
    for (i = 0; i < n; ++i)
        d += x[i] * x[i];

    return d;
}

double f_axispar(int n, double x[], double option[])
{
    int i;
    double d;

    d = 0.0;
    for (i = 0; i < n; ++i)
        d += i * (x[i] * x[i]);

    return d;
}

double f_axisrot(int n, double x[], double option[])
{
    int i, j;
    double d;

    d = 0.0;
    for (i = 0; i < n; ++i)
        for (j = 0; j < i; ++j)
            d += x[j] * x[j];

    return d;
}

double f_rosenbrock(int n, double x[], double option[])
{
    int i;
    double d, d1, d2;

    if (n == 2) {
        d1 = x[1] - (x[0] * x[0]);
        d2 = 1.0 - x[0];
        d  = (100.0 * d1 * d1) + (d2 * d2);
    }
    else {
        d = 0.0;
        for (i = 1; i < n; ++i) {
            d1 = x[i] - (x[i-1] * x[i-1]);
            d2 = 1.0 - x[i-1];
            d += (100.0 * d1 * d1) + (d2 * d2);
        }
    }
    return d;
}

double f_ackley(int n, double x[], double option[])
{
    int i;
    double a, b, c;
    double d1, d2;

    if (option) {
        a = option[0];
        b = option[1];
        c = option[2];
    }
    else {
        a = 20.0;
        b = 0.2;
        c = 2.0 * M_PI;
    }

    d1 = d2 = 0.0;
    for (i = 0; i < n; ++i) {
        d1 += x[i] * x[i];
        d2 += cos(c * x[i]);
    }

    return -a * exp(-b * sqrt(d1 / n)) - exp(d2 / n) + a + M_E;
}

double f_michalewicz(int n, double x[], double option[])
{
    int i;
    double m = 10.0;
    double sum = 0.0;

    if (option)
        m = option[0];

    for (i = 0; i < n; ++i) {
        double d1 = (i+1) * x[i] * x[i];
        double d2 = sin(d1 / M_PI);

        sum += sin(x[i]) * pow(d2, 2 * m);
    }
    return -sum;
}

double f_cos(int n, double x[], double option[])
{
    int i;
    double sum = 0.0;

    for (i = 0; i < n; ++i)
        sum += cos(x[i]);

    if (option)
        sum *= option[0];

    return sum;
}

double f_flat(int n, double x[], double option[])
{
    if (option)
        return option[0];
    return 1.0;
}

double f_sum(int n, double x[], double option[])
{
    int i;
    double sum = 0.0;

    for (i = 0; i < n; ++i)
        sum += x[i];

    if (option)
        sum *= option[0];

    return sum;
}

double f_oka1_1(int n, double x[], double option[])
{
    return cos(M_PI/12.0)*x[0] - sin(M_PI/12.0)*x[1];
}

double f_oka1_2(int n, double x[], double option[])
{
    double x1p = cos(M_PI/12.0)*x[0] - sin(M_PI/12.0)*x[1];
    double x2p = sin(M_PI/12.0)*x[0] + cos(M_PI/12.0)*x[1];

    return (sqrt(2*M_PI) - sqrt(fabs(x1p)) +
            2 * pow(fabs(x2p - 3 * cos(x1p) - 3), 1.0/3.0));
}

double f_oka2_1(int n, double x[], double option[])
{
    // There are 2 objectives, 3 decision (x) variables.
    // The x variables are in the ranges:
    //  x[0] in [-PI, PI]
    //  x[1] in [-5, 5]
    //  x[2] in [-5, 5]
    //
    return x[0];
}

double f_oka2_2(int n, double x[], double option[])
{
    // There are 2 objectives, 3 decision (x) variables.
    // The x variables are in the ranges:
    //  x[0] in [-PI, PI]
    //  x[1] in [-5, 5]
    //  x[2] in [-5, 5]
    //
    return (1.0 - (1.0 / (4.0 * M_PI * M_PI)) * pow(x[0] + M_PI, 2) +
            pow(fabs(x[1] - 5.0 * cos(x[0])), 1.0/3.0) +
            pow(fabs(x[2] - 5.0 * sin(x[0])), 1.0/3.0));
}

double f_vlmop2_1(int n, double x[], double option[])
{
    // x variables must be in the range [-2,2].
    int i;
    double sum1 = 0.0;

    for (i = 0; i < 2; ++i)
        sum1 += pow(x[i] - (1.0/sqrt(2.0)), 2);

    return 1 - exp(-sum1);
}

double f_vlmop2_2(int n, double x[], double option[])
{
    // x variables must be in the range [-2,2].
    int i;
    double sum2 = 0.0;

    for (i = 0; i < 2; ++i)
        sum2 += pow(x[i] + (1.0/sqrt(2.0)), 2.0);

    return 1 - exp(-sum2);
}

double f_vlmop3_1(int n, double x[], double option[])
{
    return 0.5*(x[0]*x[0]+x[1]*x[1]) + sin(x[0]*x[0]+x[1]*x[1]);
}

double f_vlmop3_2(int n, double x[], double option[])
{
    return (pow(3*x[0]-2*x[1]+4.0, 2.0) / 8.0 +
            pow(x[0]-x[1]+1, 2.0)/27.0 + 15.0);
}

double f_vlmop3_3(int n, double x[], double option[])
{
    return (1.0 / (x[0]*x[0]+x[1]*x[1]+1.0) -
            1.1 * exp(-(x[0]*x[0]) - (x[1]*x[1])));
}

double f_kno1_1(int n, double x[], double option[])
{
    // There are 2 objectives, 2 decision (x) variables.
    // The x variables are in the range [0,3]
    //
    double f;
    double g;
    double c;

    c = x[0]+x[1];
    f = 20 - (11.0 + 3.0*sin((5.0*c)*(0.5*c)) +
              3.0*sin(4.0*c) + 5.0*sin(2.0*c+2.0));
    g = (M_PI/2.0)*(x[0]-x[1]+3.0)/6.0;

    return 20-(f*cos(g));
}

double f_kno1_2(int n, double x[], double option[])
{
    // There are 2 objectives, 2 decision (x) variables.
    // The x variables are in the range [0,3]
    //
    double f;
    double g;
    double c;

    c = x[0]+x[1];
    f = 20-( 11+3*sin((5*c)*(0.5*c)) + 3*sin(4*c) + 5 *sin(2*c+2));
    g = (M_PI/2.0)*(x[0]-x[1]+3.0)/6.0;

    return 20-(f*sin(g));
}

double f_dtlz1a_1(int n, double x[], double option[])
{
    // There are 2 objectives, 6 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i;
    double g = 0.0;
    assert(n == 6);

    for (i = 1; i < n; ++i) {
        // Note this is 20*PI in Deb's DTLZ1.
        g += (x[i] - 0.5) * (x[i] - 0.5) - cos(2 * M_PI * (x[i] - 0.5));
    }
    g += n-1;
    g *= 100;

    return 0.5 * (    x[0]) * (1 + g);
}

double f_dtlz1a_2(int n, double x[], double option[])
{
    // There are 2 objectives, 6 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i;
    double g = 0.0;
    assert(n == 6);

    for (i = 1; i < n; ++i) {
        // Note this is 20*PI in Deb's DTLZ1.
        g += (x[i] - 0.5) * (x[i] - 0.5) - cos(2 * M_PI * (x[i] - 0.5));
    }
    g += n-1;
    g *= 100;

    return 0.5 * (1 - x[0]) * (1 + g);
}

double f_dtlz2a_1(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i;
    double alph = 1.0;
    double g = 0.0;
    assert(n == 8);

    for(i = 2; i < n; ++i)
        g += (x[i] - 0.5) * (x[i] - 0.5);

    return ((1 + g) *
            cos(pow(x[0], alph) * M_PI/2) *
            cos(pow(x[1], alph) * M_PI/2));
}

double f_dtlz2a_2(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i;
    double alph = 1.0;
    double g = 0.0;
    assert(n == 8);

    if (option)
        alph = option[0];

    for (i = 2; i < n; ++i)
        g += (x[i] - 0.5) * (x[i] - 0.5);

    return ((1 + g) *
            cos(pow(x[0], alph) * M_PI/2) *
            sin(pow(x[1], alph) * M_PI/2));
}

double f_dtlz2a_3(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i;
    double alph = 1.0;
    double g = 0.0;
    assert(n == 8);

    if (option)
        alph = option[0];

    for (i = 2; i < n; ++i)
        g += (x[i] - 0.5) * (x[i] - 0.5);

    return (1 + g) * sin(pow(x[0], alph) * M_PI/2);
}

double f_dtlz4a_1(int n, double x[], double option[])
{
    int i;
    double g = 0.0;
    for(i=2;i<8;i++)
        g+=(x[i]-0.5)*(x[i]-0.5);

    return (1 + g)*cos(pow(x[0],100.0)*M_PI/2)*cos(pow(x[1],100.0)*M_PI/2);
}

double f_dtlz4a_2(int n, double x[], double option[])
{
    int i;
    double g = 0.0;
    for(i=2;i<8;i++)
        g+=(x[i]-0.5)*(x[i]-0.5);

    return (1 + g)*cos(pow(x[0],100.0)*M_PI/2)*sin(pow(x[1],100.0)*M_PI/2);
}

double f_dtlz4a_3(int n, double x[], double option[])
{
    int i;
    double g = 0.0;
    for(i=2;i<8;i++)
        g+=(x[i]-0.5)*(x[i]-0.5);

    return (1 + g)*sin(pow(x[0],100.0)*M_PI/2);
}

double f_dtlz7a_1(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    assert(n == 8);
    return x[0];
}

double f_dtlz7a_2(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    assert(n == 8);
    return x[1];
}

double f_dtlz7a_3(int n, double x[], double option[])
{
    // There are 3 objectives, 8 decision (x) variables.
    // The x variables are in the range [0,1].
    //
    int i, nobjs = 3;
    double g, h, sum;
    assert(n == 8);

    g = 0.0;
    for (i = 2; i < n; ++i)
        g += x[i];

    g *= 9.0 / (n - nobjs + 1);
    g += 1.0;

    sum = 0.0;
    for (i = 0; i < nobjs - 1; ++i)
        sum += ( x[i]/(1.0+g) * (1.0+sin(3*M_PI*x[i])) );
    h = nobjs - sum;

    return (1 + g) * h;
}
