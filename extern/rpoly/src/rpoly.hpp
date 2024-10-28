// rpoly_ak1.cpp - Program for calculating the roots of a polynomial of real coefficients.
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
// 27 May 2014
//
// The sub-routines listed below are translations of the FORTRAN routines included in RPOLY.FOR,
// posted off the NETLIB site as TOMS/493:
//
// http://www.netlib.org/toms/493
//
// TOMS/493 is based on the Jenkins-Traub algorithm.
//
// To distinguish the routines posted below from others, an _ak1 suffix has been appended to them.
//
// Following is a list of the major changes made in the course of translating the TOMS/493 routines
// to the C++ versions posted below:
// 1) All global variables have been eliminated.
// 2) The "FAIL" parameter passed into RPOLY.FOR has been eliminated.
// 3) RPOLY.FOR solves polynomials of degree up to 100, but does not explicitly state this limit.
//     rpoly_ak1 explicitly states this limit; uses the macro name MAXDEGREE to specify this limit;
//     and does a check to ensure that the user input variable Degree is not greater than MAXDEGREE
//     (if it is, an error message is output and rpoly_ak1 terminates). If a user wishes to compute
//     roots of polynomials of degree greater than MAXDEGREE, using a macro name like MAXDEGREE provides
//     the simplest way of offering this capability.
// 4) All "GO TO" statements have been eliminated.
//
// A small main program is included also, to provide an example of how to use rpoly_ak1. In this
// example, data is input from a file to eliminate the need for a user to type data in via
// the console.

#pragma once

#include <cfloat>
#include <cmath>
#include <array>


namespace rpoly
{
    /// \brief Maximum degree of the polynomial
    static constexpr int MaxDegree = 5;

    /// \brief Maximum degree of the polynomial plus one, the number of coefficients in a degree polynomial
    static constexpr int MaxDegreeP1 = MaxDegree + 1;

    /// \brief Calculated the roots of a polynomial of real coefficients
    void rpoly_ak1(const std::array<double, MaxDegreeP1>& op,
                   int& Degree,
                   std::array<double, MaxDegree>& zeror,
                   std::array<double, MaxDegree>& zeroi);
} // namespace rpoly
