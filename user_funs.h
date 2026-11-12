#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);


//--------------------------------------------------
// Laboratorium 3

// Funkcja testowa z zewnętrzną funkcją kary.
matrix ff3Ta(matrix, matrix = NAN, matrix = NAN);

// Funkcja testowa z wewnętrzną funkcją kary.
matrix ff3Tb(matrix, matrix = NAN, matrix = NAN);

// Funkcja celu dla problemu lotu piłki.
matrix ff3R(matrix, matrix = NAN, matrix = NAN);

// Układ równań różniczkowych opisujący ruch piłki.
matrix df3(double, matrix, matrix = NAN, matrix = NAN);