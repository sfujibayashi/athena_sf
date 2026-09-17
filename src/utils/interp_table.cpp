//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file interp_table.cpp
//! \brief implements functions in class InterpTable2D an intpolated lookup table

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../coordinates/coordinates.hpp" // Coordinates
#include "interp_table.hpp"

//! A contructor that setts the size of the table with number of variables nvar
//! and dimensions nx2 x nx1 (interpolated dimensions)
InterpTable2D::InterpTable2D(const int nvar, const int nx2, const int nx1) {
  SetSize(nvar, nx2, nx1);
}

//! Set size of table
void InterpTable2D::SetSize(const int nvar, const int nx2, const int nx1) {
  nvar_ = nvar; // number of variables/tables
  nx2_ = nx2; // slower indexing dimension
  nx1_ = nx1; // faster indexing dimension
  data.NewAthenaArray(nvar, nx2, nx1);
}

//! Set the corrdinate limits for x1
void InterpTable2D::SetX1lim(Real x1min, Real x1max) {
  x1min_ = x1min;
  x1max_ = x1max;
  x1norm_ = (nx1_ - 1) / (x1max - x1min);
}

//! Set the corrdinate limits for x2
void InterpTable2D::SetX2lim(Real x2min, Real x2max) {
  x2min_ = x2min;
  x2max_ = x2max;
  x2norm_ = (nx2_ - 1) / (x2max - x2min);
}

void InterpTable2D::GetX1lim(Real &x1min, Real &x1max) {
  x1min = x1min_;
  x1max = x1max_;
}


void InterpTable2D::GetX2lim(Real &x2min, Real &x2max) {
  x2min = x2min_;
  x2max = x2max_;
}

void InterpTable2D::GetSize(int &nvar, int &nx2, int &nx1) {
  nvar = nvar_;
  nx2 = nx2_;
  nx1 = nx1_;
}

//! Bilinear interpolation
Real InterpTable2D::interpolate(int var, Real x2, Real x1) {
  Real x, y, xrl, yrl, out;
  x = (x2 - x2min_) * x2norm_;
  y = (x1 - x1min_) * x1norm_;
  int xil = static_cast<int>(x); // lower x index
  int yil = static_cast<int>(y); // lower y index
  int nx = nx2_;
  int ny = nx1_;
  // if off table, do linear extrapolation
  if (xil < 0) { // below xmin
    xil = 0;
  } else if (xil >= nx - 1) { // above xmax
    xil = nx - 2;
  }
  xrl = 1 + xil - x;  // x residual

  if (yil < 0) { // below ymin
    yil = 0;
  } else if (yil >= ny - 1) { // above ymax
    yil = ny - 2;
  }
  yrl = 1 + yil - y;  // y residual

  // Sample from the 4 nearest data points and weight appropriately
  // data is an attribute of the eos class
  out =   xrl  *  yrl  *data(var, xil , yil )
          +   xrl  *(1-yrl)*data(var, xil ,yil+1)
          + (1-xrl)*  yrl  *data(var,xil+1, yil )
          + (1-xrl)*(1-yrl)*data(var,xil+1,yil+1);
  return out;
}

//! A contructor that sets the size of the table with number of variables nvar
//! and dimensions nx3 x nx2 x nx1 (interpolated dimensions)
InterpTable3D::InterpTable3D(const int nvar, const int nx3, const int nx2, const int nx1) {
  SetSize(nvar, nx3, nx2, nx1);
}

//! Set size of table
void InterpTable3D::SetSize(const int nvar, const int nx3, const int nx2, const int nx1) {
  nvar_ = nvar;
  nx3_ = nx3; // the slowest indexing dimension
  nx2_ = nx2;
  nx1_ = nx1; // fastest dimension

  data.NewAthenaArray(nvar, nx3, nx2, nx1);
}

//! Set the corrdinate limits for x1
void InterpTable3D::SetX1lim(Real x1min, Real x1max) {
  x1min_ = x1min;
  x1max_ = x1max;
  x1norm_ = (nx1_ - 1) / (x1max - x1min);
}

//! Set the corrdinate limits for x2
void InterpTable3D::SetX2lim(Real x2min, Real x2max) {
  x2min_ = x2min;
  x2max_ = x2max;
  x2norm_ = (nx2_ - 1) / (x2max - x2min);
}

//! Set the corrdinate limits for x3
void InterpTable3D::SetX3lim(Real x3min, Real x3max) {
  x3min_ = x3min;
  x3max_ = x3max;
  x3norm_ = (nx3_ - 1) / (x3max - x3min);
}

void InterpTable3D::GetX1lim(Real &x1min, Real &x1max) {
  x1min = x1min_;
  x1max = x1max_;
}

void InterpTable3D::GetX2lim(Real &x2min, Real &x2max) {
  x2min = x2min_;
  x2max = x2max_;
}

void InterpTable3D::GetX3lim(Real &x3min, Real &x3max) {
  x3min = x3min_;
  x3max = x3max_;
}

void InterpTable3D::GetSize(int &nvar, int &nx3, int &nx2, int &nx1) {
  nvar = nvar_;
  nx3 = nx3_;
  nx2 = nx2_;
  nx1 = nx1_;
}

//! Tri-linear interpolation
Real InterpTable3D::interpolate(int var, Real x3, Real x2, Real x1) {
  Real x, y, z, xrl, yrl, zrl, out;
  x = (x3 - x3min_) * x3norm_;
  y = (x2 - x2min_) * x2norm_;
  z = (x1 - x1min_) * x1norm_;
  int xil = static_cast<int>(x); // lower x index
  int yil = static_cast<int>(y); // lower y index
  int zil = static_cast<int>(z); // lower z index
  int nx = nx3_;
  int ny = nx2_;
  int nz = nx1_;
  // if off table, do linear extrapolation
  if (xil < 0) { // below xmin
    xil = 0;
  } else if (xil >= nx - 1) { // above xmax
    xil = nx - 2;
  }
  xrl = 1 + xil - x;  // x residual

  if (yil < 0) { // below ymin
    yil = 0;
  } else if (yil >= ny - 1) { // above ymax
    yil = ny - 2;
  }
  yrl = 1 + yil - y;  // y residual

  if (zil < 0) { // below zmin
    zil = 0;
  } else if (zil >= nz - 1) { // above zmax
    zil = nz - 2;
  }
  zrl = 1 + zil - z;  // z residual

  // Sample from the 8 nearest data points and weight appropriately
  // data is an attribute of the eos class
  out =        xrl  *   yrl *   zrl * data(var,xil  ,yil  ,zil  )
          +    xrl  *(1-yrl)*   zrl * data(var,xil  ,yil+1,zil  )
          + (1-xrl) *   yrl *   zrl * data(var,xil+1,yil  ,zil  )
          + (1-xrl) *(1-yrl)*   zrl * data(var,xil+1,yil+1,zil  )
          +    xrl  *   yrl *(1-zrl)* data(var,xil  ,yil  ,zil+1)
          +    xrl  *(1-yrl)*(1-zrl)* data(var,xil  ,yil+1,zil+1)
          + (1-xrl) *   yrl *(1-zrl)* data(var,xil+1,yil  ,zil+1)
          + (1-xrl) *(1-yrl)*(1-zrl)* data(var,xil+1,yil+1,zil+1);

  return out;
}
