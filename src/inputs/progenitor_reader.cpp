//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file progenitor_reader.cpp
//! \brief Reader for standardized 1D progenitor profiles stored in HDF5.

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "../athena.hpp"
#include "../defs.hpp"
#include "progenitor_reader.hpp"

#ifdef HDF5OUTPUT

#include <hdf5.h>

namespace {

//----------------------------------------------------------------------------------------
// Check whether a dataset exists.

bool DatasetExists(hid_t file, const std::string &name) {
  return (H5Lexists(file, name.c_str(), H5P_DEFAULT) > 0);
}

//----------------------------------------------------------------------------------------
// Read a rank-1 HDF5 dataset.
// HDF5 performs conversion to native double; values are then cast to Real.

std::vector<Real> Read1DRealDataset(hid_t file,
                                    const std::string &filename,
                                    const std::string &name,
                                    std::size_t expected_size = 0) {
  hid_t dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Could not open dataset '" << name << "' in file '"
        << filename << "'." << std::endl;
    ATHENA_ERROR(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  if (dataspace < 0) {
    H5Dclose(dataset);

    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Could not obtain dataspace for dataset '" << name
        << "' in file '" << filename << "'." << std::endl;
    ATHENA_ERROR(msg);
  }

  const int rank = H5Sget_simple_extent_ndims(dataspace);
  if (rank != 1) {
    H5Sclose(dataspace);
    H5Dclose(dataset);

    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Dataset '" << name << "' in file '" << filename
        << "' must have rank 1, but rank is " << rank << "."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  hsize_t dims[1];
  H5Sget_simple_extent_dims(dataspace, dims, nullptr);

  const std::size_t n = static_cast<std::size_t>(dims[0]);

  if (expected_size != 0 && n != expected_size) {
    H5Sclose(dataspace);
    H5Dclose(dataset);

    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Dataset '" << name << "' in file '" << filename
        << "' has length " << n
        << ", while expected length is " << expected_size << "."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  std::vector<double> buffer(n);

  const herr_t status =
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
              H5P_DEFAULT, buffer.data());

  H5Sclose(dataspace);
  H5Dclose(dataset);

  if (status < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Failed to read dataset '" << name << "' from file '"
        << filename << "'." << std::endl;
    ATHENA_ERROR(msg);
  }

  std::vector<Real> result(n);
  for (std::size_t i = 0; i < n; ++i) {
    result[i] = static_cast<Real>(buffer[i]);
  }

  return result;
}

//----------------------------------------------------------------------------------------
// Check that all values in an array are finite.

void CheckFinite(const std::vector<Real> &data,
                 const std::string &filename,
                 const std::string &name) {
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (!std::isfinite(data[i])) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
          << "Non-finite value found in dataset '" << name
          << "' in file '" << filename << "' at index " << i << "."
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}

//----------------------------------------------------------------------------------------
// Check radial ordering.

void CheckMonotonic(const std::vector<Real> &data,
                    const std::string &filename,
                    const std::string &name) {
  for (std::size_t i = 1; i < data.size(); ++i) {
    if (!(data[i] > data[i-1])) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
          << "Dataset '" << name << "' in file '" << filename
          << "' must be strictly increasing (center -> surface)." << std::endl
          << "At indices " << i-1 << " and " << i << ": "
          << data[i-1] << ", " << data[i] << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}

//----------------------------------------------------------------------------------------
// Check quantities that must be positive.

void CheckPositive(const std::vector<Real> &data,
                   const std::string &filename,
                   const std::string &name) {
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (!(data[i] > 0.0)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
          << "Dataset '" << name << "' in file '" << filename
          << "' contains a non-positive value at index " << i
          << ": " << data[i] << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn ProgenitorProfile ReadProgenitorProfile(const std::string &filename)
//! \brief Read and validate a standardized progenitor profile.

ProgenitorProfile ReadProgenitorProfile(const std::string &filename) {
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Could not open progenitor file '" << filename << "'."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // Required datasets
  const char *required[] = {
      "/mass_face",
      "/radius_face",
      "/rho",
      "/press",
      "/csound"
  };

  for (const char *name : required) {
    if (!DatasetExists(file, name)) {
      H5Fclose(file);

      std::stringstream msg;
      msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
          << "Required dataset '" << name << "' is missing from file '"
          << filename << "'." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // At least one rotation variable must be available.
  const bool has_jrot = DatasetExists(file, "/jrot");
  const bool has_omega = DatasetExists(file, "/omega");

  if (!has_jrot && !has_omega) {
    H5Fclose(file);

    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "File '" << filename << "' must contain at least one of"
        << " '/jrot' or '/omega'." << std::endl;
    ATHENA_ERROR(msg);
  }

  ProgenitorProfile profile;

  // Read mass first to determine N+1.
  profile.mass_face = Read1DRealDataset(file, filename, "/mass_face");

  // n = ncell
  const std::size_t nface = profile.mass_face.size();

  if (nface < 3) {
    H5Fclose(file);

    std::stringstream msg;
    msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
        << "Progenitor profile in file '" << filename
        << "' must contain at least two points." << std::endl;
    ATHENA_ERROR(msg);
  }

  const std::size_t ncell = nface - 1;

  profile.nface = static_cast<int>(nface);
  profile.ncell = static_cast<int>(ncell);
  
  // Allocate derived cell-center coordinates
  profile.mass.resize(ncell);
  profile.radius.resize(ncell);

  profile.radius_face = Read1DRealDataset(file, filename, "/radius_face", nface);
  profile.rho    = Read1DRealDataset(file, filename, "/rho", ncell);
  profile.press  = Read1DRealDataset(file, filename, "/press", ncell);
  profile.csound = Read1DRealDataset(file, filename, "/csound", ncell);

  // Rotation
  profile.has_jrot = has_jrot;
  profile.has_omega = has_omega;

  if (has_jrot) {
    profile.jrot = Read1DRealDataset(file, filename, "/jrot", ncell);
  }

  if (has_omega) {
    profile.omega = Read1DRealDataset(file, filename, "/omega", ncell);
  }

  // Optional thermodynamic/composition fields
  profile.has_ye = DatasetExists(file, "/ye");
  if (profile.has_ye) {
    profile.ye = Read1DRealDataset(file, filename, "/ye", ncell);
  }

  profile.has_temp = DatasetExists(file, "/temp");
  if (profile.has_temp) {
    profile.temp = Read1DRealDataset(file, filename, "/temp", ncell);
  }

  H5Fclose(file);

  // Basic validation
  CheckFinite(profile.mass_face, filename, "/mass_face");
  CheckFinite(profile.radius_face, filename, "/radius_face");
  CheckFinite(profile.rho, filename, "/rho");
  CheckFinite(profile.press, filename, "/press");
  CheckFinite(profile.csound, filename, "/csound");

  if (profile.has_jrot) {
    CheckFinite(profile.jrot, filename, "/jrot");
  }
  if (profile.has_omega) {
    CheckFinite(profile.omega, filename, "/omega");
  }
  if (profile.has_ye) {
    CheckFinite(profile.ye, filename, "/ye");
  }
  if (profile.has_temp) {
    CheckFinite(profile.temp, filename, "/temp");
  }

  // Standard format is center -> surface.
  CheckMonotonic(profile.mass_face, filename, "/mass_face");
  CheckMonotonic(profile.radius_face, filename, "/radius_face");

  CheckPositive(profile.rho, filename, "/rho");
  CheckPositive(profile.press, filename, "/press");
  CheckPositive(profile.csound, filename, "/csound");

  for(int i=0; i<ncell; i++){
    profile.mass[i] = 0.5*(profile.mass_face[i]+profile.mass_face[i+1]);
    profile.radius[i] = std::pow(0.5*(std::pow(profile.radius_face[i],3)+std::pow(profile.radius_face[i+1],3)),1.0/3.0);
  }

  return profile;
}

#else  // HDF5OUTPUT

ProgenitorProfile ReadProgenitorProfile(const std::string &filename) {
  std::stringstream msg;
  msg << "### FATAL ERROR in ReadProgenitorProfile" << std::endl
      << "Cannot read progenitor file '" << filename << "' because Athena++ "
      << "was compiled without HDF5 support." << std::endl
      << "Reconfigure Athena++ with -hdf5." << std::endl;
  ATHENA_ERROR(msg);

  return ProgenitorProfile();
}

#endif  // HDF5OUTPUT
