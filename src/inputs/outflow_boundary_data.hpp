class OutflowBoundaryData {
public:
  OutflowBoundaryData(const std::string &filename);
  ~OutflowBoundaryData();
  
  Real InterpolateRho(Real time, Real theta) const;
  Real InterpolatePress(Real time, Real theta) const;
  // ...
  
private:
  int nt_, nth_;
  
  std::vector<Real> time_;
  std::vector<Real> theta_;
  
  AthenaArray<Real> rho_;
  AthenaArray<Real> press_;
  AthenaArray<Real> vx_;
  AthenaArray<Real> vy_;
  AthenaArray<Real> vz_;
};
