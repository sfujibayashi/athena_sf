class DynamicMetricDriver {
 public:
  virtual ~DynamicMetricDriver() = default;

  virtual void Update() = 0;
  virtual void UpdateBlackHoleMass(int stage) = 0;

  virtual Real GetBlackHoleMassAccretionRate() const = 0;
  
};
