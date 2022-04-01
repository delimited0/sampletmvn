#ifndef TUVN
#define TUVN

class Tuvn {
public:
  Tuvn() {}
  
  virtual double rtuvn(const double mean, const double sd, 
                       const double lower, const double upper) = 0;
};

class LG2015Tuvn : public Tuvn {
public:
  LG2015Tuvn() : Tuvn() {}
  
  double rtuvn(const double mean, const double sd, 
               const double lower, const double upper) override;
};

class BE2017Tuvn : public Tuvn {
public:
  BE2017Tuvn() : Tuvn() {}
  
  double rtuvn(const double mean, const double sd, 
               const double lower, const double upper) override;  
};

double rtuvn_lg2015(double mean, double sd, double lower, double upper);
double rtuvn_be2017(double mean, double sd, double lower, double upper);

#endif