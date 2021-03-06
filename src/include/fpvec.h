#ifndef FPVEC_H
#define FPVEC_H
#include <vector>

using std::vector;
using vpvec = vector<uint32_t>;

extern double to_double (const vpvec &);
extern double to_double (const uint32_t *, int);
extern void set (vpvec &, double);

extern QString to_string (vpvec, int prec);
extern vpvec from_string (QString, int sz);

struct invalid_decimal_string
{
};

extern vpvec add (const vpvec &srca, const vpvec &srcb);
extern vpvec sub (const vpvec &srca, const vpvec &srcb);
extern vpvec div1 (const vpvec &src, int32_t d);
extern vpvec mul1 (const vpvec &src, uint32_t d);

#endif
