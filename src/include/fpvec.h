#ifndef FPVEC_H
#define FPVEC_H
#include <vector>

using std::vector;
using vpvec = vector<uint32_t>;

extern double to_double (const vpvec &);
extern void set (vpvec &, double);

extern vpvec add (const vpvec &srca, const vpvec &srcb);
extern vpvec sub (const vpvec &srca, const vpvec &srcb);
extern vpvec div1 (const vpvec &src, int32_t d);
extern vpvec mul1 (const vpvec &src, uint32_t d);

#endif
