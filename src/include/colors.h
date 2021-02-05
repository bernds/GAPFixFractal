#ifndef COLORS_H
#define COLORS_H

#include <QVector>
#include <cstdint>

extern QVector<uint32_t> interpolate_colors (const QVector<uint32_t> &src, int steps,
					     bool narrow_blacks = false, bool narrow_whites = false, int nfactor = 1);

#endif
