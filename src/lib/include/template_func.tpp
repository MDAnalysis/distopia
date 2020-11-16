#ifndef XDIST_DISTOPIA_TEMPLATE_H
#define XDIST_DISTOPIA_TEMPLATE_H

template <typename T>
T SingleTemplatePairwiseDistance(const T *coords1, const T *coords2,
                                     const T *box) {
  T dx = 0.0;

  for (unsigned char i = 0; i < 3; ++i) {
    T rij = coords1[i] - coords2[i];
    T adj = round(rij / box[i]);
    rij -= adj * box[i];
    dx += rij * rij;
  }
  // sqrt must be defined for T
  return sqrt(dx);
}

#endif //XDIST_DISTOPIA_TEMPLATE_H