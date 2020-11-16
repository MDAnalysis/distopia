#ifndef XDIST_DISTOPIA_TEMPLATE_H
#define XDIST_DISTOPIA_TEMPLATE_H


template <typename Numeric>
// can have more than one type if required
// the word numeric is just descriptive and is often denoted as T
Numeric SingleTemplatePairwiseDistance(const Numeric *coords1, const Numeric *coords2,
                                 const Numeric *box) {
  // check to make sure Numeric is actually numeric otherwise dont compile
  // NOTE this constexpr is checked at COMPILE time.
  static_assert(std::is_arithmetic<Numeric>::value, "Requires numeric type");
  // okay now we are good.
  Numeric dx = 0.0;
  for (unsigned char i = 0; i < 3; ++i) {
    Numeric rij = coords1[i] - coords2[i];
    Numeric adj = round(rij / box[i]);
    rij -= adj * box[i];
    dx += rij * rij;
  }
  // sqrt must be defined for T or linker will have a fit
  return sqrt(dx);
}

#endif // XDIST_DISTOPIA_TEMPLATE_H