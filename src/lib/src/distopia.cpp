//
// Created by Richard Gowers on 8/13/20.
//
#include "distopia.h"

// Orthogonal boxes
#define BOXTYPE 1
#include "distopia.tpl"
#undef BOXTYPE
#undef MIC
#undef SINGLEDISTANCE
#undef SINGLEANGLE

// No periodic boundaries
#define BOXTYPE 2
#include "distopia.tpl"
#undef BOBOXTYPE