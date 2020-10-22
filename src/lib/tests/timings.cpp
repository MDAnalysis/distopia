#include <iostream>
#include <math.h>
#include <chrono>

#include "distopia.h"  // a fancy approach
#include "vanilla.h"  // a naive approach
#include "calc_distances.h"  // from mdanalysis
#include "distancekernels.h"  // from mdtraj
#include "anglekernels.h"  // from mdtraj

bool loadHeader(FILE* fp, int* Ncoords, float* box) {
  // header format:
  // natoms
  // boxx, boxy, boxz
  char tmp[1024];

  if (!fgets(tmp, 1024, fp))
    abort();

  *Ncoords = strtol(tmp, NULL, 10);

  fgets(tmp, 1024, fp);
  char* next = tmp;
  for (unsigned i=0; i<3; ++i)
    *box++ = strtof(next, &next);

  return true;
}

bool loadCoords(FILE* fp, int Ncoords, float* coords) {
  char tmp[4096];

  for (unsigned int i=0; i<Ncoords; ++i) {
    fgets(tmp, 4096, fp);
    char* next = tmp;
    for (unsigned char j=0; j<3; ++j)
      *coords++ = strtof(next, &next);
  }

  return true;
}

#define TOL 0.0005 // we should pay attention to this tol
static bool verify(const float* ref, const float* other, unsigned int Ncoords) {
  for (unsigned int i=0; i<Ncoords; ++i)
    if (fabs(ref[i] - other[i]) > TOL) {
      printf("wrong at pos %d\n", i);
      return false;
    }
  return true;
}

int main(int argc, char* argv[]) {
  // usage: file.in
  if (argc <=1) {
    printf("Too few arguments, please supply a coordinate file as a command line argument.\n");
    return(0);
  }
  char* fname = argv[1];

  float box[3];
  float *coords, *coords1, *coords2, *coords3, *results;
  int Ncoords=0;

  std::cout << "\nBEGIN TIMINGS\n";

  FILE* fp = fopen(fname, "r");
  if (!fp)
      return 1;
   if(!loadHeader(fp, &Ncoords, box))
       return 2;

  coords = (float*) malloc(Ncoords * 3 * sizeof(float));
  results = (float*) malloc(Ncoords * sizeof(float) /2);

  loadCoords(fp, Ncoords, coords);

  std::cout << "Read " << Ncoords << " coordinates \n";


  // DISTANCES
  // split coordinates in half
  if (Ncoords % 2 != 0) {
    std::cout << "Ncoords "<< Ncoords << "are not able to be split into 2 \n";
    return 1;
  }

  std::cout << "\nDISTANCES\n";

  coords1 = coords;
  coords2 = coords + (3*Ncoords/2);
  int Nresults = Ncoords/2;

  std::chrono::steady_clock::time_point t1, t2;
  std::chrono::duration<double> dt;

  t1 = std::chrono::steady_clock::now();

  VanillaCalcBonds(coords1, coords2, box, Nresults, results);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);
  std::cout << "Regular calc_bonds:    " << dt.count() << "\n";
  std::cout << "per result regular:    " << dt.count()/Nresults << "\n";


  float* ref_results = (float*) malloc(sizeof(float) * Nresults);
  memcpy(ref_results, results, sizeof(float) * Nresults);

  t1 = std::chrono::steady_clock::now();

  CalcBondsOrtho(coords1, coords2, box, Nresults, results);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);

  std::cout << "XMM calc_bonds:        " << dt.count() << "\n";
  std::cout << "per result XMM:        " << dt.count()/Nresults << "\n";


  if (!verify(ref_results, results, Nresults))
    std::cout << "XMM result wrong!\n";
  else
    std::cout << "XMM Results verified\n";

  t1 = std::chrono::steady_clock::now();

  _calc_bond_distance_ortho((coordinate*)coords1,
                            (coordinate*)coords2, Nresults, box, results);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);

  std::cout << "MDA calc_bonds:        " << dt.count() << "\n";
  std::cout << "per result MDA:        " << dt.count()/Nresults << "\n";


  if (!verify(ref_results, results, Nresults))
    std::cout << "MDA result wrong!\n";
  else
    std::cout << "MDA Results verified\n";

  t1 = std::chrono::steady_clock::now();

  dist_mic(coords1,
           coords2, box, results, Nresults);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);

  std::cout << "MDtraj calc_bonds:     " << dt.count() << "\n";
  std::cout << "per result MDtraj:     " << dt.count()/Nresults << "\n";


  if (!verify(ref_results, results, Nresults))
    std::cout << "MDtraj result wrong!\n";
  else
    std::cout << "MDtraj Results verified\n";

  // ANGLES
  // split coordinates in three

  std::cout << "\nANGLES\n";

  if (Ncoords % 3 != 0) {
    std::cout << "Ncoords "<< Ncoords << "are not able to be split into 3 \n";
    return 1;
  }

  coords1 = coords;
  coords2 = coords + (3*Ncoords/3);
  coords3 = coords + (6*Ncoords/3);
  Nresults = Ncoords/3;
  
  // these are not strictly nessecary (should we keep an explicit handle)
  results = (float*) realloc(results, Nresults * sizeof(float));
  ref_results = (float*) realloc(ref_results, Nresults * sizeof(float));


  t1 = std::chrono::steady_clock::now();
  // seems sensitive to roundoff
  VanillaCalcAngles(coords1, coords2, coords3, box, Nresults, results);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);
  std::cout << "Regular calc_angles:    " << dt.count() << "\n";
  std::cout << "per result calc_angles: " << dt.count()/Nresults << "\n";
  memcpy(ref_results, results, sizeof(float) * Nresults);

  t1 = std::chrono::steady_clock::now();

  _calc_angle_ortho((coordinate*)coords1,
                            (coordinate*)coords2, (coordinate*)coords3, Nresults, box, results);
  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);
  std::cout << "MDA calc_angles:        " << dt.count() << "\n";
  std::cout << "per result MDA:         " << dt.count()/Nresults << "\n";
  
  if (!verify(ref_results, results, Nresults)) {
    std::cout << "MDA result wrong!\n";
  }
  else {
    std::cout << "MDA Results verified\n";
  }

  t1 = std::chrono::steady_clock::now();

  angle_mic(coords, triplets, box, results, 1, Nresults);

  t2 = std::chrono::steady_clock::now();

  dt = (t2 - t1);

  std::cout << "MDtraj calc_bonds:     " << dt.count() << "\n";
  std::cout << "per result MDtraj:     " << dt.count()/Nresults << "\n";
  

  return 0;
}