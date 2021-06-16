
inline double square(double val) {
  return val * val;
}

inline double distanceSquared(double xDiff, double yDiff, double zDiff) {
  return square(xDiff) + square(yDiff) + square(zDiff);
}

inline double calculateLJFactor(double dr2, uniform double sigmasquare, uniform double epsilon24) {
  double invdr2 = 1. / dr2;
  double lj2 = sigmasquare * invdr2;
  double lj6 = lj2 * lj2 * lj2;
  double lj12 = lj6 * lj6;
  double lj12m6 = lj12 - lj6;
  double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
  return fac;
}

export void SoAFunctorPairISPC(uniform int64 numParticles1, uniform int64 numParticles2,
                               uniform int64* uniform ownedStatePtr1, uniform int64* uniform ownedStatePtr2,
                               uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                               uniform double* uniform x1ptr, uniform double* uniform x2ptr,
                               uniform double* uniform y1ptr, uniform double* uniform y2ptr,
                               uniform double* uniform z1ptr, uniform double* uniform z2ptr,
                               uniform double* uniform fx1ptr, uniform double* uniform fx2ptr,
                               uniform double* uniform fy1ptr, uniform double* uniform fy2ptr,
                               uniform double* uniform fz1ptr, uniform double* uniform fz2ptr) {
  for(uniform int i = 0; i < numParticles1; ++i) {
    if(ownedStatePtr1[i] == 0) {
      continue;
    }
    double fxacc = 0;
    double fyacc = 0;
    double fzacc = 0;
    foreach (j = 0 ... numParticles2) {
      if (ownedStatePtr2[j] == 0) {
      continue;
      }
      double drx = x1ptr[i] - x2ptr[j];
      double dry = y1ptr[i] - y2ptr[j];
      double drz = z1ptr[i] - z2ptr[j];
      double dr2 = distanceSquared(drx, dry, drz);

      if (dr2 <= cutoffsquare) {
        double fac = calculateLJFactor(dr2, sigmasquare, epsilon24);

        double fx = drx * fac;
        double fy = dry * fac;
        double fz = drz * fac;

        fx2ptr[j] -= fx;
        fy2ptr[j] -= fy;
        fz2ptr[j] -= fz;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
      }
    }
    fx1ptr[i] += reduce_add(fxacc);
    fy1ptr[i] += reduce_add(fyacc);
    fz1ptr[i] += reduce_add(fzacc);
  }
}

export void SoAFunctorPairISPC2(uniform int64 numParticles1, uniform int64 numParticles2,
    uniform int64* uniform ownedStatePtr1, uniform int64* uniform ownedStatePtr2,
    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
    uniform double* uniform x1ptr, uniform double* uniform x2ptr,
uniform double* uniform y1ptr, uniform double* uniform y2ptr,
uniform double* uniform z1ptr, uniform double* uniform z2ptr,
uniform double* uniform fx1ptr, uniform double* uniform fx2ptr,
uniform double* uniform fy1ptr, uniform double* uniform fy2ptr,
uniform double* uniform fz1ptr, uniform double* uniform fz2ptr) {

  for(uniform int i = 0; i < numParticles1; ++i) {
    if(ownedStatePtr1[i] == 0) {
      continue;
    }
    double fxacc = 0;
    double fyacc = 0;
    double fzacc = 0;
    foreach (j = 0 ... numParticles2) {
      if (ownedStatePtr2[j] == 0) {
        continue;
      }
      double drx = x1ptr[i] - x2ptr[j];
      double dry = y1ptr[i] - y2ptr[j];
      double drz = z1ptr[i] - z2ptr[j];
      double dr2 = distanceSquared(drx, dry, drz);

      if (dr2 <= cutoffsquare) {
        double fac = calculateLJFactor(dr2, sigmasquare, epsilon24);

        double fx = drx * fac;
        double fy = dry * fac;
        double fz = drz * fac;

        fx2ptr[j] -= fx;
        fy2ptr[j] -= fy;
        fz2ptr[j] -= fz;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
      }
    }
    fx1ptr[i] += reduce_add(fxacc);
    fy1ptr[i] += reduce_add(fyacc);
    fz1ptr[i] += reduce_add(fzacc);
  }
}

export void SoAFunctorVerletISPC(uniform size_t indexFirst, uniform size_t* uniform neighborList,
                                 uniform size_t numNeighbors,
                                 uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                                 uniform double* uniform x1ptr, uniform double* uniform y1ptr,
                                 uniform double* uniform z1ptr, uniform int64* uniform ownedStates,
                                 uniform double* uniform fx1ptr, uniform double* uniform fy1ptr,
                                 uniform double* uniform fz1ptr) {
  double fxacc = 0;
  double fyacc = 0;
  double fzacc = 0;
  foreach (i = 0 ... numNeighbors) {
    size_t index = neighborList[i];
    if (ownedStates[index] == 0) {
      continue;
    }
    double drx = x1ptr[index] - x1ptr[indexFirst];
    double dry = y1ptr[index] - y1ptr[indexFirst];
    double drz = z1ptr[index] - z1ptr[indexFirst];
    double dr2 = distanceSquared(drx, dry, drz);

    if (dr2 <= cutoffsquare) {
      double fac = calculateLJFactor(dr2, sigmasquare, epsilon24);

      double fx = drx * fac;
      double fy = dry * fac;
      double fz = drz * fac;

//      fx1ptr[index] -= fx;
//      fy1ptr[index] -= fy;
//      fz1ptr[index] -= fz;

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
    }
  }
  fx1ptr[indexFirst] += reduce_add(fxacc);
  fy1ptr[indexFirst] += reduce_add(fyacc);
  fz1ptr[indexFirst] += reduce_add(fzacc);
}