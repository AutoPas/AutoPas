
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

inline void SoAFunctorPairISPC(uniform int64 numParticles1, uniform int64 numParticles2,
                               uniform int64* uniform ownedStatePtr1, uniform int64* uniform ownedStatePtr2,
                               uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                               uniform double* uniform x1ptr, uniform double* uniform x2ptr,
                               uniform double* uniform y1ptr, uniform double* uniform y2ptr,
                               uniform double* uniform z1ptr, uniform double* uniform z2ptr,
                               uniform double* uniform fx1ptr, uniform double* uniform fx2ptr,
                               uniform double* uniform fy1ptr, uniform double* uniform fy2ptr,
                               uniform double* uniform fz1ptr, uniform double* uniform fz2ptr, uniform bool newton3) {
  for(uniform int i = 0; i < numParticles1; ++i) {
    if(ownedStatePtr1[i] == 0) {
      continue;
    }
    double fxacc = 0;
    double fyacc = 0;
    double fzacc = 0;
    foreach (j = 0 ... numParticles2) {
      // A double mask that is later multiplied to fac seems to be considerably faster than if or cif and a
      // little bit faster than adding this condition to the cutoff condition
      double mask = ownedStatePtr2[j] == 0 ? 0 : 1;

      double drx = x1ptr[i] - x2ptr[j];
      double dry = y1ptr[i] - y2ptr[j];
      double drz = z1ptr[i] - z2ptr[j];
      double dr2 = distanceSquared(drx, dry, drz);

      if (dr2 <= cutoffsquare) {
        double fac = mask * calculateLJFactor(dr2, sigmasquare, epsilon24);

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

inline void SoAFunctorVerletISPC(uniform size_t indexFirst, uniform size_t* uniform neighborList,
                                 uniform size_t numNeighbors,
                                 uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                                 uniform double* uniform x1ptr, uniform double* uniform y1ptr,
                                 uniform double* uniform z1ptr, uniform int64* uniform ownedStates,
                                 uniform double* uniform fx1ptr, uniform double* uniform fy1ptr,
                                 uniform double* uniform fz1ptr, uniform bool newton3) {
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

      if(newton3) {
        fx1ptr[index] -= fx;
        fy1ptr[index] -= fy;
        fz1ptr[index] -= fz;
      }

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
    }
  }
  fx1ptr[indexFirst] += reduce_add(fxacc);
  fy1ptr[indexFirst] += reduce_add(fyacc);
  fz1ptr[indexFirst] += reduce_add(fzacc);
}

inline void SoAFunctorSingleISPC(uniform int64 numParticles, uniform int64* uniform ownedStatePtr,
                                    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                                    uniform double* uniform xptr, uniform double* uniform yptr,
                                    uniform double* uniform zptr, uniform double* uniform fxptr,
                                    uniform double* uniform fyptr, uniform double* uniform fzptr, uniform bool newton3) {
  for(uniform size_t i = 0; i < numParticles; ++i) {
    if(ownedStatePtr[i] == 0) {
      continue;
    }

    double fxacc = 0;
    double fyacc = 0;
    double fzacc = 0;

    foreach (j = i + 1 ... numParticles) {
      if (ownedStatePtr[j] == 0) {
        continue;
      }
      double drx = xptr[i] - xptr[j];
      double dry = yptr[i] - yptr[j];
      double drz = zptr[i] - zptr[j];
      double dr2 = distanceSquared(drx, dry, drz);

      if (dr2 <= cutoffsquare) {
        double fac = calculateLJFactor(dr2, sigmasquare, epsilon24);

        double fx = drx * fac;
        double fy = dry * fac;
        double fz = drz * fac;

        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
      }
    }
    fxptr[i] += reduce_add(fxacc);
    fyptr[i] += reduce_add(fyacc);
    fzptr[i] += reduce_add(fzacc);
  }
}



// export declarations

export void SoAFunctorPairISPC_N3(uniform int64 numParticles1, uniform int64 numParticles2,
    uniform int64* uniform ownedStatePtr1, uniform int64* uniform ownedStatePtr2,
    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
    uniform double* uniform x1ptr, uniform double* uniform x2ptr,
    uniform double* uniform y1ptr, uniform double* uniform y2ptr,
    uniform double* uniform z1ptr, uniform double* uniform z2ptr,
    uniform double* uniform fx1ptr, uniform double* uniform fx2ptr,
    uniform double* uniform fy1ptr, uniform double* uniform fy2ptr,
    uniform double* uniform fz1ptr, uniform double* uniform fz2ptr) {
  SoAFunctorPairISPC(numParticles1, numParticles2, ownedStatePtr1, ownedStatePtr2, sigmasquare, epsilon24, cutoffsquare,
                     x1ptr, x2ptr, y1ptr, y2ptr, z1ptr, z2ptr, fx1ptr, fx2ptr, fy1ptr, fy2ptr, fz1ptr, fz2ptr, true);
}

export void SoAFunctorPairISPC_NoN3(uniform int64 numParticles1, uniform int64 numParticles2,
    uniform int64* uniform ownedStatePtr1, uniform int64* uniform ownedStatePtr2,
    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
    uniform double* uniform x1ptr, uniform double* uniform x2ptr,
    uniform double* uniform y1ptr, uniform double* uniform y2ptr,
    uniform double* uniform z1ptr, uniform double* uniform z2ptr,
    uniform double* uniform fx1ptr, uniform double* uniform fx2ptr,
    uniform double* uniform fy1ptr, uniform double* uniform fy2ptr,
    uniform double* uniform fz1ptr, uniform double* uniform fz2ptr) {
  SoAFunctorPairISPC(numParticles1, numParticles2, ownedStatePtr1, ownedStatePtr2, sigmasquare, epsilon24, cutoffsquare,
                     x1ptr, x2ptr, y1ptr, y2ptr, z1ptr, z2ptr, fx1ptr, fx2ptr, fy1ptr, fy2ptr, fz1ptr, fz2ptr, false);
}

export void SoAFunctorVerletISPC_N3(uniform size_t indexFirst, uniform size_t* uniform neighborList,
                                    uniform size_t numNeighbors, uniform double sigmasquare, uniform double epsilon24,
                                    uniform double cutoffsquare, uniform double* uniform x1ptr,
                                    uniform double* uniform y1ptr, uniform double* uniform z1ptr,
                                    uniform int64* uniform ownedStates, uniform double* uniform fx1ptr,
                                    uniform double* uniform fy1ptr, uniform double* uniform fz1ptr) {
  SoAFunctorVerletISPC(indexFirst, neighborList, numNeighbors, sigmasquare, epsilon24, cutoffsquare, x1ptr, y1ptr,
                       z1ptr, ownedStates, fx1ptr, fy1ptr, fz1ptr, true);
}

export void SoAFunctorVerletISPC_NoN3(uniform size_t indexFirst, uniform size_t* uniform neighborList,
                                    uniform size_t numNeighbors, uniform double sigmasquare, uniform double epsilon24,
                                    uniform double cutoffsquare, uniform double* uniform x1ptr,
                                    uniform double* uniform y1ptr, uniform double* uniform z1ptr,
                                    uniform int64* uniform ownedStates, uniform double* uniform fx1ptr,
                                    uniform double* uniform fy1ptr, uniform double* uniform fz1ptr) {
  SoAFunctorVerletISPC(indexFirst, neighborList, numNeighbors, sigmasquare, epsilon24, cutoffsquare, x1ptr, y1ptr,
                       z1ptr, ownedStates, fx1ptr, fy1ptr, fz1ptr, false);
}

export void SoAFunctorSingleISPC_N3(uniform int64 numParticles, uniform int64* uniform ownedStatePtr,
                                    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                                    uniform double* uniform xptr, uniform double* uniform yptr,
                                    uniform double* uniform zptr, uniform double* uniform fxptr,
                                    uniform double* uniform fyptr, uniform double* uniform fzptr) {
  SoAFunctorSingleISPC(numParticles, ownedStatePtr, sigmasquare, epsilon24, cutoffsquare, xptr, yptr, zptr,
                       fxptr, fyptr, fzptr, true);
}

export void SoAFunctorSingleISPC_NoN3(uniform int64 numParticles, uniform int64* uniform ownedStatePtr,
                                    uniform double sigmasquare, uniform double epsilon24, uniform double cutoffsquare,
                                    uniform double* uniform xptr, uniform double* uniform yptr,
                                    uniform double* uniform zptr, uniform double* uniform fxptr,
                                    uniform double* uniform fyptr, uniform double* uniform fzptr) {
  SoAFunctorSingleISPC(numParticles, ownedStatePtr, sigmasquare, epsilon24, cutoffsquare, xptr, yptr, zptr,
                       fxptr, fyptr, fzptr, false);
}