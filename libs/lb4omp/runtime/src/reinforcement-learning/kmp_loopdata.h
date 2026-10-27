#pragma once

typedef struct {
    int n;              // #iterations (for reinforcement learning agent)
    int p;              // #threads (for reinforcement learning agent)
    int autoSearch;     // flag: continue automatic search or not
    int cDLS;           // current DLS
    int bestDLS;        // Best DLS
    int searchTrials;   // number of trials to find the best DLS
    int cChunk;         // current chunk size
    int bestChunk;      // chunk size of the best DLS technique
    double cTime;       // current DLS time
    double bestTime;    // loop time of the best DLS
    double cLB;         // load imbalance of the current DLS
    double bestLB;      // load imbalance of the best DLS
    double cRobust;     // current Robustness
    double cStdDev;     // standard deviation of the current DLS
    double bestStdDev;  // standard deviation of the best DLS
    double cSkew;  // skewness of the current DLS
    double bestSkew;// skewness of the best DLS
    double cKurt;   // kurtosis of the current DLS
    double bestKurt;// kurtosis of the best DLS
    double cCOV;    // coefficient of variation of the current DLS
    double bestCOV; // coefficient of variation of the best DLS
} LoopData;
