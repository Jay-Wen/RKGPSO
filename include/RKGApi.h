#pragma once

#ifdef RKG_BUILD
#define RKG_API extern "C" __declspec(dllexport)
#else
#define RKG_API extern "C" __declspec(dllimport)
#endif
struct RKGConfig {
    int multiScale, temporalGate, kneeLeader, kneeAttraction, minimumArchive;
    double c3, gamma, confidenceThreshold, stabilityScale, kneeFraction, epsilon;
};

RKG_API unsigned RKG_Version();

RKG_API int RKG_Archive(const double* objectives,int n,int capacity,int* indices);

RKG_API int RKG_UpdateArchive(const double* archiveObjectives,int archiveCount,
    const double* candidateObjectives,int candidateCount,int capacity,
    int* sourceKinds,int* sourceIndices,double* outputObjectives);
RKG_API int RKG_Detect(const double* objectives,int n,const RKGConfig* config,int hadPrevious,
    double* previousAndCentroid,double* scores,int* indices,double* confidence,int* gate);
RKG_API int RKG_Weights(const double* objectives,int n,double* weights);
RKG_API int RKG_Sample(const double* weights,int n,double uniform);
RKG_API double RKG_Attraction(const RKGConfig* config,int gate,double confidence,double progress);
RKG_API double RKG_Velocity(double inertia,double oldVelocity,double c1,double r1,double personal,
    double c2,double r2,double leader,double c3,double r3,double knee,double current,double limit);

RKG_API int RKG_UpdateParticleState(int dimensions,double inertia,const double* oldVelocity,
    double c1,const double* r1,const double* personal,double c2,const double* r2,const double* leader,
    double c3,const double* r3,const double* knee,const double* current,const double* limits,
    double* nextVelocity,double* nextPosition);

typedef int (__cdecl *RKGStep)(void* context,int iteration,double inertia);
RKG_API int RKG_Run(int generations,void* context,RKGStep step);

typedef double (__cdecl *RKGUniform)(void* context);
RKG_API int RKG_Leaders(const double* objectives,int n,const double* scores,const int* kneeIds,int kneeCount,
    int gate,int useKneeLeader,double epsilon,int population,int* leaders,int* kneeLeaders,void* context,RKGUniform random);

RKG_API int RKG_UpdatePersonal(const double* archive,int n,const double* current,const double* personal,
    double epsilon,void* context,RKGUniform random);
