#define RKG_BUILD
#include "RKGApi.h"
#include "RKGBOPSO.h"
#include <stdexcept>
namespace {
std::vector<rkg::Point> points(const double* values,int n) {
    if(n<0 || (n>0 && !values)) throw std::invalid_argument("objectives");
    std::vector<rkg::Point> p;for(int i=0;i<n;++i) p.push_back(rkg::Point(values[2*i],values[2*i+1]));return p;
}
rkg::Config config(const RKGConfig* value) {
    if(!value) throw std::invalid_argument("config");
    rkg::Config c;c.multiScale=value->multiScale!=0;c.temporalGate=value->temporalGate!=0;
    c.kneeLeader=value->kneeLeader!=0;c.kneeAttraction=value->kneeAttraction!=0;c.minimumArchive=value->minimumArchive;
    c.c3=value->c3;c.gamma=value->gamma;c.confidenceThreshold=value->confidenceThreshold;
    c.stabilityScale=value->stabilityScale;c.kneeFraction=value->kneeFraction;c.epsilon=value->epsilon;return c;
}
}
unsigned RKG_Version(){return 1;}
int RKG_Archive(const double* values,int n,int capacity,int* indices) {
    try {if(n>0&&!indices)return -1;auto ids=rkg::archiveIndices(points(values,n),capacity);
        std::copy(ids.begin(),ids.end(),indices);return (int)ids.size();}catch(...){return -1;}
}
int RKG_UpdateArchive(const double* archiveValues,int archiveCount,const double* candidateValues,
    int candidateCount,int capacity,int* sourceKinds,int* sourceIndices,double* outputValues) {
    try {
        if(archiveCount<0||candidateCount<0||capacity<2||
           (archiveCount>0&&!archiveValues)||(candidateCount>0&&!candidateValues)||
           !sourceKinds||!sourceIndices||!outputValues)return -1;
        std::vector<rkg::Point> merged=points(archiveValues,archiveCount);
        std::vector<rkg::Point> candidates=points(candidateValues,candidateCount);
        merged.insert(merged.end(),candidates.begin(),candidates.end());
        std::vector<int> keep=rkg::archiveIndices(merged,capacity);
        for(size_t i=0;i<keep.size();++i){
            int id=keep[i];
            sourceKinds[i]=id<archiveCount?0:1;
            sourceIndices[i]=id<archiveCount?id:id-archiveCount;
            outputValues[2*i]=merged[id].cost;
            outputValues[2*i+1]=merged[id].loss;
        }
        return static_cast<int>(keep.size());
    }catch(...){return -1;}
}
int RKG_Detect(const double* values,int n,const RKGConfig* c,int hadPrevious,double* center,
    double* scores,int* indices,double* confidence,int* gate) {
    try {if(!center||!confidence||!gate||(n>0&&(!scores||!indices)))return -1;
        auto k=rkg::detect(points(values,n),config(c),hadPrevious!=0,rkg::Point(center[0],center[1]));
        std::copy(k.score.begin(),k.score.end(),scores);std::copy(k.indices.begin(),k.indices.end(),indices);
        center[0]=k.centroid.cost;center[1]=k.centroid.loss;*confidence=k.confidence;*gate=k.gate?1:0;return (int)k.indices.size();
    }catch(...){return -1;}
}
int RKG_Weights(const double* values,int n,double* weights) {
    try {if(n>0&&!weights)return -1;auto w=rkg::crowdingWeights(points(values,n));std::copy(w.begin(),w.end(),weights);return n;}
    catch(...){return -1;}
}
int RKG_Sample(const double* values,int n,double uniform) {
    try {if(n<=0||!values)return -1;return rkg::sample(std::vector<double>(values,values+n),uniform);}catch(...){return -1;}
}
double RKG_Attraction(const RKGConfig* c,int gate,double confidence,double progress) {
    try {rkg::Knee k;k.gate=gate!=0;k.confidence=confidence;return rkg::attraction(config(c),k,progress);}
    catch(...){return std::numeric_limits<double>::quiet_NaN();}
}
double RKG_Velocity(double inertia,double oldVelocity,double c1,double r1,double personal,
    double c2,double r2,double leader,double c3,double r3,double knee,double current,double limit) {
    return rkg::velocity(inertia,oldVelocity,c1,r1,personal,c2,r2,leader,c3,r3,knee,current,limit);
}
int RKG_UpdateParticleState(int dimensions,double inertia,const double* oldVelocity,
    double c1,const double* r1,const double* personal,double c2,const double* r2,const double* leader,
    double c3,const double* r3,const double* knee,const double* current,const double* limits,
    double* nextVelocity,double* nextPosition) {
    try {
        if(dimensions<1||!oldVelocity||!r1||!personal||!r2||!leader||!r3||!knee||
           !current||!limits||!nextVelocity||!nextPosition)return -1;
        for(int i=0;i<dimensions;++i){
            if(!rkg::finite(limits[i])||limits[i]<0)return -1;
            nextVelocity[i]=rkg::velocity(inertia,oldVelocity[i],c1,r1[i],personal[i],
                c2,r2[i],leader[i],c3,r3[i],knee[i],current[i],limits[i]);
            nextPosition[i]=current[i]+nextVelocity[i];
        }
        return dimensions;
    }catch(...){return -1;}
}
int RKG_Run(int generations,void* context,RKGStep step) {
    if(generations<1||!step)return -1;
    try {
        
        for(int t=0;t<generations;++t){
            double u=double(t)/generations;
            double inertia=0.9-2.25*u*u*u/7.0+6.75*u*u/7.0-8.0*u/7.0;
            if(step(context,t,inertia)!=1)return t;
        }
        return generations;
    }catch(...){return -1;}
}
int RKG_Leaders(const double* values,int n,const double* scores,const int* kneeIds,int kneeCount,
    int gate,int useKneeLeader,double epsilon,int population,int* leaders,int* kneeLeaders,void* context,RKGUniform random) {
    try {
        if(n<1||!scores||kneeCount<0||population<1||!leaders||!kneeLeaders||!random||(kneeCount>0&&!kneeIds))return -1;
        auto weights=rkg::crowdingWeights(points(values,n));std::vector<double> knees(n,0);
        for(int i=0;i<kneeCount;++i){int j=kneeIds[i];if(j<0||j>=n)return -1;knees[j]=scores[j]+epsilon;}
        for(int i=0;i<population;++i){
            leaders[i]=rkg::sample(gate&&useKneeLeader?knees:weights,random(context));
            kneeLeaders[i]=gate?rkg::sample(knees,random(context)):leaders[i];
            if(leaders[i]<0||kneeLeaders[i]<0)return -1;
        }
        return population;
    }catch(...){return -1;}
}
int RKG_UpdatePersonal(const double* values,int n,const double* current,const double* personal,
    double epsilon,void* context,RKGUniform random) {
    try {
        if(n<1||!current||!personal||!random)return -1;
        auto archive=points(values,n);rkg::Point a(current[0],current[1]),b(personal[0],personal[1]);
        bool replace=rkg::valid(a)&&(!rkg::valid(b)||rkg::dominates(a,b));
        if(rkg::valid(a)&&rkg::valid(b)&&!rkg::dominates(a,b)&&!rkg::dominates(b,a)){
            rkg::Bounds bounds(archive);auto sparse=rkg::crowdingWeights(archive);
            int target=rkg::sample(sparse,random(context));if(target<0)return -1;
            auto q=bounds.normalize(archive[target]),x=bounds.normalize(a),y=bounds.normalize(b);
            double angle=std::atan2(q.loss,q.cost);
            double da=std::fabs(std::atan2(x.loss,x.cost)-angle),db=std::fabs(std::atan2(y.loss,y.cost)-angle);
            replace=da<db||(std::fabs(da-db)<epsilon&&random(context)<0.5);
        }
        return replace?1:0;
    }catch(...){return -1;}
}
