#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <string>

namespace rkg {
struct Config {
    bool enabled;             
    bool multiScale;          
    bool temporalGate;        
    bool kneeLeader;         
    bool kneeAttraction;      
    int archiveCapacity;      
    int minimumArchive;       
    int plotEvery;            
    double c3;                
    double gamma;             
    double confidenceThreshold;
    double stabilityScale;    
    double kneeFraction;      
    double epsilon;           
    Config() : enabled(true), multiScale(true), temporalGate(true), kneeLeader(true),
        kneeAttraction(true), archiveCapacity(100), minimumArchive(7), plotEvery(10),
        c3(1.0), gamma(2.0), confidenceThreshold(0.8), stabilityScale(0.1),
        kneeFraction(0.2), epsilon(1e-12) {}
};

inline std::string variantTag(const Config& c) {
    std::string tag;
    if(!c.multiScale) tag+="NoMultiScale";
    if(!c.temporalGate) {if(!tag.empty()) tag+="_";tag+="NoTemporalGate";}
    if(!c.kneeLeader) {if(!tag.empty()) tag+="_";tag+="NoKneeLeader";}
    if(!c.kneeAttraction) {if(!tag.empty()) tag+="_";tag+="NoKneeAttraction";}
    return tag.empty()?"Full":tag;
}
struct Point {
    double cost, loss; 
    Point(double c = 0, double l = 0) : cost(c), loss(l) {}
};
inline double clamp(double x, double lo, double hi) {
    return (std::max)(lo, (std::min)(hi, x));
}
inline bool finite(double x) {
    return x == x && x <= (std::numeric_limits<double>::max)()
        && x >= -(std::numeric_limits<double>::max)();
}
inline bool valid(const Point& p) {
    return finite(p.cost) && finite(p.loss) && p.cost >= 0 && p.loss >= 0 && p.loss <= 1;
}
inline bool dominates(const Point& a, const Point& b) {
    
    return a.cost <= b.cost && a.loss <= b.loss && (a.cost < b.cost || a.loss < b.loss);
}
struct Bounds {
    double c0, c1, l0, l1;
    Bounds() : c0(0), c1(0), l0(0), l1(0) {}
    explicit Bounds(const std::vector<Point>& p) : c0(0), c1(0), l0(0), l1(0) {
        if (p.empty()) return;
        c0 = c1 = p[0].cost; l0 = l1 = p[0].loss;
        for (size_t i = 1; i < p.size(); ++i) {
            c0 = (std::min)(c0, p[i].cost); c1 = (std::max)(c1, p[i].cost);
            l0 = (std::min)(l0, p[i].loss); l1 = (std::max)(l1, p[i].loss);
        }
    }
    Point normalize(const Point& p, double eps = 1e-12) const {
        return Point((p.cost-c0)/(std::max)(eps,c1-c0),
                     (p.loss-l0)/(std::max)(eps,l1-l0));
    }
};

inline std::vector<int> nonDominated(const std::vector<Point>& p) {
    std::vector<int> ids;
    for (size_t i=0; i<p.size(); ++i) {
        if (!valid(p[i])) continue;
        bool keep = true;
        for (size_t j=0; j<p.size(); ++j) {
            if (j==i || !valid(p[j])) continue;
            if (dominates(p[j],p[i]) || (j<i && p[j].cost==p[i].cost && p[j].loss==p[i].loss)) {
                keep=false; break;
            }
        }
        if (keep) ids.push_back(static_cast<int>(i));
    }
    std::sort(ids.begin(),ids.end(),[&](int a,int b) {
        return p[a].cost < p[b].cost || (p[a].cost==p[b].cost && p[a].loss<p[b].loss);
    });
    return ids;
}

inline std::vector<double> crowding(const std::vector<Point>& p) {
    std::vector<double> d(p.size(),0);
    if (p.empty()) return d;
    for (int objective=0; objective<2; ++objective) {
        std::vector<int> ids;
        for (size_t i=0;i<p.size();++i) ids.push_back(static_cast<int>(i));
        std::sort(ids.begin(),ids.end(),[&](int a,int b) {
            double x=objective==0?p[a].cost:p[a].loss;
            double y=objective==0?p[b].cost:p[b].loss;
            return x==y?a<b:x<y;
        });
        d[ids.front()]=d[ids.back()]=(std::numeric_limits<double>::infinity)();
        double range=objective==0?p[ids.back()].cost-p[ids.front()].cost:
            p[ids.back()].loss-p[ids.front()].loss;
        if (range<=0) continue;
        for (size_t j=1;j+1<ids.size();++j)
            d[ids[j]]+=(objective==0?p[ids[j+1]].cost-p[ids[j-1]].cost:
                p[ids[j+1]].loss-p[ids[j-1]].loss)/range;
    }
    return d;
}
inline std::vector<int> archiveIndices(const std::vector<Point>& p,int capacity) {
    std::vector<int> ids=nonDominated(p);
    capacity=(std::max)(2,capacity); 
    while (ids.size()>static_cast<size_t>(capacity)) {
        std::vector<Point> q;
        for (size_t i=0;i<ids.size();++i) q.push_back(p[ids[i]]);
        std::vector<double> d=crowding(q);
        size_t remove=1;
        for (size_t i=1;i+1<ids.size();++i) if(d[i]<d[remove]) remove=i;
        ids.erase(ids.begin()+remove); 
    }
    return ids;
}
struct Knee {
    std::vector<double> score;
    std::vector<int> indices;
    Point centroid; 
    double confidence;
    bool gate;
    Knee() : confidence(0),gate(false) {}
};

inline Knee detect(const std::vector<Point>& p,const Config& cfg,
                   bool hadPrevious,const Point& previous) {
    Knee k; k.score.assign(p.size(),0);
    if (p.size()<3) return k;
    Bounds bounds(p);
    if(bounds.c1-bounds.c0<=cfg.epsilon || bounds.l1-bounds.l0<=cfg.epsilon) return k;
    std::vector<int> ids;
    for(size_t i=0;i<p.size();++i) ids.push_back(static_cast<int>(i));
    std::sort(ids.begin(),ids.end(),[&](int a,int b){return p[a].cost<p[b].cost;});
    std::vector<int> scales; scales.push_back(1);
    if(cfg.multiScale) { scales.push_back(2); scales.push_back(3); }
    for(int i=1;i+1<static_cast<int>(ids.size());++i) {
        std::vector<double> distances;
        for(size_t s=0;s<scales.size();++s) {
            int h=scales[s];
            if(i-h<0 || i+h>=static_cast<int>(ids.size())) continue;
            Point a=bounds.normalize(p[ids[i-h]],cfg.epsilon);
            Point b=bounds.normalize(p[ids[i+h]],cfg.epsilon);
            Point q=bounds.normalize(p[ids[i]],cfg.epsilon);
            double dx=b.cost-a.cost,dy=b.loss-a.loss;
            double length=std::sqrt(dx*dx+dy*dy);
            if(length<=cfg.epsilon) continue;
           
            double distance=std::fabs(dx*(q.loss-a.loss)-dy*(q.cost-a.cost))/length;
            distances.push_back(distance/(length+cfg.epsilon));
        }
        if(distances.empty()) continue;
        std::sort(distances.begin(),distances.end());
        size_t n=distances.size();
        k.score[ids[i]]=n%2?distances[n/2]:0.5*(distances[n/2-1]+distances[n/2]);
    }
    std::vector<int> order=ids;
    std::sort(order.begin(),order.end(),[&](int a,int b) {
        return k.score[a]==k.score[b]?a<b:k.score[a]>k.score[b];
    });
    int count=(std::max)(1,static_cast<int>(std::ceil(p.size()*cfg.kneeFraction)));
    double sum=0;
    for(size_t i=0;i<order.size() && static_cast<int>(k.indices.size())<count;++i) {
        int j=order[i];
        if(k.score[j]<=cfg.epsilon) break;
        k.indices.push_back(j); sum+=k.score[j];
        k.centroid.cost+=k.score[j]*p[j].cost;
        k.centroid.loss+=k.score[j]*p[j].loss;
    }
    if(sum<=cfg.epsilon) return k;
    k.centroid.cost/=sum; k.centroid.loss/=sum;
    if(hadPrevious) {
        Point now=bounds.normalize(k.centroid,cfg.epsilon), old=bounds.normalize(previous,cfg.epsilon);
        double dx=now.cost-old.cost,dy=now.loss-old.loss;
        k.confidence=std::exp(-std::sqrt(dx*dx+dy*dy)/(std::max)(cfg.stabilityScale,cfg.epsilon));
    }
    k.gate=static_cast<int>(p.size())>=cfg.minimumArchive &&
        (!cfg.temporalGate || (hadPrevious && k.confidence>=cfg.confidenceThreshold));
    return k;
}

inline std::vector<double> crowdingWeights(const std::vector<Point>& p) {
    std::vector<double> w=crowding(p);
    double largest=1;
    for(size_t i=0;i<w.size();++i) if(finite(w[i])) largest=(std::max)(largest,w[i]);
    for(size_t i=0;i<w.size();++i) w[i]=finite(w[i])?(std::max)(w[i],1e-12):2*largest;
    return w;
}

inline int sample(const std::vector<double>& weights,double uniform) {
    double total=0;
    for(size_t i=0;i<weights.size();++i) total+=(std::max)(0.0,weights[i]);
    if(total<=0) return -1;
    double target=clamp(uniform,0,1-1e-15)*total,cumulative=0;
    int last=-1;
    for(size_t i=0;i<weights.size();++i) if(weights[i]>0) {
        cumulative+=weights[i]; last=static_cast<int>(i);
        if(target<cumulative) return last;
    }
    return last;
}
inline double attraction(const Config& cfg,const Knee& k,double progress) {
    if(!cfg.kneeAttraction || !k.gate) return 0;
    
    return cfg.c3*std::pow(clamp(progress,0,1),cfg.gamma)*(cfg.temporalGate?k.confidence:1.0);
}

inline double velocity(double inertia,double oldVelocity,double c1,double r1,double personal,
    double c2,double r2,double leader,double c3,double r3,double knee,double current,double limit) {
    return clamp(inertia*oldVelocity+c1*r1*(personal-current)+c2*r2*(leader-current)+
                 c3*r3*(knee-current),-limit,limit);
}
struct State {
    Config config;
    bool active, hadPrevious;
    int iteration;
    unsigned long long evaluations, attempts, accepted; 
    Point previous;
    Knee knee;
    std::vector<int> leaders,kneeLeaders; 
    std::vector<std::vector<double> > vx,vy,vr;
    double strength; 
    State() : active(false),hadPrevious(false),iteration(0),evaluations(0),attempts(0),accepted(0),strength(0) {}
};
} 
