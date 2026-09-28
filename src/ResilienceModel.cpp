#define RKG_BUILD
#include "ResilienceApi.h"
#include "RKGApi.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <windows.h>
using namespace std;

static double nanValue(){return numeric_limits<double>::quiet_NaN();}

static bool licenseAvailable(){return true;}
static void invalidateProbabilities(RSSegment* s){
    if(s) s->p_0=s->p_1=s->p_2=s->p_3=nanValue();
}
static double clamp01(double x){return (std::max)(0.0,(std::min)(1.0,x));}
unsigned RS_Version(){return 2;}
const char* RS_ModelStatus(){return "paper_partial_alignment_legacy_accessibility_recovery_integral";}
double RS_SampleOffset(){return licenseAvailable()?10:nanValue();}
double RS_SampleStep(){return licenseAvailable()?20:nanValue();}

double RS_EffectivePGA(int type,double pga){return licenseAvailable()?(pga<0.05?(type==3?0.5:0.2):pga):nanValue();}

double RS_Lognormal(double pga,double median,double beta){
    if(!licenseAvailable())return nanValue();
    if(!std::isfinite(pga)||pga<0||median<=0||beta<=0)return nanValue();
    if(pga==0)return 0;
    return .5*std::erfc(-(std::log(pga)-std::log(median))/(beta*std::sqrt(2.0)));
}
static void probabilities(RSSegment& s,double q1,double q2,double q3){
    if(!std::isfinite(q1)||!std::isfinite(q2)||!std::isfinite(q3)){s.p_0=s.p_1=s.p_2=s.p_3=nanValue();return;}
    
    q1=clamp01(q1);q2=(std::min)(q1,clamp01(q2));q3=(std::min)(q2,clamp01(q3));
    s.p_0=1-q1;s.p_1=q1-q2;s.p_2=q2-q3;s.p_3=q3;
}
void RS_BridgeFragility(RSSegment* segment){
    if(!licenseAvailable()){invalidateProbabilities(segment);return;}
    RSSegment& s=*segment;double p=s.m_ctrl.PGA;
    if(s.bridgeHeight<=50) probabilities(s,RS_Lognormal(p,.20,.42),RS_Lognormal(p,.26,.40),RS_Lognormal(p,.38,.37));
    else if(s.bridgeHeight<=100) probabilities(s,RS_Lognormal(p,.31,.58),RS_Lognormal(p,.39,.55),RS_Lognormal(p,.50,.52));
    else {p*=1.0125;probabilities(s,RS_Lognormal(p,.60,.60),RS_Lognormal(p,.90,.60),RS_Lognormal(p,1.10,.60));}

}
void RS_Fragility(RSSegment* segment){
    if(!licenseAvailable()){invalidateProbabilities(segment);return;}
    RSSegment& s=*segment;
    if(s.Type==1){RS_BridgeFragility(&s);return;}
   
    s.p_1=0;double h=s.m_dem.DesH-s.m_dem.EarthH,p=s.m_ctrl.PGA*9.8;
    if(s.Type==2) RS_TunnelFragility(p,h,&s.p_0,&s.p_2,&s.p_3);
    else if(s.Type==3 && h>=0) RS_FillFragility(p,h,&s.p_0,&s.p_2,&s.p_3);
    else if(s.Type==3) RS_CutFragility(p,h,&s.p_0,&s.p_2,&s.p_3);
    else s.p_0=s.p_1=s.p_2=s.p_3=nanValue();
}



void RS_TunnelFragility(double PGA,double dertaH,double* none,double* moderate,double* extensive)
{
    if(!licenseAvailable()){
        if(none)*none=nanValue();if(moderate)*moderate=nanValue();if(extensive)*extensive=nanValue();return;
    }
    double& NonePro=*none; double& ModeratePro=*moderate; double& ExtensivePro=*extensive;

    if (dertaH >= -80)
    {
        
        if (PGA < 0.29 * 9.8)
        {
            NonePro = 1;
            ModeratePro = ExtensivePro = 0;
        }
        
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = 0.91;
            ModeratePro = 0.064;
            ExtensivePro = 0.026;

        }
       
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = 0.686;
            ModeratePro = 0.174;
            ExtensivePro = 0.14;
        }
    }
    else if (dertaH < -80 && dertaH >= -200)
    {
        double Slope = (-dertaH - 80) / 120;
       
        if (PGA < 0.29 * 9.8)
        {
            NonePro = 1;
            ModeratePro = ExtensivePro = 0;
        }
        
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = (0.938 - 0.91) * Slope + 0.91;
            ModeratePro = (0.062 - 0.064) * Slope + 0.064;
            ExtensivePro = -0.026 * Slope + 0.026;
        }
        
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = (0.739 - 0.686) * Slope + 0.686;
            ModeratePro = (0.221 - 0.174) * Slope + 0.174;
            ExtensivePro = (0.04 - 0.14) * Slope + 0.14;
        }
    }
    else if (dertaH < -200 && dertaH >= -300)
    {
        double Slope = (-dertaH - 200) / 100;
       
        if (PGA < 0.29 * 9.8)
        {
            NonePro = 1;
            ModeratePro = ExtensivePro = 0;
        }
       
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = (0.98 - 0.938) * Slope + 0.938;
            ModeratePro = (0.02 - 0.062) * Slope + 0.062;
            ExtensivePro = 0;
        }
        
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = (0.968 - 0.739) * Slope + 0.739;
            ModeratePro = (0.024 - 0.221) * Slope + 0.221;
            ExtensivePro = (0.008 - 0.04) * Slope + 0.04;
        }
    }
    else
    {
        double Slope = (-dertaH - 300) / 100;
        
        if (PGA < 0.3 * 9.8)
        {
            NonePro = 1;
            ModeratePro = ExtensivePro = 0;
        }
        
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = (1 - 0.98) * Slope + 0.98;
            ModeratePro = -0.02 * Slope + 0.02;
            ExtensivePro = 0;
        }
      
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = (0.979 - 0.968) * Slope + 0.968;
            ModeratePro = (0.018 - 0.024) * Slope + 0.024;
            ExtensivePro = (0.003 - 0.008) * Slope + 0.008;
        }
    }
}
void RS_FillFragility(double PGA,double dertaH,double* none,double* moderate,double* extensive)
{
    if(!licenseAvailable()){
        if(none)*none=nanValue();if(moderate)*moderate=nanValue();if(extensive)*extensive=nanValue();return;
    }
    double& NonePro=*none; double& ModeratePro=*moderate; double& ExtensivePro=*extensive;

    if (dertaH < 4 && dertaH >= 2)
    {
        double Slope = (dertaH - 2) / 2;
        
        if (PGA < 0.29 * 9.8)
        {
            NonePro = (0.949 - 0.967) * Slope + 0.967;
            ModeratePro = (0.038 - 0.024) * Slope + 0.024;
            ExtensivePro = (0.013 - 0.009) * Slope + 0.009;
        }
  
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = (0.882 - 0.915) * Slope + 0.915;
            ModeratePro = (0.075 - 0.052) * Slope + 0.052;
            ExtensivePro = (0.043 - 0.033) * Slope + 0.033;
        }
        /*else*/
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = (0.808 - 0.852) * Slope + 0.852;
            ModeratePro = (0.113 - 0.085) * Slope + 0.085;
            ExtensivePro = (0.079 - 0.063) * Slope + 0.063;
        }
    }
    
    else 
    {
        double Slope = (dertaH - 4) / 2;
        
        if (PGA < 0.29 * 9.8)
        {
            NonePro = (0.934 - 0.949) * Slope + 0.949;
            ModeratePro = (0.043 - 0.038) * Slope + 0.038;
            ExtensivePro = (0.023 - 0.013) * Slope + 0.013;
        }
        
        else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
        {
            NonePro = (0.853 - 0.882) * Slope + 0.882;
            ModeratePro = (0.095 - 0.075) * Slope + 0.075;
            ExtensivePro = (0.052 - 0.043) * Slope + 0.043;
        }
       
        else if (PGA >= 0.4 * 9.8)
        {
            NonePro = (0.776 - 0.808) * Slope + 0.808;
            ModeratePro = (0.132 - 0.113) * Slope + 0.113;
            ExtensivePro = (0.092 - 0.079) * Slope + 0.079;
        }
    }
}
void RS_CutFragility(double PGA,double dertaH,double* none,double* moderate,double* extensive)
{
    if(!licenseAvailable()){
        if(none)*none=nanValue();if(moderate)*moderate=nanValue();if(extensive)*extensive=nanValue();return;
    }
    double& NonePro=*none; double& ModeratePro=*moderate; double& ExtensivePro=*extensive;

    double Slope = (dertaH + 2) / (-6 + 2);
    
    if (PGA < 0.29 * 9.8)
    {
        NonePro = (0.903 - 1) * Slope + 1;
        ModeratePro = (0.068 - 0) * Slope;
        ExtensivePro = (0.029 - 0) * Slope;
    }
   
    else if (PGA >= 0.29 * 9.8 && PGA < 0.4 * 9.8)
    {
        NonePro = (0.812 - 1) * Slope + 1;
        ModeratePro = (0.118 - 0) * Slope;
        ExtensivePro = (0.07 - 0) * Slope;
    }
    
    else if (PGA >= 0.4 * 9.8)
    {
        NonePro = (0.731 - 1) * Slope + 1;
        ModeratePro = (0.151 - 0) * Slope;
        ExtensivePro = (0.118 - 0) * Slope;
    }
}

namespace {
using RecoveryCore = void (__cdecl *)(RSSegment*);

RecoveryCore loadRecoveryCore()
{
    
    HMODULE engine = GetModuleHandleW(L"RailwayEngine.dll");
    if (!engine) {
        static int moduleAnchor = 0;
        HMODULE self = nullptr;
        wchar_t path[MAX_PATH] = {};
        GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCWSTR>(&moduleAnchor), &self);
        if (self && GetModuleFileNameW(self, path, MAX_PATH)) {
            wchar_t* slash = wcsrchr(path, L'\\');
            if (slash) wcscpy_s(slash + 1, MAX_PATH - (slash + 1 - path), L"RailwayEngine.dll");
            engine = LoadLibraryExW(path, nullptr,
                LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
            if (!engine) engine = LoadLibraryW(path);
        }
    }
    return engine ? reinterpret_cast<RecoveryCore>(GetProcAddress(engine, "RE_RecoveryCore")) : nullptr;
}
}

void RS_Recovery(RSSegment* segment)
{
    if(!licenseAvailable()){
        if(segment)segment->T_1=segment->T_2=segment->T_3=nanValue();return;
    }
    static RecoveryCore core = loadRecoveryCore();
    if (core) {
        core(segment);
        return;
    }
   
    if (segment) segment->T_1 = segment->T_2 = segment->T_3 = nanValue();
}
void RS_Remaining(RSSegment* segment){
    if(!licenseAvailable()){
        if(segment)segment->fun_0=segment->fun_1=segment->fun_2=segment->fun_3=nanValue();return;
    }
    RSSegment& s=*segment;s.fun_0=1;s.fun_1=s.Type==2?.86:.83;
    s.fun_2=(s.Type==1?.37:s.Type==2?.40:.34)*.36;s.fun_3=0;
}


double RS_SegmentResilience(RSSegment* segment)
{
if(!licenseAvailable())return nanValue();
RSSegment& temp=*segment;
    double r_seg = 0;
    if (temp.Type == 1) {
		temp.r_1 = log((1.0 - temp.fun_1) / 0.2) / temp.T_1;
        temp.r_2 = log((1.0 - temp.fun_2) / 0.2) / temp.T_2;
        temp.r_3 = log((1.0 - temp.fun_3) / 0.2) / temp.T_3;

        double T = 180;
        r_seg = temp.p_0 * 1 + temp.p_1 * (1 - (1 - temp.fun_1) * (1 - exp(-temp.r_1 * temp.T_1)) / (temp.r_1 * T))
            + temp.p_2 * (1 - (1 - temp.fun_2) * (1 - exp(-temp.r_2 * temp.T_2)) / (temp.r_2 * T))
            + temp.p_3 * (1 - (1 - temp.fun_3) * (1 - exp(-temp.r_3 * temp.T_3)) / (temp.r_3 * T));
    }
    else{
        temp.r_1 = 0;
        temp.r_2 = log((1.0 - temp.fun_2) / 0.2) / temp.T_2;
        temp.r_3 = log((1.0 - temp.fun_3) / 0.2) / temp.T_3;
        double T = 180;
        r_seg = temp.p_0 * 1
            + temp.p_2 * (1 - (1 - (temp.fun_2 + temp.fun_1)/2) * (1 - exp(-temp.r_2 * temp.T_2)) / (temp.r_2 * T))
            + temp.p_3 * (1 - (1 - temp.fun_3) * (1 - exp(-temp.r_3 * temp.T_3)) / (temp.r_3 * T));
    }
    

return r_seg;
}
double RS_LineResilience(const RSPoint* points,int n){
    if(!licenseAvailable())return nanValue();
    if(!points||n<1)return nanValue();
  
    double sum=0,weakest=numeric_limits<double>::infinity();
    for(int i=0;i<n;++i){double v=points[i].r;if(!std::isfinite(v)||v<0)return nanValue();if(v==0){weakest=0;continue;}sum+=1/v;weakest=(std::min)(weakest,v);}
    
    if(weakest==0)return 0;
    return 0.2*n/sum+(1-0.2)*weakest;
}
double RS_DesignHeight(const RSProfilePoint* points,int n,double cml)
{
    if(!licenseAvailable())return nanValue();

    if (n <= 0)
        return 0.0;

    if (n == 1)
        return points[0].height;

    
    if (cml <= points[0].mileage)
        return points[0].height;

    if (cml >= points[n - 1].mileage)
        return points[n - 1].height;

    for (int i = 0; i < n - 1; ++i)
    {
        double c1 = points[i].mileage;
        double c2 = points[i + 1].mileage;

        if (cml >= c1 && cml <= c2)
        {
            double h1 = points[i].height;
            double h2 = points[i + 1].height;

            if (fabs(c2 - c1) < 1e-8)
                return h1;

            double ratio = (cml - c1) / (c2 - c1);
            return h1 + ratio * (h2 - h1);
        }
    }

    return points[n - 1].height;
}
static vector<RSRange> MergeRanges(vector<RSRange> ranges)
{
    vector<RSRange> merged;
    if (ranges.empty()) return merged;

    sort(ranges.begin(), ranges.end(), [](const RSRange& a,const RSRange& b){return a.St<b.St;});
    merged.push_back(ranges[0]);

    for (size_t i = 1; i < ranges.size(); ++i)
    {
        if (ranges[i].St <= merged.back().Ed)
        {
            merged.back().Ed = max(merged.back().Ed, ranges[i].Ed);
        }
        else
        {
            merged.push_back(ranges[i]);
        }
    }
    return merged;
}
int RS_RoadRanges(const RSRange* input,int count,double lineLength,RSRange* out,int capacity)
{
if(!licenseAvailable())return -2;
vector<RSRange> occupied; if(count>0) occupied.assign(input,input+count);
   
    vector<RSRange> merged = MergeRanges(occupied);

    vector<RSRange> subgrade;
    double lineSt = 0.0;
    double lineEd = lineLength;

    double cur = lineSt;

    
    for (size_t i = 0; i < merged.size(); ++i)
    {
        if (cur < merged[i].St)
        {
            RSRange r;
            r.St = cur;
            r.Ed = merged[i].St;
            subgrade.push_back(r);

        }
        cur = max(cur, merged[i].Ed);
    }
    if (cur < lineEd) {
        RSRange r;
        r.St = cur;
        r.Ed = lineEd;
        subgrade.push_back(r);
    }


    if(!out) return static_cast<int>(subgrade.size());
    if(capacity<static_cast<int>(subgrade.size())) return -1;
    std::copy(subgrade.begin(),subgrade.end(),out);
    return static_cast<int>(subgrade.size());
}
