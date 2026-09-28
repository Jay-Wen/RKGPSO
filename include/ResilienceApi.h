#pragma once

#ifdef RKG_BUILD
#define RS_API extern "C" __declspec(dllexport)
#else
#define RS_API extern "C" __declspec(dllimport)
#endif
struct RSRange { double St,Ed; };
struct RSPoint { double cml,r; };
struct RSProfilePoint { double mileage,height; };
struct RSSegment {
    int Type; 
    double bridgeHeight; 
    struct { double x,y,lc; } pt; 
    struct { double PGA; } m_ctrl; 
    struct { double DesH,EarthH; } m_dem; 
    double p_0,p_1,p_2,p_3;
    double fun_0,fun_1,fun_2,fun_3;
    double T_1,T_2,T_3,r_1,r_2,r_3;
};
RS_API unsigned RS_Version();

RS_API double RS_SampleOffset();
RS_API double RS_SampleStep();
RS_API double RS_EffectivePGA(int type,double pga);

RS_API int RS_RoadRanges(const RSRange* occupied,int count,double lineLength,RSRange* out,int capacity);
RS_API double RS_DesignHeight(const RSProfilePoint* points,int count,double mileage);
RS_API double RS_Lognormal(double pga,double average,double dispersion);
RS_API void RS_BridgeFragility(RSSegment* segment);
RS_API void RS_TunnelFragility(double pgaMS2,double height,double* none,double* moderate,double* extensive);
RS_API void RS_FillFragility(double pgaMS2,double height,double* none,double* moderate,double* extensive);
RS_API void RS_CutFragility(double pgaMS2,double height,double* none,double* moderate,double* extensive);
RS_API void RS_Recovery(RSSegment* segment);
RS_API void RS_Remaining(RSSegment* segment);
RS_API double RS_SegmentResilience(RSSegment* segment);
RS_API double RS_LineResilience(const RSPoint* points,int count);

RS_API void RS_Fragility(RSSegment* segment);
RS_API const char* RS_ModelStatus();
