#pragma once


typedef unsigned (__cdecl *RailwayEngineVersion)();

typedef int (__cdecl *RailwayEngineLicenseStatus)();

typedef unsigned (__cdecl *RailwayEngineExpiryDate)();
typedef int (__cdecl *RailwayEngineRun)(const wchar_t* configPath);
