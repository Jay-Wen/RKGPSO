#include <windows.h>
#include <iostream>
#include <string>
#include <filesystem>
#include "RailwayEngineApi.h"
#include "RKGApi.h"
namespace fs = std::filesystem;

struct ConsolePause {
    bool enabled;
    ~ConsolePause() {
        if(enabled) { std::cout << "\nPress Enter to close this window..." << std::flush; std::cin.get(); }
    }
};
static bool ConfigExists(const fs::path& path) {
    std::error_code error;
    return fs::is_regular_file(path,error);
}
int wmain(int argc,wchar_t** argv) {
    SetConsoleOutputCP(CP_UTF8);
    ConsolePause pause{argc==1};
    wchar_t exe[32768];
    DWORD length=GetModuleFileNameW(NULL,exe,32768);
    if(!length || length>=32768) {std::cerr<<"Unable to determine the program directory.\n";return 2;}
    const fs::path dir=fs::path(exe).parent_path();
    const bool check=argc==2 && std::wstring(argv[1])==L"--check";
    fs::path config;
    if(argc==1) {
        
        config=dir/L"run.ini";
        if(!ConfigExists(config)) config=dir.parent_path()/L"run.ini";
    } else if(!check) {
       
        if(argc==2 && argv[1][0]!=L'-') config=argv[1];
        else if(argc==3 && std::wstring(argv[1])==L"--config") config=argv[2];
        else {std::cerr<<"Usage: double-click the EXE, or run rkgpso.exe <run.ini path>, rkgpso.exe --config <run.ini path>, or rkgpso.exe --check\n";return 2;}
    }
    if(!check) {
        if(!ConfigExists(config)) {
            std::cerr<<"run.ini was not found. Place it beside the EXE or one level above the bin directory, or specify it on the command line.\nSearched path: "<<config.u8string()<<"\n";return 2;
        }
        std::error_code error;
        config=fs::absolute(config,error);
        if(error){std::cerr<<"Unable to resolve the configuration path.\n";return 2;}
        std::cout<<"Configuration file: "<<config.u8string()<<"\n";
    }
    if(argc==1) {
        
        std::wstring command=L"\""+std::wstring(exe)+L"\" --config \""+config.wstring()+L"\"";
        STARTUPINFOW startup={};startup.cb=sizeof(startup);
        PROCESS_INFORMATION process={};
        std::cout<<"Starting calculation. This window will remain open when the run finishes.\n"<<std::flush;
        if(!CreateProcessW(exe,&command[0],NULL,NULL,TRUE,0,NULL,NULL,&startup,&process)) {
            std::cerr<<"Unable to start the calculation process. Windows error code="<<GetLastError()<<"\n";return 3;
        }
        CloseHandle(process.hThread);
        DWORD waited=WaitForSingleObject(process.hProcess,INFINITE),code=0;
        if(waited!=WAIT_OBJECT_0 || !GetExitCodeProcess(process.hProcess,&code)) {
            std::cerr<<"Unable to read the calculation process exit status. Windows error code="<<GetLastError()<<"\n";
            CloseHandle(process.hProcess);return 3;
        }
        CloseHandle(process.hProcess);
        std::cout<<"The calculation process has finished. Exit code="<<code<<" (0x"<<std::hex<<code<<std::dec<<")\n";
        if(code!=0)std::cout<<"The run did not complete successfully. Check rkgpso_console.log and the generated rkgpso_crash.txt, if present, in the project RESULT directory.\n";
        return static_cast<int>(code);
    }
    
    SetDllDirectoryW(dir.c_str());
    HMODULE library=LoadLibraryW((dir/L"RailwayEngine.dll").c_str());
    if(!library){std::cerr<<"Failed to load RailwayEngine.dll. Windows error code="<<GetLastError()<<"\n";return 3;}
    auto version=(RailwayEngineVersion)GetProcAddress(library,"RE_Version");
    auto engineLicense=(RailwayEngineLicenseStatus)GetProcAddress(library,"RE_LicenseStatus");
    auto engineExpiry=(RailwayEngineExpiryDate)GetProcAddress(library,"RE_ExpiryDate");
    auto run=(RailwayEngineRun)GetProcAddress(library,"RE_Run");
    if(!version||!engineLicense||!engineExpiry||!run||version()!=1||RKG_Version()!=1){std::cerr<<"DLL ABI version mismatch.\n";return 3;}
    if(engineLicense()!=0) {
        const unsigned date=engineExpiry();
        std::cerr<<"RailwayEngine.dll has expired. Last valid date: "
            <<date/10000<<"-"<<(date/100)%100<<"-"<<date%100<<".\n";
        return 5;
    }
    if(check) {std::cout<<"RKGAlgorithm ABI=1; RailwayEngine ABI=1; DLL load and license OK\n";return 0;}
    int result=run(config.c_str());
    std::cout<<"RKG-BOPSO exit code: "<<result<<"\n";
   
    return result;
}
