#pragma once
#ifndef SIMHANDLER_H
#define SIMHANDLER_H

#include "pch.h"

#define DLL_API __declspec(dllexport)
//----------------------------------------------------------------------
extern "C" {
	DLL_API void    SIM_Initialize(void* (*callback)(int, const void*));
	DLL_API void    SIM_Shutdown();
	DLL_API int64_t SIM_HandleMessage(Hashes::SimMessageHashes sim_msg_id, int msg_length, char* msg);
	DLL_API char*   SIM_BeginSave(uint32_t* out_size, int x, int y);
	DLL_API void    SIM_EndSave();
	DLL_API void    SIM_DebugCrash();

	DLL_API void    ConduitTemperatureManager_Initialize();
	DLL_API void    ConduitTemperatureManager_Shutdown();
	DLL_API int     ConduitTemperatureManager_Add(float, float, int, int, float, float, bool);
	DLL_API void    ConduitTemperatureManager_Set(int, float, float, int);
	DLL_API void    ConduitTemperatureManager_Remove(int handle);
	DLL_API void*   ConduitTemperatureManager_Update(float dt, BuildingTemperatureInfo* building_conductivity_data);
	DLL_API void    ConduitTemperatureManager_Clear();
}
#endif