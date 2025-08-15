#pragma once
#ifndef CLASS_CONDUIT_H
#define CLASS_CONDUIT_H

#include "ClassSim.h"
#include "SimHandler.h"
#include <threads.h>

//----- define struct -----
#pragma pack(push, 1)
struct alignas(1) ConduitTemperatureUpdateData
{
	int    numEntries;
	float* temperatures;
	int    numFrozenHandles;
	int*   frozenHandles;
	int    numMeltedHandles;
	int*   meltedHandles;
};
#pragma pack(pop)
class ConduitTemperatureManager
{
public:
	struct Data
	{
		float  temperature;                    // Content
		float  thermalConductivity;            // Content
		float  heatCapacity;                   // Content
		Handle conduitBuildingTemperatureHandle;
		float  conduitHeatCapacity;            // Conduit
		float  conduitThermalConductivity;     // Conduit
		bool   conduitInsulated;               // Conduit
		float  lowStateTransitionTemperature;  // Content
		float  highStateTransitionTemperature; // Content
	};
	CompactedVector<ConduitTemperatureManager::Data> data;
	std::vector<int>   meltedContentHandles;
	std::vector<int>   frozenContentHandles;
	std::vector<float> temperatures;           // Content
	ConduitTemperatureUpdateData updateData;
	int curReleaseListIdx = 0;
	std::vector<Handle> frameDelayedReleasedHandles[2];
	mtx_t dataMutex;

	Handle Add(float, float, uint32_t, int, float, float, bool);
	void   ReleaseQueuedHandles();
	void   Set(int, float, float, uint32_t);
	ConduitTemperatureUpdateData* Update(float dt, BuildingTemperatureInfo* building_temperature_info);
};

//----- define Variable -----
extern std::unique_ptr<ConduitTemperatureManager> gConduitTemperatureManager;

#endif