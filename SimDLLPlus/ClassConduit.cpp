#include "ClassBase.h"
#include "ClassConduit.h"

std::unique_ptr<ConduitTemperatureManager> gConduitTemperatureManager;

//----- (CompactedVector) --------------------------------------------
template <typename T>
Handle CompactedVector<T>::Free(int handle)
{
	int     handleVal = handle & 0xFFFFFF;
	uint8_t handleVer = (handle >> 24) + 1;
	if (handleVal >= 0x100000)
		return handle;

	this->handles.versions[handleVal] = handleVer;
	this->handles.freeHandles.push_back(handleVal | (handleVer << 24));
	int handleItm = this->handles.items[handleVal];
	this->handles.items[handleVal] = 0;
	// Replace the target handle with the last handle
	int posEnd = (int)(this->data.size() - 1);
	if (handleItm < posEnd) {
		this->data[handleItm] = this->data[posEnd];
		int handleVal = this->dataHandleIndices[posEnd] & 0xFFFFFF;
		this->handles.items[handleVal] = handleItm;
		this->dataHandleIndices[handleItm] = this->dataHandleIndices[posEnd];
	}
	this->data.pop_back();
	this->dataHandleIndices.pop_back();
	return 0xFFFFFF;
}

template <typename T>
Handle CompactedVector<T>::AddData(T* new_data)
{
	Handle result = -1;
	if (this->handles.freeHandles.size()) {
		result = this->handles.freeHandles.back();
		this->handles.freeHandles.pop_back();
		this->handles.items[result & 0xFFFFFF] = (Handle)this->data.size();
	}
	else {
		this->handles.versions.push_back(0);
		result = this->handles.items.size() & 0xFFFFFF;
		this->handles.items.push_back((Handle)this->data.size());
	}
	this->dataHandleIndices.push_back(result);
	this->data.push_back(*new_data);

	return result;
}

template <typename T>
T* CompactedVector<T>::GetData(int handle)
{
	if (!HANDLE_AVAILABLE(handle, this->handles.versions)) {
		ASSERT_TEXT("Illegal Handle or Version");
		LOGGER_PRINT("Handle %d, Version %d\n", handle, this->handles.versions[handle & 0xFFFFFF]);
	}
	return &this->data[this->handles.items[handle & 0xFFFFFF]];
}

//----- (ConduitTemperatureManager) ----------------------------------
Handle ConduitTemperatureManager::Add(float contents_temperature, float contents_mass, uint32_t contents_element_hash,
	int conduit_structure_temperature_handle, float conduit_heat_capacity, float conduit_thermal_conductivity, bool conduit_insulated)
{
	mtx_lock(&this->dataMutex);

	Element& elem = gElements[GetElementIndex(contents_element_hash)];
	if ((contents_temperature <= 0.0 && contents_mass > 0.0) || contents_temperature > SIM_MAX_TEMPERATURE) {
		ASSERT_TEXT("Assert Temperature failed");
		contents_temperature = elem.defaultValues.temperature;
	}
	ConduitTemperatureManager::Data data = {
		.temperature                      = contents_temperature,
		.thermalConductivity              = elem.thermalConductivity,
		.heatCapacity                     = contents_mass * elem.specificHeatCapacity,
		.conduitBuildingTemperatureHandle = conduit_structure_temperature_handle,
		.conduitHeatCapacity              = conduit_heat_capacity,
		.conduitThermalConductivity       = conduit_thermal_conductivity,
		.conduitInsulated                 = conduit_insulated,
		.lowStateTransitionTemperature    = elem.lowTempTransitionIdx  != 0xFFFF ? elem.lowTemp  : 0,
		.highStateTransitionTemperature   = elem.highTempTransitionIdx != 0xFFFF ? elem.highTemp : FLT_MAX
	};
	Handle result = this->data.AddData(&data);

	mtx_unlock(&this->dataMutex);

	LOGGER_PRINT("ConduitTemperatureManager::%s done.\n", __func__);
	return result;
};

void ConduitTemperatureManager::ReleaseQueuedHandles()
{
	LOGGER_PRINT2("ConduitTemperatureManager::%s.\n", __func__);
	mtx_lock(&this->dataMutex);

	this->curReleaseListIdx ^= 1;
	for (Handle handle: this->frameDelayedReleasedHandles[this->curReleaseListIdx]){
		this->data.Free(handle);
	}
	this->frameDelayedReleasedHandles[this->curReleaseListIdx].clear();

	mtx_unlock(&this->dataMutex);
	LOGGER_PRINT2("ConduitTemperatureManager::%s done.\n", __func__);
}

void ConduitTemperatureManager::Set(int handle, float contents_temperature, float contents_mass, uint32_t contents_element_hash)
{
	mtx_lock(&this->dataMutex);

	if (HANDLE_AVAILABLE(handle, this->data.handles.versions)) {
		ConduitTemperatureManager::Data* p_data = gConduitTemperatureManager->data.GetData(handle);
		Element& elem = gElements[GetElementIndex(contents_element_hash)];
		if ((contents_temperature <= 0.0 && contents_mass > 0.0) || contents_temperature > SIM_MAX_TEMPERATURE) {
			ASSERT_TEXT("Assert Temperature failed");
			contents_temperature = elem.defaultValues.temperature;
		}

		p_data->temperature         = contents_temperature;
		p_data->heatCapacity        = contents_mass * elem.specificHeatCapacity;
		p_data->thermalConductivity = elem.thermalConductivity;
		p_data->lowStateTransitionTemperature  = elem.lowTempTransitionIdx  != 0xFFFF ? elem.lowTemp  : 0;
		p_data->highStateTransitionTemperature = elem.highTempTransitionIdx != 0xFFFF ? elem.highTemp : FLT_MAX;
	}
	else {
		ASSERT_TEXT("Illegal Handle or Version");
		LOGGER_PRINT("Handle %d, Version %d\n", handle, this->data.handles.versions[handle & 0xFFFFFF]);
	}

	mtx_unlock(&this->dataMutex);

	//LOGGER_PRINT("ConduitTemperatureManager::%s done.\n", __func__);
}

ConduitTemperatureUpdateData* ConduitTemperatureManager::Update(float dt, BuildingTemperatureInfo* building_temperature_info)
{
	mtx_lock(&this->dataMutex);

	this->meltedContentHandles.clear();
	this->frozenContentHandles.clear();
	this->temperatures.resize(this->data.handles.items.size());

	if (building_temperature_info) {
		for (int i = 0; i < this->data.data.size(); i++) {
			ConduitTemperatureManager::Data& conduitData = this->data.data[i];

			int handleVal = this->data.dataHandleIndices[i] & 0xFFFFFF;
			int conduitHandle = conduitData.conduitBuildingTemperatureHandle & 0xFFFFFF;
			if (conduitHandle >= 0x100000)
				continue;

			float temperaturePipe = building_temperature_info[conduitHandle].temperature;
			if (conduitData.heatCapacity <= 0.0001 || conduitData.conduitHeatCapacity <= 0.0001 || temperaturePipe <= 0.0) {
				this->temperatures[handleVal] = conduitData.temperature;
				continue;
			}
			ASSERT_TEMP(1, conduitData.temperature);
			ASSERT_TEMP(1, temperaturePipe);

			float TC;
			if (conduitData.conduitInsulated)
				TC = MIN_F(conduitData.thermalConductivity, conduitData.conduitThermalConductivity);
			else
				TC = 0.5f * (conduitData.thermalConductivity + conduitData.conduitThermalConductivity);
			float heatTrans   = (conduitData.temperature - temperaturePipe) * TC * 50 * dt * gSimData->debugProperties.buildingToBuildingTemperatureScale; // 0.001
			float tempMax     = MAX_F(conduitData.temperature, temperaturePipe);
			float tempMin     = MIN_F(conduitData.temperature, temperaturePipe);
			float tempDeltaCT = CLAMP_F(conduitData.temperature - heatTrans / conduitData.heatCapacity,        tempMax, tempMin) - conduitData.temperature;
			float tempDeltaPi = CLAMP_F(temperaturePipe             + heatTrans / conduitData.conduitHeatCapacity, tempMax, tempMin) - temperaturePipe;
			heatTrans         = MIN_F(fabsf(tempDeltaCT) * conduitData.heatCapacity, fabsf(tempDeltaPi) * conduitData.conduitHeatCapacity);
			if (conduitData.temperature < temperaturePipe)
				heatTrans = -heatTrans;

			float tempFinCT = CLAMP_F(conduitData.temperature - heatTrans / conduitData.heatCapacity       , SIM_MAX_TEMPERATURE, 0);
			float tempFinPi = CLAMP_F(temperaturePipe             + heatTrans / conduitData.conduitHeatCapacity, SIM_MAX_TEMPERATURE, 0);
			if ((tempFinCT - tempFinPi) * (conduitData.temperature - temperaturePipe) < 0) {
				tempFinCT = (conduitData.heatCapacity * conduitData.temperature + conduitData.conduitHeatCapacity * temperaturePipe) / (conduitData.heatCapacity + conduitData.conduitHeatCapacity);
			}

			ASSERT_TEMP(1, tempFinCT);
			// ASSERT_TEMP(1, tempFinCD);
			this->temperatures[handleVal] = tempFinCT;
			conduitData.temperature       = tempFinCT;

			if (fabsf(heatTrans) > 1e-6) {
				ModifyBuildingEnergyMessage msg = {
					.handle         =  conduitData.conduitBuildingTemperatureHandle,
					.deltaKJ        = (conduitData.temperature - tempFinCT) * conduitData.heatCapacity,
					.minTemperature = tempMin,
					.maxTemperature = tempMax };
				gSim->simFrameManager.HandleMessage(Hashes::ModifyBuildingEnergy, 16, (char*)&msg);
			}

			if (tempFinCT < conduitData.lowStateTransitionTemperature - 3)
				this->frozenContentHandles.push_back(handleVal);
			if (tempFinCT > conduitData.highStateTransitionTemperature + 3)
				this->meltedContentHandles.push_back(handleVal);
		}
	}
	else {
		for (int i = 0; i < this->data.data.size(); i++) {
			int handleVal = this->data.dataHandleIndices[i] & 0xFFFFFF;
			this->temperatures[handleVal] = this->data.data[i].temperature;
		}
	}

	this->updateData.numEntries       = (int)this->temperatures.size();
	this->updateData.temperatures     =      this->temperatures.data();
	this->updateData.numFrozenHandles = (int)this->frozenContentHandles.size();
	this->updateData.frozenHandles    =      this->frozenContentHandles.data();
	this->updateData.numMeltedHandles = (int)this->meltedContentHandles.size();
	this->updateData.meltedHandles    =      this->meltedContentHandles.data();

	mtx_unlock(&this->dataMutex);

	LOGGER_PRINT("ConduitTemperatureManager::%s done. Entries:%d; Frozen:%d; Melted:%d\n", __func__, this->updateData.numEntries, this->updateData.numFrozenHandles, this->updateData.numMeltedHandles);
	return &this->updateData;
}

//----- (EA_ConduitTemperatureManager Handler) -----------------------
void ConduitTemperatureManager_Initialize()
{
	LOGGER_PRINT("%s\n", __func__);
	gConduitTemperatureManager = std::make_unique<ConduitTemperatureManager>();
	LOGGER_PRINT("%s done. gConduitTemperatureManager: %lld\n", __func__, (uint64_t)gConduitTemperatureManager.get());
}

void ConduitTemperatureManager_Shutdown()
{
	LOGGER_PRINT("%s\n", __func__);
	gConduitTemperatureManager.reset();
	LOGGER_PRINT("%s done. gConduitTemperatureManager: %lld\n", __func__, (uint64_t)gConduitTemperatureManager.get());
}

int ConduitTemperatureManager_Add(float contents_temperature, float contents_mass, int contents_elem_hash,
	int conduit_structure_temperature_handle, float conduit_heat_capacity, float conduit_thermal_conductivity, bool conduit_insulated)
{
	LOGGER_PRINT("%s\n", __func__);
	return gConduitTemperatureManager->Add(contents_temperature, contents_mass, contents_elem_hash, 
		conduit_structure_temperature_handle, conduit_heat_capacity, conduit_thermal_conductivity, conduit_insulated);
}

void ConduitTemperatureManager_Set(int handle, float contents_temperature, float contents_mass, int contents_element_hash)
{
	LOGGER_PRINT("%s\n", __func__);
	gConduitTemperatureManager->Set(handle, contents_temperature, contents_mass, contents_element_hash);
	LOGGER_PRINT("%s done.\n", __func__);
}

void ConduitTemperatureManager_Remove(int handle)
{
	LOGGER_PRINT("%s\n", __func__);
	mtx_lock(&gConduitTemperatureManager->dataMutex);

	ConduitTemperatureManager::Data* Data = gConduitTemperatureManager->data.GetData(handle);
	Data->conduitHeatCapacity              = -1;
	Data->conduitThermalConductivity       = -1;
	Data->conduitBuildingTemperatureHandle = -1;
	gConduitTemperatureManager->frameDelayedReleasedHandles[gConduitTemperatureManager->curReleaseListIdx].push_back(handle);

	mtx_unlock(&gConduitTemperatureManager->dataMutex);
	LOGGER_PRINT("%s done.\n", __func__);
}

void* ConduitTemperatureManager_Update(float dt, BuildingTemperatureInfo* building_temperatures)
{
	LOGGER_PRINT("%s. gConduitTemperatureManager: %lld\n", __func__, (uint64_t)gConduitTemperatureManager.get());
	return (void*)gConduitTemperatureManager->Update(dt, building_temperatures);
}

void ConduitTemperatureManager_Clear()
{
	LOGGER_PRINT("%s\n", __func__);
	mtx_lock(&gConduitTemperatureManager->dataMutex);

	gConduitTemperatureManager->data.dataHandleIndices  .clear();
	gConduitTemperatureManager->data.data               .clear();
	gConduitTemperatureManager->data.handles.items      .clear();
	gConduitTemperatureManager->data.handles.freeHandles.clear();
	gConduitTemperatureManager->data.handles.versions   .clear();

	mtx_unlock(&gConduitTemperatureManager->dataMutex);
	LOGGER_PRINT("%s done.\n", __func__);
}