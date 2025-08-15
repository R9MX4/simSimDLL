#include "ClassBase.h"
#include "ClassConduit.h"
#include "SimHandler.h"
#include <Windows.h>
#include <dbghelp.h>

#pragma comment(lib, "Dbghelp.lib")

#define SAVE_VERSION 14

//----- (Decompile Basic) ----------------------------------------------
GameDataUpdate* PrepareGameDataUpdate(const GameData* gameData)
{
#define ASSIGN_GAMEDATA(_NUM, _INFO, _INFO2) \
	gGameDataUpdate._NUM  = (int)gameData->_INFO2.size(); \
	gGameDataUpdate._INFO =      gameData->_INFO2.data();

	LOGGER_PRINT2("%s, ptr %llu\n", __func__, (uint64_t)gameData);
	gGameDataUpdate.numFramesProcessed = gameData->numFramesProcessed;
	gGameDataUpdate.elementIdx         = gameData->cells->elementIdx  .data();
	gGameDataUpdate.temperature        = gameData->cells->temperature .data();
	gGameDataUpdate.mass               = gameData->cells->mass        .data();
	gGameDataUpdate.properties         = gameData->cells->properties  .data();
	gGameDataUpdate.insulation         = gameData->cells->insulation  .data();
	gGameDataUpdate.strengthInfo       = gameData->cells->strengthInfo.data();
	gGameDataUpdate.radiation          = gameData->cells->radiation   .data();
	gGameDataUpdate.diseaseIdx         = gameData->cells->diseaseIdx  .data();
	gGameDataUpdate.diseaseCount       = gameData->cells->diseaseCount.data();

	ASSIGN_GAMEDATA(numSolidInfo						, solidInfo							, solidInfo						);
	ASSIGN_GAMEDATA(numLiquidChangeInfo					, liquidChangeInfo					, liquidChangeInfo				);
	ASSIGN_GAMEDATA(numSolidSubstanceChangeInfo			, solidSubstanceChangeInfo			, solidSubstanceChangeInfo		);
	ASSIGN_GAMEDATA(numSubstanceChangeInfo				, substanceChangeInfo				, substanceChangeInfo			);
	ASSIGN_GAMEDATA(numCallbackInfo						, callbackInfo						, callbackInfo					);
	ASSIGN_GAMEDATA(numSpawnFallingLiquidInfo			, spawnFallingLiquidInfo			, spawnFallingLiquidInfo		);
	ASSIGN_GAMEDATA(numDigInfo							, digInfo							, digInfo						);
	ASSIGN_GAMEDATA(numSpawnOreInfo						, spawnOreInfo						, spawnOreInfo					);
	ASSIGN_GAMEDATA(numSpawnFXInfo						, spawnFXInfo						, spawnFXInfo					);
	ASSIGN_GAMEDATA(numUnstableCellInfo					, unstableCellInfo					, unstableCellInfo				);
	ASSIGN_GAMEDATA(numWorldDamageInfo					, worldDamageInfo					, worldDamageInfo				);
	ASSIGN_GAMEDATA(numBuildingTemperatureInfo			, buildingTemperatureInfo			, buildingTemperatureInfo		);
	ASSIGN_GAMEDATA(numMassConsumedCallbacks			, massConsumedCallbacks				, massConsumedCallbacks			);
	ASSIGN_GAMEDATA(numMassEmittedCallbacks				, massEmittedCallbacks				, massEmittedCallbacks			);
	ASSIGN_GAMEDATA(numDiseaseConsumedCallbacks			, diseaseConsumedCallbacks			, diseaseConsumedCallbacks		);
	ASSIGN_GAMEDATA(numComponentStateChangedMessages	, componentStateChangedMessages		, componentStateChangedMessages	);
	ASSIGN_GAMEDATA(numRemovedMassEntries				, removedMassEntries				, consumedMassInfo				);
	ASSIGN_GAMEDATA(numEmittedMassEntries				, emittedMassEntries				, emittedMassInfo				);
	ASSIGN_GAMEDATA(numElementChunkInfos				, elementChunkInfos					, elementChunkInfo				);
	ASSIGN_GAMEDATA(numElementChunkMeltedInfos			, elementChunkMeltedInfos			, elementChunkMeltedInfo		);
	ASSIGN_GAMEDATA(numBuildingOverheatInfos			, buildingOverheatInfos				, buildingOverheatInfo			);
	ASSIGN_GAMEDATA(numBuildingNoLongerOverheatedInfos	, buildingNoLongerOverheatedInfos	, buildingNoLongerOverheatedInfo);
	ASSIGN_GAMEDATA(numBuildingMeltedInfos				, buildingMeltedInfos				, buildingMeltedInfo			);
	ASSIGN_GAMEDATA(numCellMeltedInfos					, cellMeltedInfos					, cellMeltedInfo				);
	ASSIGN_GAMEDATA(numDiseaseEmittedInfos				, diseaseEmittedInfos				, diseaseEmittedInfo			);
	ASSIGN_GAMEDATA(numDiseaseConsumedInfos				, diseaseConsumedInfos				, diseaseConsumedInfo			);
	ASSIGN_GAMEDATA(numRadiationConsumedCallbacks		, radiationConsumedCallbacks		, radiationConsumedCallbacks	);

	for (LiquidChangeInfo info : gameData->liquidChangeInfo)
		if (info.cellIdx < 0 || info.cellIdx >= gSimData->numGameCells)
			ASSERT_TEXT("LiquidChangeInfo cell index over range");

	gGameDataUpdate.accumulatedFlow                  = gameData->accumulatedFlow.get();
	gGameDataUpdate.propertyTextureFlow              = gameData->propertyTextureFlow.get();
	gGameDataUpdate.propertyTextureLiquid            = gameData->propertyTextureLiquid.get();
	gGameDataUpdate.propertyTextureLiquidData        = gameData->propertyTextureLiquidData.get();
	gGameDataUpdate.propertyTextureExposedToSunlight = gameData->propertyTextureExposedToSunlight.get();

	return &gGameDataUpdate;
}

uint16_t GetElementIndex(uint32_t hash)
{
	auto it = gElementIndices.find(hash);
	if (it == gElementIndices.end()) {
		if (hash) LOGGER_PRINT("%s: Not find element: %d\n", __func__, hash);
		return 0xFFFF;
	}
	else return it->second;
}

void DestroyElementsTable()
{
	gElementNames				.clear();
	gElements					.clear();
	gElementTemperatureData		.clear();
	gElementPostProcessData		.clear();
	gElementLiquidData			.clear();
	gElementPhysicsData			.clear();
	gElementPressureData		.clear();
	gElementPropertyTextureData	.clear();
	gElementLightAbsorptionData	.clear();
	gElementRadiationData		.clear();
	gElementStateData			.clear();
	gElementIndices				.clear();
}

void CleanUp()
{
	if (gSim.get()) {
		if (gSim->mThread.joinable()) {
			gFrameSync.GameSync_CleanUp();
			gSim->mThread.join();
		}
		gSim = nullptr;
	}
	gGameMessageHandler = NULL;
	gSimData.reset();
	gFrameSync.clear();
	DestroyElementsTable();
	gDisease.reset();
}

void CellRead(CellSOA* cells, int cell_idx, BinaryBufferReader* reader, int saved_version)
{
	int hash;
	*reader >> hash;
	if      (hash == 41996153 ) hash = -1736594426; // Cuprite
	else if (hash == 596754987) hash =  1282846257; // MaficRock
	cells->elementIdx[cell_idx] = GetElementIndex(hash);
	//LOGGER_PRINT("Cell Index %d, element: %u->%u\n", cell_idx, hash, cells->elementIdx[cell_idx]);

	*reader >> cells->temperature[cell_idx] >> cells->mass[cell_idx];
	if (std::isinf(cells->temperature[cell_idx]) || std::isnan(cells->temperature[cell_idx])) cells->temperature[cell_idx] = 293;
	if (std::isinf(cells->mass[cell_idx])        || std::isnan(cells->mass       [cell_idx])) cells->mass[cell_idx] = 100;
	if (saved_version >= 14) {
		*reader >> cells->radiation[cell_idx];
		if (std::isinf(cells->radiation[cell_idx]) || std::isnan(cells->radiation[cell_idx])) cells->radiation[cell_idx] = 0;
	}
	if (cells->elementIdx[cell_idx] == 0xFFFF) {
		cells->elementIdx [cell_idx] = GetElementIndex(0x2D39BF75u);
		cells->mass       [cell_idx] = 0.0;
		cells->temperature[cell_idx] = 0.0;
		cells->radiation  [cell_idx] = 0.0;
	}
}

void DoLoadTimeStateTransition(SimData* simData, const int cell)
{
	ElementTemperatureData* data = &gElementTemperatureData[simData->updatedCells->elementIdx[cell]];
	if (simData->updatedCells->temperature[cell] < (data->lowTemp - 3.0) && data->lowTempTransitionIdx != 0xFFFF) {
		simData->updatedCells->temperature[cell] += 1.5;
		simData->updatedCells->elementIdx [cell] = data->lowTempTransitionIdx;
	}
	else if (simData->updatedCells->temperature[cell] > (data->highTemp + 3.0) && data->highTempTransitionIdx != 0xFFFF)
	{
		simData->updatedCells->temperature[cell] -= 1.5;
		simData->updatedCells->elementIdx [cell] = data->highTempTransitionIdx;
	}
	ASSERT_TEMP(simData->updatedCells->mass[cell], simData->updatedCells->temperature[cell]);
}

//----- (SIM Message) --------------------------------------------------
GameDataUpdate* PrepareGameData(BinaryBufferReader* reader)
{
	if (!gSim.get()) return 0;

	gFrameSync.GameSync_PrepareGameData(reader);
	return PrepareGameDataUpdate(gFrameSync.mGameData.get());
}

int CreateElementsTable(BinaryBufferReader* reader)
{
	DestroyElementsTable();

	int numElements;
	*reader >> numElements;
	gElementNames				.resize(numElements);
	gElements					.resize(numElements);
	gElementTemperatureData		.resize(numElements);
	gElementPostProcessData		.resize(numElements);
	gElementLiquidData			.resize(numElements);
	gElementPhysicsData			.resize(numElements);
	gElementPressureData		.resize(numElements);
	gElementPropertyTextureData	.resize(numElements);
	gElementLightAbsorptionData	.resize(numElements);
	gElementRadiationData		.resize(numElements);
	gElementStateData			.resize(numElements);

	for (int i = 0; i < numElements; i++) {
		reader->ReadBytes(sizeof(Element), &gElements[i]);

		auto it = gElementIndices.find(gElements[i].id);
		if (it == gElementIndices.end())
			gElementIndices.insert({ gElements[i].id, gElements[i].elementsTableIdx });
		else
			ASSERT_TEXT("Multiple Element ID");
	}

	for (int i = 0; i < numElements; i++) {
		reader->ReadString(&gElementNames[i]);
		LOGGER_PRINT("Element No.%3d: %s\t, id:%11d, def Temp:%.2f\n", gElements[i].elementsTableIdx, gElementNames[i].c_str(), gElements[i].id, gElements[i].defaultValues.temperature);
	}

	for (int i = 0; i < numElements; i++) {
		gElementTemperatureData[i].state                               = gElements[i].state;
		gElementTemperatureData[i].lowTempTransitionIdx                = gElements[i].lowTempTransitionIdx;
		gElementTemperatureData[i].highTempTransitionIdx               = gElements[i].highTempTransitionIdx;
		gElementTemperatureData[i].lowTempTransitionOreIdx             = gElements[i].lowTempTransitionOreID  == 0x2D39BF75 ? -1 : GetElementIndex(gElements[i].lowTempTransitionOreID);
		gElementTemperatureData[i].highTempTransitionOreIdx            = gElements[i].highTempTransitionOreID == 0x2D39BF75 ? -1 : GetElementIndex(gElements[i].highTempTransitionOreID);
		gElementTemperatureData[i].specificHeatCapacity                = gElements[i].specificHeatCapacity;
		gElementTemperatureData[i].thermalConductivity                 = gElements[i].thermalConductivity;
		gElementTemperatureData[i].defaultMass                         = gElements[i].defaultValues.mass;
		gElementTemperatureData[i].massAreaScale                       = (gElements[i].state & 3) == 1 ? 1 : 0.001f;
		gElementTemperatureData[i].gasSurfaceAreaMultiplier            = gElements[i].gasSurfaceAreaMultiplier;
		gElementTemperatureData[i].liquidSurfaceAreaMultiplier         = gElements[i].liquidSurfaceAreaMultiplier;
		gElementTemperatureData[i].solidSurfaceAreaMultiplier          = gElements[i].solidSurfaceAreaMultiplier;
		gElementTemperatureData[i].lowTemp                             = gElements[i].lowTemp;
		gElementTemperatureData[i].highTemp                            = gElements[i].highTemp;
		gElementTemperatureData[i].lowTempTransitionOreMassConversion  = gElements[i].lowTempTransitionOreMassConversion;
		gElementTemperatureData[i].highTempTransitionOreMassConversion = gElements[i].highTempTransitionOreMassConversion;

		gElementPostProcessData[i].sublimateIndex       = gElements[i].sublimateIndex;
		gElementPostProcessData[i].convertIndex         = gElements[i].convertIndex;
		gElementPostProcessData[i].state                = gElements[i].state;
		gElementPostProcessData[i].sublimateRate        = gElements[i].sublimateRate;
		gElementPostProcessData[i].sublimateEfficiency  = gElements[i].sublimateEfficiency;
		gElementPostProcessData[i].sublimateProbability = gElements[i].sublimateProbability;
		gElementPostProcessData[i].offGasPercentage     = gElements[i].offGasPercentage;
		gElementPostProcessData[i].molarMass            = gElements[i].molarMass;
		gElementPostProcessData[i].strength             = gElements[i].strength;
		gElementPostProcessData[i].maxMass              = gElements[i].maxMass;
		gElementPostProcessData[i].minHorizontalFlow    = gElements[i].minHorizontalFlow;
		gElementPostProcessData[i].sublimateFX          = (CellProperties)gElements[i].sublimateFX;

		gElementLiquidData[i].state             = gElements[i].state;
		gElementLiquidData[i].flow              = gElements[i].flow;
		gElementLiquidData[i].viscosity         = gElements[i].viscosity;
		gElementLiquidData[i].minHorizontalFlow = gElements[i].minHorizontalFlow;
		gElementLiquidData[i].minVerticalFlow   = gElements[i].minVerticalFlow;
		gElementLiquidData[i].maxMass           = gElements[i].maxMass;

		gElementPhysicsData[i].temperature = gElements[i].defaultValues.temperature;
		gElementPhysicsData[i].mass        = gElements[i].defaultValues.mass;
		gElementPhysicsData[i].pressure    = gElements[i].defaultValues.pressure;

		gElementPressureData[i].state       = gElements[i].state;
		gElementPressureData[i].flow        = gElements[i].flow;
		gElementPressureData[i].molarVolume = gElements[i].molarMass > 0 ? 22.4f / gElements[i].molarMass : 1;

		gElementPropertyTextureData[i].colour = gElements[i].colour;
		gElementPropertyTextureData[i].state  = gElements[i].state;

		gElementLightAbsorptionData[i].factor    = gElements[i].lightAbsorptionFactor;
		gElementLightAbsorptionData[i].massScale = (gElements[i].state & 3) == 3 ? INFINITY : 1.0f / gElements[i].defaultValues.mass;
		
		gElementRadiationData[i].factor        = gElements[i].radiationAbsorptionFactor;
		gElementRadiationData[i].rads_per_1000 = gElements[i].radiationPer1000Mass;

		gElementStateData[i].state = gElements[i].state;
	}

	LOGGER_PRINT("%s elem count: %d\n", __func__, numElements);
	return numElements;
}

int CreateElementInteractions(BinaryBufferReader* reader)
{
	gGasObliterations   .clear();
	gLiquidConversions  .clear();
	gLiquidObliterations.clear();

	int count;
	for (*reader >> count; count > 0; count--) {
		Hashes::ElementInteractionHashes handler;
		*reader >> handler;
		switch (handler)
		{
		case Hashes::GasObliteration: {
			GasObliteration info;
			reader->ReadBytes(sizeof(GasObliteration), &info);
			gGasObliterations.push_back(info);
			break;
		}
		case Hashes::LiquidConversion: {
			LiquidConversion info;
			reader->ReadBytes(sizeof(LiquidConversion), &info);
			gLiquidConversions.push_back(info);
			break;
		}
		case Hashes::LiquidObliteration: {
			LiquidObliteration info;
			reader->ReadBytes(sizeof(LiquidObliteration), &info);
			gLiquidObliterations.push_back(info);
			break;
		}
		default:
			break;
		}
	}

	return 0;
}

int CreateElementInteractionsLocked(BinaryBufferReader* reader)
{
	mtx_lock(&gFrameSync.mSimMutex);
	int ElementInteractions = CreateElementInteractions(reader);
	mtx_unlock(&gFrameSync.mSimMutex);
	return ElementInteractions;
}

int CreateDiseaseTable(BinaryBufferReader* reader)
{
	gDisease = std::make_unique<Disease>(reader);
	return 0;
}

GameDataUpdate* Start(BinaryBufferReader* reader)
{
	if (gSimData->initSettleThermalBoundaries) {
		gSimData->cells->CopyFrom(gSimData->updatedCells.get());
		gSimData->SettleThermalBoundaries(gSimData->cells.get(), gSimData->updatedCells.get());
	}
	gSimData->simEvents->substanceChangeInfo           .clear();
	gSimData->simEvents->spawnLiquidInfo               .clear();
	gSimData->simEvents->spawnOreInfo                  .clear();
	gSimData->simEvents->unstableCellInfo              .clear();
	gSimData->simEvents->elementChunkMeltedInfo        .clear();
	gSimData->simEvents->buildingMeltedInfo            .clear();
	gSimData->simEvents->buildingOverheatInfo          .clear();
	gSimData->simEvents->buildingNoLongerOverheatedInfo.clear();
	gSimData->simEvents->cellMeltedInfo                .clear();
	gSimData->simEvents->callbackInfo                  .clear();
	gSimData->simEvents->worldDamageInfo               .clear();
	gSimData->simEvents->massConsumedCallbacks         .clear();
	gSimData->simEvents->radiationConsumedCallbacks    .clear();
	gSimData->simEvents->massEmittedCallbacks          .clear();
	gSimData->simEvents->diseaseConsumedCallbacks      .clear();
	gSimData->simEvents->spawnFXInfo                   .clear();
	gSimData->simEvents->componentStateChangedMessages .clear();
	gSimData->simEvents->digInfo                       .clear();

	gFrameSync.initGameData(gSimData->width - 2, gSimData->height - 2, gSimData->width, gSimData->height, gSimData->updatedCells.get());
	gSimData->InitializeBoundary();
	gSimData->cells->CopyFrom(gSimData->updatedCells.get());
	gSim->Start();
	return PrepareGameDataUpdate(gFrameSync.mGameData.get());
}

SimData* Load(BinaryBufferReader* reader)
{
	std::string head = "SIMSAVE";
	char tmp[9] = {0};
	reader->ReadBytes(8, tmp);
	if (head != tmp) {
		LOGGER_PRINT("%s Not SIMSAVE, Quit!!!\n", __func__);
		return 0;
	}
	int saved_version, base_width, base_height, world_x = 0, world_y = 0;
	uint8_t saved_options = 0;
	*reader >> saved_version;
	if (saved_version >= 15) return 0;
	*reader >> base_width >> base_height;
	if (saved_version >= 14) *reader >> world_x >> world_y;
	if (saved_version >= 13) *reader >> saved_options;
	gSimData->ApplySaveSettings(saved_version);
	LOGGER_PRINT("%s: Save Version %d, Option %d, base WxH %dx%d, world WxH %dx%d\n", __func__, saved_version, saved_options, base_width, base_height, world_x, world_y);

	for (int locY = 0; locY < base_height; locY++) {
		for (int locX = 0; locX < base_width; locX++) {
			int cell = locX + world_x + (world_y + locY) * gSimData->width;
			CellRead(gSimData->updatedCells.get(), cell, reader, saved_version);
			if (gSimData->updatedCells->temperature[cell] <= 0)
				gSimData->updatedCells->temperature[cell] = 293;
			if (gSimData->updatedCells->radiation  [cell] <= 0)
				gSimData->updatedCells->radiation  [cell] = 0;
			if (gSimData->updatedCells->elementIdx [cell] == gSimData->voidElementIdx ||
				gSimData->updatedCells->elementIdx [cell] == gSimData->vacuumElementIdx)
			{
				gSimData->updatedCells->temperature[cell] = 0;
				gSimData->updatedCells->mass       [cell] = 0;
				gSimData->updatedCells->radiation  [cell] = 0;
			}
			DoLoadTimeStateTransition(gSimData.get(), cell);
		}
	}

	if (saved_version >= 9) {
		for (int locY = 0; locY < base_height; locY++) {
			for (int locX = 0; locX < base_width; locX++) {
				int cell = locX + world_x + (world_y + locY) * gSimData->width;
				uint32_t disease_hash;
				*reader >> disease_hash >> gSimData->updatedCells->diseaseCount[cell];
				if (disease_hash) {
					gSimData->updatedCells->diseaseIdx[cell] = gDisease->GetDiseaseIndex(disease_hash);
				}
				else {
					gSimData->updatedCells->diseaseIdx  [cell] = 0xFF;
					gSimData->updatedCells->diseaseCount[cell] = 0;
				}
			}
		}
	}

	if (saved_version >= 8) {
		for (int locY = 0; locY < base_height; locY++) {
			int cell = world_x + (world_y + locY) * gSimData->width;
			reader->ReadBytes(sizeof(float) * base_width, &gSimData->biomeTemperature[cell]);
		}
	}
	//LOGGER_PRINT("%s: Dump biomeTemperature\n", __func__);
	//for (int locY = 0; locY < base_height; locY++) {
	//	for (int locX = 0; locX < base_width; locX++) {
	//		int cell = locX + world_x + (world_y + locY) * gSimData->width;
	//		LOGGER_PRINT("%4.1f ", gSimData->biomeTemperature[cell]);
	//	}
	//	LOGGER_PRINT("\n");
	//}
	if (saved_version < 11)
		gSimData->initSettleThermalBoundaries = true;

	return gSimData.get();
}

int AllocateCells(BinaryBufferReader* reader)
{
	int  width, height;
	bool radiation_enabled, headless;
	*reader >> width >> height >> radiation_enabled >> headless;
	gSimData = std::make_unique<SimData>(width + 2, height + 2, radiation_enabled, headless);

	LOGGER_PRINT("%s width: %d, height: %d, Radiaton:%s, Headless:%s\n", __func__, width, height, radiation_enabled ? "Yes" : "No", headless ? "True" : "False");
	return 0;
}

int DefineWorldOffsets(BinaryBufferReader* reader)
{
	int numWorlds;
	*reader >> numWorlds;
	gSimData->worlds.resize(numWorlds);
	for (int i = 0; i < numWorlds; i++) {
		*reader >> gSimData->worlds[i].offsetX;
		*reader >> gSimData->worlds[i].offsetY;
		*reader >> gSimData->worlds[i].width;
		*reader >> gSimData->worlds[i].height;
	}
	return 0;
}

int ClearUnoccupiedCells()
{
	int cellCnt = gSimData->width * gSimData->height;

	gSimData->cells->elementIdx         .assign(cellCnt, gSimData->vacuumElementIdx);
	gSimData->cells->mass               .assign(cellCnt,  0);
	gSimData->cells->temperature        .assign(cellCnt,  0);
	gSimData->cells->diseaseCount       .assign(cellCnt,  0);
	gSimData->cells->diseaseIdx         .assign(cellCnt, -1);
	gSimData->cells->radiation          .assign(cellCnt,  0);
	gSimData->updatedCells->elementIdx  .assign(cellCnt, gSimData->vacuumElementIdx);
	gSimData->updatedCells->mass        .assign(cellCnt,  0);
	gSimData->updatedCells->temperature .assign(cellCnt,  0);
	gSimData->updatedCells->diseaseCount.assign(cellCnt,  0);
	gSimData->updatedCells->diseaseIdx  .assign(cellCnt, -1);
	gSimData->updatedCells->radiation   .assign(cellCnt,  0);
	memset(gSimData->biomeTemperature.get(), 0, cellCnt * sizeof(float));
	return 0;
}

int SetSavedOptions(BinaryBufferReader* reader)
{
	uint8_t clear, set;
	*reader >> clear >> set;
	if (gSimData.get()) {
		gSimData->savedOptions &= ~clear;
		gSimData->savedOptions |=  set;
	}
	return 0;
}

//----- (SIM Handler) --------------------------------------------------
void SIM_Initialize(void* (*callback)(int, const void*))
{
	if (gLogger == NULL) {
		INIT_TIMER();
		LOGGER_INIT();
		LOGGER_PRINT("%s Init Logger\n\n", __func__);

		LOGGER_INIT2();
		LOGGER_PRINT2("%s Init Logger\n\n", __func__);
	}

	//Timer::Initialize();
	HANDLE CurrentProcess = GetCurrentProcess();
	SymInitialize(CurrentProcess, NULL, false);
	//SetUnhandledExceptionFilter((LPTOP_LEVEL_EXCEPTION_FILTER)SimDLLUnhandledExceptionFilter);
	CleanUp();
	gGameMessageHandler = callback;
	gFrameSync.clear();

	gSim = nullptr;
	gSim = std::make_unique<Sim>(&gFrameSync);
}

void SIM_Shutdown() {
	LOGGER_PRINT("%s\n", __func__);

	CleanUp();
	SymCleanup(GetCurrentProcess());
}

int64_t SIM_HandleMessage(Hashes::SimMessageHashes sim_msg_id, int msg_length, char* msg)
{
	LOGGER_PRINT("%s Command:%d,%d,%lld\n", __func__, sim_msg_id, msg_length, (int64_t)msg);

	if (!gSim)
		return 0;
	if (gSim->simFrameManager.HandleMessage(sim_msg_id, msg_length, msg))
		return 0;

	BinaryBufferReader reader(msg_length, msg);
	switch (sim_msg_id) {
	case Hashes::PrepareGameData:                        return (int64_t)PrepareGameData(&reader);
	case Hashes::Elements_CreateTable:                   return CreateElementsTable(&reader);
	case Hashes::Elements_CreateInteractions:            return CreateElementInteractionsLocked(&reader);
	case Hashes::Disease_CreateTable:                    return CreateDiseaseTable(&reader);
	case Hashes::SetWorldZones:                          return SimData::SetWorldZones(&reader);
	case Hashes::SimData_InitializeFromCells:            return SimData::InitializeFromCells(&reader);
	case Hashes::SimData_ResizeAndInitializeVacuumCells: return SimData::ResizeAndInitializeVacuumCells(&reader);
	case Hashes::SimData_FreeCells:                      return SimData::FreeGridCells(&reader);
	case Hashes::Load:                 return (int64_t)Load(&reader);
	case Hashes::Start:                return (int64_t)Start(&reader);
	case Hashes::AllocateCells:        return AllocateCells(&reader);
	case Hashes::DefineWorldOffsets:   return DefineWorldOffsets(&reader);
	case Hashes::ClearUnoccupiedCells: return ClearUnoccupiedCells();
	case Hashes::ToggleProfiler:       return 0; //std::streambuf::sync
	case Hashes::SetSavedOptions:      return SetSavedOptions(&reader);
	default:
		ASSERT_TEXT("Unknow SIM_Handle");
		return 0;
	}
}

char* SIM_BeginSave(uint32_t* out_size, int x, int y)
{
	LOGGER_PRINT("%s. x = %d, y = %d\n", __func__, x, y);

	if (SaveBuffer) {
		delete(SaveBuffer);
		SaveBuffer = 0;
	}

	int cellCnt  = gSimData->width * gSimData->height;
	*out_size = 28 * cellCnt + 29;
	SaveBuffer = new Buffer(*out_size);
	BinaryBufferWriter writer(SaveBuffer);
	LOGGER_PRINT("%s. SaveBuffer = %lld, size = %d\n", __func__, (uint64_t)SaveBuffer, *out_size);

	writer.WriteBytes(8, (void*)"SIMSAVE");
	writer << SAVE_VERSION;
	writer << gSimData->width << gSimData->height << x << y << gSimData->savedOptions;

	for (int cell = 0; cell < cellCnt; cell++) {
		writer << gElements[gSimData->updatedCells->elementIdx [cell]].id
			   <<           gSimData->updatedCells->temperature[cell]
			   <<           gSimData->updatedCells->mass       [cell]
			   <<           gSimData->updatedCells->radiation  [cell];
	}
	for (int cell = 0; cell < cellCnt; cell++) {
		uint32_t disease_hash = 0;
		if (gSimData->updatedCells->diseaseIdx[cell] != 0xFF)
			disease_hash = gDisease->diseases[gSimData->updatedCells->diseaseIdx[cell]].hashID;
		writer << disease_hash << gSimData->updatedCells->diseaseCount[cell];
	}

	writer.WriteBytes(sizeof(float) * cellCnt, gSimData->biomeTemperature.get());

	LOGGER_PRINT("%s done. SaveBuffer data:%lld\n", __func__, (uint64_t)SaveBuffer->mData);
	return SaveBuffer->mData;
}

void SIM_EndSave()
{
	LOGGER_PRINT("%s\n", __func__);

	if (SaveBuffer) {
		delete(SaveBuffer);
		SaveBuffer = 0;
	}
}

void SIM_DebugCrash() {
	LOGGER_PRINT("%s\n", __func__);
}

//----- (Sim Class) ----------------------------------------------------
Sim::Sim(FrameSync* fs) : SimBase()
{
	this->mName          = "SimThread";
	this->mThread        = std::thread();
	this->mExitRequested = false;
	this->frameSync      = fs;
}

void Sim::Main()
{
	while (!this->mExitRequested) {
		LOGGER_PRINT2("--Sim %s Update, frameSync %lld\n", __func__, (uint64_t)this->frameSync);
		mtx_lock(&this->frameSync->mSimMutex);
		LOGGER_PRINT2("--Sim %s Update, gConduitTemperatureManager %lld\n", __func__, (uint64_t)gConduitTemperatureManager.get());
		if (gConduitTemperatureManager.get())
			gConduitTemperatureManager->ReleaseQueuedHandles();

		int64_t num_frames_processed = this->simFrameManager.BeginFrameProcessing();
		LOGGER_PRINT2("--Sim %s frame:%lld\n", __func__, num_frames_processed);
		memset(gSimData->flow.get(), 0, sizeof(Vector4 < float>) * gSimData->width * gSimData->height);
		if (num_frames_processed <= 0) {
			gSimData->UpdateComponentsDataListOnly();
		}
		else {
			for (int i = 0; i < num_frames_processed; i++) {
				LOGGER_PRINT2("--Sim %s process:%d/%lld\n", __func__, i, num_frames_processed);
				float Frame = this->simFrameManager.ProcessNextFrame(gSimData.get());
				LOGGER_PRINT2("--Sim %s NextFrame:%f\n", __func__, Frame);
				if (Frame <= 0) {
					this->CopyUpdatedCellsToCells(gSimData.get());
					gSimData->UpdateComponentsDataListOnly();
					continue;
				}

				Frame += this->elapsedSeconds;
				int frameCnt = (int)MAX_F(1, 5 * Frame);
				this->elapsedSeconds = min(0, Frame - frameCnt * 0.2f);
				while (frameCnt) {
					LOGGER_PRINT2("--Sim %s UpdateData, frame: %d\n", __func__, frameCnt);
					this->UpdateData(gSimData.get());
					frameCnt--;
				}
			}
		}
		LOGGER_PRINT2("--Sim %s Start EndFrameProcessing\n", __func__);
		this->simFrameManager.EndFrameProcessing();
		this->CopySimDataToGame(gSimData.get(), frameSync->mSimData.get(), frameSync->mGameData.get(), (int)num_frames_processed);
		mtx_unlock(&this->frameSync->mSimMutex);
		this->frameSync->SimSync();
		LOGGER_PRINT2("--Sim %s. mSimReadyToSwap:%s, mExitRequested:%s\n\n", __func__, this->frameSync->mSimReadyToSwap ? "T" : "F", this->mExitRequested ? "T" : "F");
	}
	LOGGER_PRINT2("--Sim %s done\n", __func__);
}

void Sim::Start()
{
	if (this->mThread.joinable()) {
		if (!this->mExitRequested)
			this->mExitRequested = true;
		this->mThread.join();
	}
	this->mExitRequested = false;

	this->mThread = std::thread(&Sim::Main, this);
}