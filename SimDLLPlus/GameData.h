#pragma once
#ifndef GAMEDATA_H
#define GAMEDATA_H

#include "global.h"
#include "ClassSim.h"
//----- define struct -----
struct SolidInfo
{
	int cellIdx;
	int solid; // Must use int instead of bool, otherwise the compiler may assign 7F0X instead of 000X
};
struct LiquidChangeInfo
{
	int cellIdx;
};
struct SolidSubstanceChangeInfo
{
	int cellIdx;
};
#pragma pack(push, 1)
struct alignas(1) GameDataUpdate
{
	int             numFramesProcessed;
	const uint16_t* elementIdx;
	const float*    temperature;
	const float*    mass;
	const uint8_t*  properties;
	const uint8_t*  insulation;
	const uint8_t*  strengthInfo;
	const float*    radiation;
	const uint8_t*  diseaseIdx;
	const int*      diseaseCount;

	int numSolidInfo;
	const SolidInfo* solidInfo;
	int numLiquidChangeInfo;
	const LiquidChangeInfo* liquidChangeInfo;
	int numSolidSubstanceChangeInfo;
	const SolidSubstanceChangeInfo* solidSubstanceChangeInfo;
	int numSubstanceChangeInfo;
	const SubstanceChangeInfo* substanceChangeInfo;
	int numCallbackInfo;
	const CallbackInfo* callbackInfo;
	int numSpawnFallingLiquidInfo;
	const SpawnFallingLiquidInfo* spawnFallingLiquidInfo;
	int numDigInfo;
	const SpawnOreInfo* digInfo;
	int numSpawnOreInfo;
	const SpawnOreInfo* spawnOreInfo;
	int numSpawnFXInfo;
	const SpawnFXInfo* spawnFXInfo;
	int numUnstableCellInfo;
	const UnstableCellInfo* unstableCellInfo;
	int numWorldDamageInfo;
	const WorldDamageInfo* worldDamageInfo;
	int numBuildingTemperatureInfo;
	const BuildingTemperatureInfo* buildingTemperatureInfo;
	int numMassConsumedCallbacks;
	const MassConsumedCallback* massConsumedCallbacks;
	int numMassEmittedCallbacks;
	const MassEmittedCallback* massEmittedCallbacks;
	int numDiseaseConsumedCallbacks;
	const DiseaseConsumedCallback* diseaseConsumedCallbacks;
	int numComponentStateChangedMessages;
	const ComponentStateChangedMessage* componentStateChangedMessages;
	int numRemovedMassEntries;
	const ConsumedMassInfo* removedMassEntries;
	int numEmittedMassEntries;
	const EmittedMassInfo* emittedMassEntries;
	int numElementChunkInfos;
	const ElementChunkInfo* elementChunkInfos;
	int numElementChunkMeltedInfos;
	const MeltedInfo* elementChunkMeltedInfos;
	int numBuildingOverheatInfos;
	const MeltedInfo* buildingOverheatInfos;
	int numBuildingNoLongerOverheatedInfos;
	const MeltedInfo* buildingNoLongerOverheatedInfos;
	int numBuildingMeltedInfos;
	const MeltedInfo* buildingMeltedInfos;
	int numCellMeltedInfos;
	const CellMeltedInfo* cellMeltedInfos;
	int numDiseaseEmittedInfos;
	const DiseaseEmittedInfo* diseaseEmittedInfos;
	int numDiseaseConsumedInfos;
	const DiseaseConsumedInfo* diseaseConsumedInfos;
	int numRadiationConsumedCallbacks;
	const RadiationConsumedCallback* radiationConsumedCallbacks;

	const float*          accumulatedFlow;
	const Vector2<float>* propertyTextureFlow;
	const void*           propertyTextureLiquid;
	const void*           propertyTextureLiquidData;
	const uint8_t*        propertyTextureExposedToSunlight;
};
#pragma pack(pop)
//----- define class -----
class GameData
{
public:
	int width;
	int height;
	int numFramesProcessed;
	std::vector<SubstanceChangeInfo>		substanceChangeInfo;
	std::vector<CallbackInfo>				callbackInfo;
	std::vector<SpawnFallingLiquidInfo>		spawnFallingLiquidInfo;
	std::vector<SpawnOreInfo>				spawnOreInfo;
	std::vector<SpawnOreInfo>				digInfo;
	std::vector<UnstableCellInfo>			unstableCellInfo;
	std::vector<WorldDamageInfo>			worldDamageInfo;
	std::vector<BuildingTemperatureInfo>	buildingTemperatureInfo;
	std::vector<MassConsumedCallback>		massConsumedCallbacks;
	std::vector<MassEmittedCallback>		massEmittedCallbacks;
	std::vector<DiseaseConsumedCallback>	diseaseConsumedCallbacks;
	std::vector<RadiationConsumedCallback>	radiationConsumedCallbacks;
	std::vector<SpawnFXInfo>				spawnFXInfo;
	std::vector<SolidInfo>					solidInfo;
	std::vector<LiquidChangeInfo>			liquidChangeInfo;
	std::vector<SolidSubstanceChangeInfo>	solidSubstanceChangeInfo;
	std::vector<ConsumedMassInfo>			consumedMassInfo;
	std::vector<EmittedMassInfo>			emittedMassInfo;
	std::vector<ElementChunkInfo>			elementChunkInfo;
	std::vector<DiseaseEmittedInfo>			diseaseEmittedInfo;
	std::vector<DiseaseConsumedInfo>		diseaseConsumedInfo;
	std::vector<MeltedInfo>	elementChunkMeltedInfo;
	std::vector<MeltedInfo>	buildingMeltedInfo;
	std::vector<MeltedInfo>	buildingOverheatInfo;
	std::vector<MeltedInfo>	buildingNoLongerOverheatedInfo;
	std::vector<ComponentStateChangedMessage> componentStateChangedMessages;
	std::vector<CellMeltedInfo>	cellMeltedInfo;
	std::unique_ptr<CellSOA>	cells;
	std::unique_ptr<Vector4<float>[], std::default_delete<Vector4<float>[]>> flow;
	std::unique_ptr<float[]>	accumulatedFlow;
	std::unique_ptr<uint8_t[]>	visibleGrid;
	std::unique_ptr<Vector2<float>[], std::default_delete<Vector2<float>[]>> propertyTextureFlow;
	std::unique_ptr<uint32_t[]>	propertyTextureLiquid;
	std::unique_ptr<uint32_t[]>	propertyTextureLiquidData;
	std::unique_ptr<uint8_t[]>	propertyTextureExposedToSunlight;

	//void swapVisibleGrid(std::unique_ptr<uint8_t[]>* other);
	explicit GameData(int width, int height) {
		this->width              = width;
		this->height             = height;
		this->numFramesProcessed = 0;

		this->substanceChangeInfo        = std::vector<SubstanceChangeInfo>();
		this->callbackInfo               = std::vector<CallbackInfo>();
		this->spawnFallingLiquidInfo     = std::vector<SpawnFallingLiquidInfo>();
		this->spawnOreInfo               = std::vector<SpawnOreInfo>();
		this->digInfo                    = std::vector<SpawnOreInfo>();
		this->unstableCellInfo           = std::vector<UnstableCellInfo>();
		this->worldDamageInfo            = std::vector<WorldDamageInfo>();
		this->buildingTemperatureInfo    = std::vector<BuildingTemperatureInfo>();
		this->massConsumedCallbacks      = std::vector<MassConsumedCallback>();
		this->massEmittedCallbacks       = std::vector<MassEmittedCallback>();
		this->diseaseConsumedCallbacks   = std::vector<DiseaseConsumedCallback>();
		this->radiationConsumedCallbacks = std::vector<RadiationConsumedCallback>();
		this->spawnFXInfo                = std::vector<SpawnFXInfo>();
		this->solidInfo                  = std::vector<SolidInfo>();
		this->liquidChangeInfo           = std::vector<LiquidChangeInfo>();
		this->solidSubstanceChangeInfo   = std::vector<SolidSubstanceChangeInfo>();
		this->consumedMassInfo           = std::vector<ConsumedMassInfo>();
		this->emittedMassInfo            = std::vector<EmittedMassInfo>();
		this->elementChunkInfo           = std::vector<ElementChunkInfo>();
		this->diseaseEmittedInfo         = std::vector<DiseaseEmittedInfo>();
		this->diseaseConsumedInfo        = std::vector<DiseaseConsumedInfo>();
		this->elementChunkMeltedInfo         = std::vector<MeltedInfo>();
		this->buildingMeltedInfo             = std::vector<MeltedInfo>();
		this->buildingOverheatInfo           = std::vector<MeltedInfo>();
		this->buildingNoLongerOverheatedInfo = std::vector<MeltedInfo>();
		this->componentStateChangedMessages  = std::vector<ComponentStateChangedMessage>();
		this->cellMeltedInfo                 = std::vector<CellMeltedInfo>();

		this->cells                           .reset();
		this->flow                            .reset();
		this->accumulatedFlow                 .reset();
		this->visibleGrid                     .reset();
		this->propertyTextureFlow             .reset();
		this->propertyTextureLiquid           .reset();
		this->propertyTextureLiquidData       .reset();
		this->propertyTextureExposedToSunlight.reset();

		int cellCount = height * width;
		this->cells           = std::make_unique<CellSOA>         (cellCount);
		this->flow            = std::make_unique<Vector4<float>[]>(cellCount);
		this->accumulatedFlow = std::make_unique<float[]>         (cellCount);
		this->visibleGrid     = std::make_unique<uint8_t[]>       (cellCount);
		this->propertyTextureFlow              = std::make_unique<Vector2<float>[]>(cellCount);
		this->propertyTextureLiquid            = std::make_unique<uint32_t[]>      (cellCount);
		this->propertyTextureLiquidData        = std::make_unique<uint32_t[]>      (cellCount);
		this->propertyTextureExposedToSunlight = std::make_unique<uint8_t[]>       (cellCount);
	}
};
#endif