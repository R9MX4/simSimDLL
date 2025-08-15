#pragma once
#ifndef CLASS_SIM_H
#define CLASS_SIM_H

#include "global.h"
#include "ClassCell.h"
//----- define struct -----
struct Region
{
	Vector2<int> minimum;
	Vector2<int> maximum;
};
struct ActiveRegion
{
	Region region;
	float currentSunlightIntensity;
	float currentCosmicRadiationIntensity;
};
struct DebugProperties
{
	float buildingTemperatureScale;				// 100
	float buildingToBuildingTemperatureScale;	// 0.001
	float biomeTemperatureLerpRate;				// 0.001
	bool  isDebugEditing;
};
struct UnstableOffsetInfo //DoUnstableCheckWithDiagonals::__l2::UnstableOffsetInfo
{
	int   offset;
	float massPercent;
	float srcMinMassScale;
	bool  checkAbove;
};
struct FloodFillInfo
{
	uint16_t x;
	uint16_t y;
	uint16_t depth;
};
struct InitializeFromCellsMessage
{
	struct Cell // Sim.Cell.Write(BinaryWriter)
	{
		uint16_t elementIdx;
		uint8_t  properties;
		uint8_t  insulation;
		uint8_t  strengthInfo;
		uint8_t  pad[3];
		float temperature;
		float mass;
	};

	struct DiseaseCell // Sim.DiseaseCell.Write(BinaryWriter)
	{
		uint8_t diseaseIdx;
		uint8_t infestationTickCount;
		uint8_t pad[2];
		int     count;
		float   accumulatedError;
	};

	int  width;
	int  height;
	int  randSeed;
	bool radiation_enabled;
	bool headless;
	Cell*        cells;
	float*       biomeTemperature;
	DiseaseCell* disease;
};
struct ResizeAndInitializeVacuumCellsMessage
{
	int grid_width;
	int grid_height;
	int width;
	int height;
	int x_offset;
	int y_offset;
};

//----- info struct -----
struct ConsumedMassInfo
{
	Handle   simHandle;
	uint16_t removedElemIdx;
	uint8_t  diseaseIdx;
	//	uint8_t  pad;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct EmittedMassInfo
{
	uint16_t elemIdx;
	uint8_t  diseaseIdx;
	//	uint8_t  pad;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct ElementChunkInfo
{
	float temperature;
	float deltaKJ;
};
struct BuildingTemperatureInfo
{
	Handle handle;
	float temperature;
};
struct DiseaseEmittedInfo
{
	uint8_t	diseaseIdx;
	int		count;
};
struct DiseaseConsumedInfo
{
	uint8_t diseaseIdx;
	int     count;
};

//----- data struct -----
struct ElementConsumerData
{
	float    consumptionRate;
	uint8_t  maxDepth;
	uint8_t  configuration; // 0:Element, 1:AllLiquid, 2:AllGas
	uint16_t elemIdx;
	uint8_t  offsetIdx;
	Vector2<uint16_t> cell;
};
struct ElementEmitterData
{
	float elapsedTime;
	float emitInterval;
	float emitMass;
	float emitTemperature;
	float maxPressure;
	int   emitDiseaseCount;
	Vector2<uint16_t> cell;
	uint16_t elemIdx;
	uint8_t  maxDepth;
	uint8_t  offsetIdx;
	uint8_t  blockedState;
	uint8_t  diseaseIdx;
	uint32_t blockedCBIdx;
	uint32_t unblockedCBIdx;
};
struct RadiationEmitterData
{
	int cell;
	uint16_t radiusX;
	uint16_t radiusY;
	float emitRads;  // Ridiation dose
	float emitRate;  // Frequence, small -> emit fast
	float emitSpeed; // Fixed 1 so far
	float emitDirection;
	float emitAngle;
	float emitTimer;
	float emitStepTimer;
	RadiationEmitterType emitType;
	int emitStep;
};
struct ElementChunkData
{
	struct ExternalAdjuster
	{
		float temperature;
		float heatCapacity;
		float thermalConductivity;
	};
	int   cell;
	float temperature;
	float heatCapacity;
	float thermalConductivity;
	float maxEnergyTransferScaleFactor;
	float groundTransferScale;
	float highTemp;
	float lowTemp;
	ElementChunkData::ExternalAdjuster adjuster;
};
struct BuildingHeatExchangeData
{
	float temperature;
	float overheatTemperature;
	float operatingKilowatts;
	float totalHeatCapacity;
	float perCellHeatCapacity;
	float thermalConductivity;
	float highTemp;
	Vector2<int> simMin;
	Vector2<int> simMax;
};
struct InContactBuildingData
{
	int inContactBuildingHandler;
	int cellsInContact;
};
struct BuildingToBuildingHeatExchangeData
{
	int heatExchange_handle;
	std::vector<InContactBuildingData> inContactBuildingsData;
};
struct DiseaseEmitterData
{
	int		cell;
	uint8_t	diseaseIdx;
	uint8_t	range;
	int		emitCount;
	float	emitInterval;
	float	elapsedTime;
};
struct DiseaseConsumerData
{
	float consumptionRate;
	uint8_t maxDepth;
	Vector2<uint16_t> cell;
	uint16_t propertyMask;
};

//----- CompactedVector message data -----
struct BuildingHeatExchangeMessageData
{
	uint16_t elemIdx;
	//uint8_t pad[2];
	float mass;
	float temperature;
	float thermalConductivity;
	float overheatTemperature;
	float operatingKilowatts;
	int minX;
	int minY;
	int maxX;
	int maxY;
};

//----- CompactedVector Detail data -----
struct AddBuildingHeatExchangeMessage
{
	int callbackIdx;
	BuildingHeatExchangeMessageData data;
};
struct ModifyBuildingHeatExchangeMessage
{
	int handle;
	BuildingHeatExchangeMessageData data;
};
struct RemoveBuildingHeatExchangeMessage
{
	int handle;
	int callbackIdx;
};
struct RegisterBuildingToBuildingHeatExchangeMessage
{
	int callbackIdx;
	int heatExchange_handle;
};
struct AddInContactBuildingToBuildingToBuildingHeatExchangeMessage
{
	//AddBuildingToBuildingHeatExchangeMessageData data;
	int self_handle;
	int buildingInContact;
	int cellsInContact;
};
struct RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage
{
	// RemoveBuildingToBuildingHeatExchangeMessageData data;
	int self_handle;
	int buildingInContact;
};
struct RemoveBuildingToBuildingHeatExchangeMessage
{
	int callbackIdx;
	int handle;
};
struct AddElementChunkMessage
{
	int   gameCell;
	int   callbackIdx;
	float mass;
	float temperature;
	float surfaceArea;
	float thickness;
	float groundTransferScale;
	uint16_t elementIdx;
	//uint8_t pad[2];
};
struct ModifyElementChunkMessage
{
	int   handle;
	float temperature;
	float heatCapacity;
};
struct ModifyElementChunkAdjusterMessage
{
	int   handle;
	float temperature;
	float heatCapacity;
	float thermalConductivity;
};
struct ModifyElementChunkEnergyMessage
{
	int   handle;
	float deltaKJ;
};
struct RemoveElementChunkMessage
{
	int handle;
	int callbackIdx;
};
struct AddElementConsumerMessage
{
	int gameCell;
	int callbackIdx;
	uint8_t  radius;
	uint8_t  configuration;
	uint16_t elementIdx;
};
struct ModifyElementConsumerMessage
{
	int   handle;
	int   gameCell;
	float consumptionRate;
};
struct RemoveElementConsumerMessage
{
	int handle;
	int callbackIdx;
};
struct AddElementEmitterMessage
{
	float maxPressure;
	int   callbackIdx;
	int   onBlockedCBIdx;
	int   onUnblockedCBIdx;
};
struct ModifyElementEmitterMessage
{
	int   handle;
	int   gameCell;
	float emitInterval;
	float emitMass;
	float emitTemperature;
	float maxPressure;
	int   diseaseCount;
	uint16_t elementIdx;
	uint8_t  maxDepth;
	uint8_t  diseaseIdx;
};
struct RemoveElementEmitterMessage
{
	int handle;
	int callbackIdx;
};
struct AddRadiationEmitterMessage
{
	int callbackIdx;
	int gameCell;
	uint16_t radiusX;
	uint16_t radiusY;
	float emitRads;
	float emitRate;
	float emitSpeed;
	float emitDirection;
	float emitAngle;
	int   emitType;
};
struct ModifyRadiationEmitterMessage
{
	int handle;
	int gameCell;
	int callbackIdx;
	uint16_t radiusX;
	uint16_t radiusY;
	float emitRads;
	float emitRate;
	float emitSpeed;
	float emitDirection;
	float emitAngle;
	int   emitType;
};
struct RemoveRadiationEmitterMessage
{
	int handle;
	int callbackIdx;
};
struct AddDiseaseEmitterMessage
{
	int callbackIdx;
};
struct ModifyDiseaseEmitterMessage
{
	int     handle;
	int     gameCell;
	uint8_t diseaseIdx;
	uint8_t maxDepth;
	//uint8_t pad[2];
	float   emitInterval;
	int     emitCount;
};
struct RemoveDiseaseEmitterMessage
{
	int handle;
	int callbackIdx;
};
struct AddDiseaseConsumerMessage
{
	int callbackIdx;
	int onBlockedCBIdx;
	int onUnblockedCBIdx;
};
struct ModifyDiseaseConsumerMessage
{
	int     handle;
	int     gameCell;
	uint8_t diseaseIdx;
	uint8_t maxDepth;
	//uint8_t pad[2];
	float   emitInterval;
	int     emitCount;
};
struct RemoveDiseaseConsumerMessage
{
	int handle;
	int callbackIdx;
};

//----- SimEvents struct -----
struct SubstanceChangeInfo
{
	int	     cellIdx;
	uint16_t oldElementIdx;
	uint16_t newElementIdx;
};
struct SpawnFallingLiquidInfo
{
	int      cellIdx;
	uint16_t elementIdx;
	uint8_t  diseaseIdx;
	//	uint8_t  pad;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct SpawnOreInfo
{
	int      cellIdx;
	uint16_t elemIdx;
	uint8_t  diseaseIdx;
	//	uint8_t  pad;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct UnstableCellInfo
{
	int      cellIdx;
	uint16_t elemIdx;
	uint8_t  fallingInfo;
	uint8_t  diseaseIdx;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct MeltedInfo
{
	Handle handle;
};
struct CellMeltedInfo
{
	uint32_t gameCell;
};
struct CallbackInfo
{
	int callbackIdx;
};
struct WorldDamageInfo
{
	int cellIdx;
	int damageSourceCellIdx;
};
struct MassConsumedCallback
{
	int	     callbackIdx;
	uint16_t elemIdx;
	uint8_t  diseaseIdx;
	//	uint8_t  pad;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct RadiationConsumedCallback
{
	int   callbackIdx;
	int   gameCell;
	float radiation;
};
struct MassEmittedCallback
{
	int      callbackIdx;
	uint16_t elemIdx;
	uint8_t  emitted;
	uint8_t  diseaseIdx;
	float    mass;
	float    temperature;
	int      diseaseCount;
};
struct DiseaseConsumedCallback
{
	int     callbackIdx;
	uint8_t diseaseIdx;
	//	uint8_t pad[3];
	int     diseaseCount;
};
struct SpawnFXInfo
{
	int   cellIdx;
	int   fxID;
	float rotation;
};
struct ComponentStateChangedMessage
{
	int callbackIdx;
	int simHandle;
};

//----- define SimComponent class -----
class SimComponent
{
public:
	//SimComponent_vtbl* __vftable /*VFT*/;
	virtual void Update(float, SimData*, Region*) {}
	virtual void UpdateDataListOnly(SimData*) {}
};
class BuildingHeatExchange : public CompactedVector<BuildingHeatExchangeData>, public SimComponent
{
public:
	std::vector<BuildingTemperatureInfo> temperatureInfo;

	Handle Register(          SimData* simData, AddBuildingHeatExchangeMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyBuildingHeatExchangeMessage* msg);
	void Update    (float dt, SimData* simData, Region* region);
	void UpdateDataListOnly(  SimData* simData);
	void ChangeBuildingTemperature      (int buildingHandle, float newTemperature);
	Handle GetSimHandleForBuildingHandle(int buildingHandle);
};
class BuildingToBuildingHeatExchange : public CompactedVector<BuildingToBuildingHeatExchangeData>, public SimComponent
{
public:
	Handle Register  (          SimData* simData, RegisterBuildingToBuildingHeatExchangeMessage* msg);
	void   Unregister(int handle);
	void   Add       (          SimData* simData, AddInContactBuildingToBuildingToBuildingHeatExchangeMessage* msg);
	void   Remove    (          SimData* simData, RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage* msg);
	void   Update    (float dt, SimData* simData, Region* region);
};
class ElementChunk : public CompactedVector<ElementChunkData>, public SimComponent
{
public:
	std::vector<ElementChunkInfo> chunkInfo;

	Handle Register    (          SimData* simData, AddElementChunkMessage* msg);
	void Unregister    (int handle); // RadiationEmitter::Unregister
	void Modify        (          SimData* simData, ModifyElementChunkMessage* msg);
	void ModifyAdjuster(          SimData* simData, ModifyElementChunkAdjusterMessage* msg);
	void ModifyEnergy  (          SimData* simData, ModifyElementChunkEnergyMessage* msg);
	void Update        (float dt, SimData* simData, Region* region);
	void UpdateDataListOnly      (SimData* simData);
};
class ElementConsumer : public CompactedVector<ElementConsumerData>, public SimComponent
{
public:
	std::vector<int>          visited;
	std::queue<FloodFillInfo> next;
	std::vector<int>          reachableCells;
	std::vector<ConsumedMassInfo> consumedMassInfo;

	Handle Register(          SimData* simData, AddElementConsumerMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyElementConsumerMessage* msg);
	void Update    (float dt, SimData* simData, Region* region);
};
class ElementEmitter : public CompactedVector<ElementEmitterData>, public SimComponent
{
public:
	std::vector<EmittedMassInfo> emittedMassInfo;
	std::vector<int> visitedCells;
	std::vector<int> reachableCells;
	std::queue<FloodFillInfo> next;

	Handle Register(          SimData* simData, AddElementEmitterMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyElementEmitterMessage* msg);
	void Update    (float dt, SimData* simData, Region* region);
	void UpdateDataListOnly(SimData* simData);
	EmittedMassInfo Emit   (ElementEmitterData* data, int cell_idx,            float tempEmit, SimData* simData);
	EmittedMassInfo TryEmit(ElementEmitterData* data, std::vector<int>& cells, int   offset,   SimData* simData);
};
class RadiationEmitter : CompactedVector<RadiationEmitterData>, public SimComponent
{
public:
	Handle Register(          SimData* simData, AddRadiationEmitterMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyRadiationEmitterMessage* msg);
	void Update    (float dt, SimData* simData, Region* region);
};
class DiseaseEmitter : CompactedVector<DiseaseEmitterData>, public SimComponent
{
public:
	std::vector<DiseaseEmittedInfo> emittedInfo;
	std::queue<FloodFillInfo> next;
	std::vector<int> visitedCells;
	std::vector<int> reachableCells;

	Handle Register(          SimData* simData, AddDiseaseEmitterMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyDiseaseEmitterMessage* msg);
	void Update    (float dt, SimData* simData, Region* region);
	void UpdateDataListOnly(  SimData* simData);
};
class DiseaseConsumer : CompactedVector<DiseaseConsumerData>, public SimComponent
{
public:
	std::vector<DiseaseConsumedInfo> consumedInfo;

	Handle Register(          SimData* simData, AddDiseaseConsumerMessage* msg);
	void Unregister(int handle);
	void Modify    (          SimData* simData, ModifyDiseaseConsumerMessage* msg);
};

//----- define Sim class -----
class SimData
{
	struct Timers
	{
		uint8_t stableCellTicks : 5;
	};

public:
	struct WorldOffsetData
	{
		int offsetX;
		int offsetY;
		int width;
		int height;
	};

	int width;
	int height;
	int numGameCells;
	uint8_t  savedOptions;
	uint16_t voidElementIdx;
	uint16_t vacuumElementIdx;
	uint16_t unobtaniumElementIdx;
	std::unique_ptr<CellSOA>	cells;
	std::unique_ptr<CellSOA>	updatedCells;
	std::unique_ptr<uint8_t[]>	worldZones;
	std::vector<float>			cosmicRadiationOcclusion;
	std::unique_ptr<float[]>	accumulatedFlow;
	std::unique_ptr<Vector4<float>[], std::default_delete<Vector4<float>[]>> flow; // mass move to [x:-noffset|y:+noffset|z:up|w:down]
	std::unique_ptr<float[]>	biomeTemperature;
	std::unique_ptr<SimData::Timers[]> timers;
	std::vector<ActiveRegion>  activeRegions;
	std::unique_ptr<uint8_t[]> visibleGrid;
	uint16_t tickCount;
	std::unique_ptr<SimEvents> simEvents;
	DebugProperties            debugProperties;
	std::vector<SimComponent*> components;
	ElementConsumer            elementConsumer;
	ElementEmitter             elementEmitter;
	RadiationEmitter           radiationEmitter;
	ElementChunk               elementChunk;
	BuildingHeatExchange       buildingHeatExchange;
	BuildingToBuildingHeatExchange buildingToBuildingHeatExchange;
	DiseaseEmitter  diseaseEmitter;
	DiseaseConsumer diseaseConsumer;
	bool initSettleThermalBoundaries;
	bool headless;
	bool radiationEnabled;
	float RADIATION_LINGER_RATE;
	float RADIATION_MAX_MASS;
	float RADIATION_BASE_WEIGHT;
	float RADIATION_DENSITY_WEIGHT;
	float RADIATION_CONSTRUCTED_FACTOR;
	std::vector<SimData::WorldOffsetData> worlds;
	uint32_t randomSeed;

	explicit SimData(int sim_width, int sim_height, bool radiation_enabled, bool headless);
	void	ApplySaveSettings(uint8_t newSaveOptions);
	void	ClearCell(CellAccessor* cell);
	uint8_t	GetStableTicksRemaining(uint64_t sim_cell);
	void	SetBoundary(int cell);
	void	InitializeBoundary();
	void	SettleThermalBoundaries(CellSOA * src_cells, CellSOA * dest_cells);
	void	UpdateComponents(float dt, Region * region);
	void	UpdateComponentsDataListOnly();
	static	int FreeGridCells(BinaryBufferReader * reader);
	static  int InitializeFromCells(BinaryBufferReader* reader);
	static  int	ResizeAndInitializeVacuumCells(BinaryBufferReader* reader);
	static	int SetWorldZones(BinaryBufferReader* reader);
};

class SimEvents
{
public:
	std::vector<SubstanceChangeInfo>	substanceChangeInfo;
	std::vector<SpawnFallingLiquidInfo>	spawnLiquidInfo;
	std::vector<SpawnOreInfo>			spawnOreInfo;
	std::vector<UnstableCellInfo>		unstableCellInfo;
	std::vector<MeltedInfo>	elementChunkMeltedInfo;
	std::vector<MeltedInfo>	buildingMeltedInfo;
	std::vector<MeltedInfo>	buildingOverheatInfo;
	std::vector<MeltedInfo>	buildingNoLongerOverheatedInfo;
	std::vector<CellMeltedInfo>				cellMeltedInfo;
	std::vector<CallbackInfo>				callbackInfo;
	std::vector<WorldDamageInfo>			worldDamageInfo;
	std::vector<MassConsumedCallback>		massConsumedCallbacks;
	std::vector<RadiationConsumedCallback>	radiationConsumedCallbacks;
	std::vector<MassEmittedCallback>		massEmittedCallbacks;
	std::vector<DiseaseConsumedCallback>	diseaseConsumedCallbacks;
	std::vector<SpawnFXInfo>				spawnFXInfo;
	std::vector<ComponentStateChangedMessage>	componentStateChangedMessages;
	std::vector<SpawnOreInfo>					digInfo;

	void ChangeSubstance   (const SimData* simData, uint64_t sim_cell);
	bool SpawnOre          (const SimData* simData, uint64_t sim_cell, uint16_t elemIdx, float mass, float temperature, uint8_t disease_idx, int disease_count, bool spawn_even_if_hidden);
	bool SpawnFallingLiquid(const SimData* simData, uint64_t sim_cell, uint16_t elemIdx, float mass, float temperature, uint8_t disease_idx, int disease_count, bool spawn_even_if_hidden);
	void SpawnFX           (const SimData* simData, uint64_t sim_cell, int fxid, float rotation);
};

//----- global Function -----
extern void  FloodRemoved(SimData* simData, float* mass_to_remove, int16_t remove_elem_idx, ConsumedMassInfo* removed_info,
	uint16_t x, uint16_t y, int max_depth, std::vector<int>* visited, std::queue<FloodFillInfo>* next);
#endif