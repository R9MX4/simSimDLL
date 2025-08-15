#pragma once
#ifndef CLASS_BASE_H
#define CLASS_BASE_H

#include <mutex>
#include <threads.h>
#include "ClassDisease.h"
#include "ClassSim.h"
#include "GameData.h"
//----- define struct -----
struct DigPoint
{
	int  gameCell;
	int  callbackIdx;
	bool skipEvent;
};
struct CellModification
{
	int   gameCell;
	int   callbackIdx;
	float temperature;
	float mass;
	int   diseaseCount;
	uint16_t elementIdx;
	uint8_t  replaceType; // 0:None, 1:Replace, 2:ReplaceAndDisplace
	uint8_t  diseaseIdx;
	uint8_t  addSubType;
};
struct CellEnergyModification
{
	int   gameCell;
	float kilojoules;
	float maxTemperature;
	int   id;
};
struct PipeChange
{
	int     gameCell;
	uint8_t layer;
	//uint8_t pad[3];
	float   mass;
	float   temperature;
	int     elementHash;
};
struct SetCellFloatValue
{
	int   gameCell;
	float value;
};
struct CellPropertiesChange
{
	int     gameCell;
	uint8_t properties;
	uint8_t set;
	//uint8_t pad[2];
};
struct ConsumeDisease
{
	int   gameCell;
	int   callbackIdx;
	float percentToConsume;
	int   maxToConsume;
};
struct CellDiseaseModification
{
	int     gameCell;
	uint8_t diseaseIdx;
	int     diseaseCount;
};
struct CellRadiationModification
{
	int   gameCell;
	float radiationDelta;
	int   callbackIdx;
};
struct RadiationParamsModification
{
	int   type;
	float value;
};
struct CellWorldZoneModification
{
	int     gameCell;
	uint8_t zoneID;
};
struct MassConsumption
{
	int      gameCell;
	int      callbackIdx;
	float    mass;
	uint16_t elementIdx;
	uint8_t  radius;
};
struct MassEmission
{
	int   gameCell;
	int   callbackIdx;
	float mass;
	float temperature;
	int   diseaseCount;
	uint16_t elementIdx;
	uint8_t  diseaseIdx;
};
struct ModifyBuildingEnergyMessage
{
	Handle handle;
	float  deltaKJ;
	float  minTemperature;
	float  maxTemperature;
};
struct MoveElementChunkMessage
{
	int handle;
	int gameCell;
};
struct CreateElementsTableMessage
{
	int      numElements;
	Element* elements;
};

struct Task
{
	virtual ~Task() {};
	virtual void InternalDoTask() {};
};

template <typename T>
struct UpdateCellTask : Task
{
	const T* src     = 0;
	T* dest          = 0;
	int num_rows     = 0;
	int src_y_start  = 0;
	int stride       = 0;
	int dest_y_start = 0;
	int width        = 0;

	void InternalDoTask();
};

struct UpdateCellSOATask : UpdateCellTask<CellSOA>
{
	void InternalDoTask();
};

//----- define Main struct and class -----
struct NewGameFrame
{
	float  elapsedSeconds;
	Region activeRegion;
	float  currentSunlightIntensity;
	float  currentCosmicRadiationIntensity;
};

struct SimFrameInfo
{
	template <typename T, typename U, typename V>
	struct ComponentMessages {
		std::vector<T> adds;
		std::vector<U> modifies;
		std::vector<V> removes;
	};

	float elapsedSeconds;
	std::vector<DigPoint> digPoints;
	std::vector<CellModification> cellModifications;
	std::vector<CellEnergyModification> cellEnergyModifications;
	std::vector<PipeChange> pipeChanges;
	std::vector<SetCellFloatValue> setInsulationValues;
	std::vector<SetCellFloatValue> setStrengthValues;
	std::vector<CellPropertiesChange> setCellProperties;
	std::vector<CellPropertiesChange> clearCellProperties;
	std::vector<ConsumeDisease> consumeDisease;
	std::vector<CellDiseaseModification> cellDiseaseModifications;
	std::vector<CellRadiationModification> cellRadiationModifications;
	std::vector<RadiationParamsModification> radiationParamsModifications;
	std::vector<CellWorldZoneModification> cellWorldZoneModifications;
	std::vector<MassConsumption> massConsumptionMessages;
	std::vector<MassEmission> massEmissionMessages;
	SimFrameInfo::ComponentMessages<AddBuildingHeatExchangeMessage, ModifyBuildingHeatExchangeMessage, RemoveBuildingHeatExchangeMessage> buildingHeatExchangeMessages;
	std::vector<ModifyBuildingEnergyMessage> modifyBuildingEnergyMessages;
	SimFrameInfo::ComponentMessages<RegisterBuildingToBuildingHeatExchangeMessage, RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage, RemoveBuildingToBuildingHeatExchangeMessage> buildingToBuildingHeatExchangeMessages;
	std::vector<AddInContactBuildingToBuildingToBuildingHeatExchangeMessage> addBuildingInContactToBuildingToBuildingHeatExchangeMessages;
	SimFrameInfo::ComponentMessages<AddElementChunkMessage, ModifyElementChunkMessage, RemoveElementChunkMessage> elementChunkMessages;
	std::vector<MoveElementChunkMessage> moveElementChunkMessages;
	std::vector<ModifyElementChunkEnergyMessage> modifyElementChunkEnergyMessages;
	std::vector<ModifyElementChunkAdjusterMessage> modifyElementChunkAdjusterMessages;
	SimFrameInfo::ComponentMessages<AddElementConsumerMessage, ModifyElementConsumerMessage, RemoveElementConsumerMessage> elementConsumerMessages;
	SimFrameInfo::ComponentMessages<AddElementEmitterMessage, ModifyElementEmitterMessage, RemoveElementEmitterMessage> elementEmitterMessages;
	SimFrameInfo::ComponentMessages<AddDiseaseEmitterMessage, ModifyDiseaseEmitterMessage, RemoveDiseaseEmitterMessage> diseaseEmitterMessages;
	SimFrameInfo::ComponentMessages<AddDiseaseConsumerMessage, ModifyDiseaseConsumerMessage, RemoveDiseaseConsumerMessage> diseaseConsumerMessages;
	SimFrameInfo::ComponentMessages<AddRadiationEmitterMessage, ModifyRadiationEmitterMessage, RemoveRadiationEmitterMessage> radiationEmitterMessages;
	DebugProperties debugProperties;
};

class SimFrameManager
{
public:
	mtx_t         frameMutex;
#ifdef __DEBUGED__
	mtx_t         regionMutex; // New Add
#endif
	SimFrameInfo* currentFrame;
	std::vector<SimFrameInfo*> framePool;
	std::vector<SimFrameInfo*> queuedFrames;
	std::vector<SimFrameInfo*> activeFrames;
	std::vector<SimFrameInfo*> processedFrames;
	NewGameFrame               activeRegion;
	float                      elapsedSeconds;
	std::vector<ActiveRegion>  activeRegions;
	int                        numFramesProcessed;

	explicit SimFrameManager();
	int64_t BeginFrameProcessing();
	void EndFrameProcessing();
	bool HandleMessage(Hashes::SimMessageHashes msg_id, int msg_length, char* msg);
	void NewFrame();
	void ProcessFrame(SimFrameInfo* frame, SimData* simData);
	float ProcessNextFrame(SimData* simData);
};

struct SimUpdateTask : Task
{
	SimData*     simData;
	SimEvents*   simEvents;
	Vector2<int> start;
	Vector2<int> end;
};

struct SimUpdateTemperatureTask : SimUpdateTask
{
	void InternalDoTask();
};

struct ParallelTaskQueue
{
	cnd_t mWorkCompleteCondition;
	cnd_t mTaskCondition;
	mtx_t mMutex;
	std::queue<Task*> mTasks;
	uint64_t mTaskCount;
	std::vector<std::thread> mWorkers;
	bool mShuttingDown;

	explicit ParallelTaskQueue(uint64_t threadCount);
	~ParallelTaskQueue();
	void WorkerFunc();
};

class SimBase
{
	struct CellInfo
	{
		int      cell;
		bool     impermeable;
		uint16_t element_index;
		const ElementPressureData* pressure;
		uint8_t  element_state; //state -> Element.State: Vacuum=0; Gas=1; Liquid=2; Solid=3;
	};

	//SimBase_vtbl* __vftable /*VFT*/;
	ParallelTaskQueue* taskQueue;
	std::vector<SimUpdateTask*> temperatureTasks;
	std::vector<SimEvents*>     simEvents;
	std::vector<UpdateCellTask<Vector4<float>>*> copyFlowTasks;
	std::vector<UpdateCellSOATask*>              copyToGameTasks;

public:
	SimFrameManager simFrameManager;
	float elapsedSeconds;

	explicit SimBase();
	~SimBase();
	// void DestroyTasks();
	void UpdateData(SimData* simData);
	void InitializeUpdateTasks();
	void ConsolidateEvents(SimData* simData);
	void CopyUpdatedCellsToCells(SimData* simData);
	void CopySimDataToGame(SimData* simData, GameData* new_game_data, const GameData* old_game_data, int num_frames_processed);
};

class FrameSync
{
public:
	std::unique_ptr<GameData> mSimData;
	std::unique_ptr<GameData> mGameData;
	bool  mSimReadyToSwap = false;
	bool  mInitialized    = false;
	mtx_t mMutex;
	cnd_t mGameCond;
	cnd_t mSimCond;
	mtx_t mSimMutex;

	void clear();
	void initGameData(int game_width, int game_height, int sim_width, int sim_height, CellSOA* updatedCells);
	//void GameSync(void(*syncFunction)(GameData*, GameData*));
	void GameSync_CleanUp();
	void GameSync_PrepareGameData(BinaryBufferReader* reader);
	void SimSync();
};

class Sim : public SimBase
{
public:
	FrameSync* frameSync;
	// Thread
	std::string mName;
	std::thread mThread;
	std::atomic<bool> mExitRequested;

	explicit Sim(FrameSync* fs);
	void Main();
	void Start();
};
#endif