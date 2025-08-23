#pragma once
#ifndef PCH_H
#define PCH_H

#define __DEBUGED__
#define __DEBUG_PRINT__
//#undef __DEBUG_PRINT__
#define __SIMDLL_PLUS__
//#undef __SIMDLL_PLUS__

#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <memory>

//----- typedef -----
typedef int Handle;
//----- global const -----
const float SIM_MAX_TEMPERATURE = 10000.0f;
const float SIM_MAX_RADIATION   = 9000000.0f;
const float SIM_MIN_RADIATION   = 0.0f;
//----- global variable -----
extern FILE*  gLogger;
extern FILE*  gLogger2;
extern int    gLogLevel;
extern double gGasDisplace;

extern void PRINT_TIME(FILE* _Logger);
extern void INIT_TIMER();
extern float MAX_F(float a, float b);
extern float MIN_F(float a, float b);
extern float CLAMP_F(float val, float max, float min);
//----- pre define -----
#ifdef __DEBUG_PRINT__
#define LOGGER_INIT()       fopen_s(&gLogger,  "simlog.txt", "w+");
#define LOGGER_INIT2()      fopen_s(&gLogger2, "simlogPara.txt", "w+");
#define LOGGER_PRINT(...)  (PRINT_TIME(gLogger ), fprintf(gLogger , __VA_ARGS__), std::fflush(gLogger ))
#define LOGGER_PRINT2(...) (PRINT_TIME(gLogger2), fprintf(gLogger2, __VA_ARGS__), std::fflush(gLogger2))
#define LOGGER_LEVEL(_LEVEL, ...)  if(gLogLevel >= _LEVEL) LOGGER_PRINT2(__VA_ARGS__);
#else
#define LOGGER_INIT()              {;}
#define LOGGER_INIT2()             {;}
#define LOGGER_PRINT(...)          (0)
#define LOGGER_PRINT2(...)         (0)
#define LOGGER_LEVEL(_LEVEL, ...)  {;}
#endif

#define ASSERT_TEXT(_text) LOGGER_PRINT("%s! File:%s Func:%s Line:%d\n", _text, __FILE__, __func__, __LINE__)
#define ASSERT_TEMP(_mass, _temperature) \
        if ((_temperature <= 0.0 && _mass > 0.0) || _temperature > SIM_MAX_TEMPERATURE || std::isnan(_temperature)) {ASSERT_TEXT("Assert Temperature failed");}
#define CELL_G2S(_gameCell, _width) (_gameCell % (_width - 2) +  _width      * (_gameCell / (_width - 2) + 1) + 1)
#define CELL_S2G(_simCell,  _width) (_simCell  %  _width      + (_width - 2) * (_simCell  /  _width      - 1) - 1)
#define CELL_AVAILABLE(_cell, _simData) (_cell >= 0 && _cell < _simData->numGameCells)
#define CELL_VISIABLE( _cell, _simData) (_simData->visibleGrid[_cell] || _simData->debugProperties.isDebugEditing)
#define HANDLE_AVAILABLE(_handle, _version) \
		((_handle & 0xFFFFFF) < 0x100000 && _version[_handle & 0xFFFFFF] == (uint8_t)(_handle >> 24))
//#define RAND_INT(_seed, _min,_max) (srand(_seed), _seed=rand(), _seed%(_max-_min)+_min)
//#define RAND_FLOAT_01(_seed)       (srand(_seed), _seed=rand(), _seed/double(RAND_MAX))
#define RAND_INT(_seed, _min,_max) (_seed=214013 * _seed + 2531011,  _seed%(_max-_min)+_min)
#define RAND_FLOAT_01(_seed)       (_seed=214013 * _seed + 2531011, (_seed&RAND_MAX)/double(RAND_MAX))
#define GET_STATE(_cells, _cell)   (gElements[_cells->elementIdx[_cell]].state & 3)

//----- template -----
template <typename T>
struct Vector2 { T x;	T y; };

template <typename T>
struct Vector3 { T x;	T y;	T z; };

template <typename T>
struct Vector4 { T x;	T y;	T z;	T w; };

template <typename T>
struct CompactedVector {
	struct HandleVector
	{
		std::vector<Handle>  freeHandles;
		std::vector<int>     items;
		std::vector<uint8_t> versions;
	};
	std::vector<T>      data;
	std::vector<Handle> dataHandleIndices;
	HandleVector        handles;

	Handle Free   (int handle);
	Handle AddData(            T* new_data);
	void   SetData(int handle, T* new_data);
	T*     GetData(int handle);
};

//----- enum -----
enum CellProperties : __int32
{
	None				= 0x0,
	GasImpermeable		= 0x1,
	LiquidImpermeable	= 0x2,
	SolidImpermeable	= 0x4,
	Unbreakable			= 0x8,
	Transparent			= 0x10,
	Opaque				= 0x20,
	NotifyOnMelt		= 0x40,
	ConstructedTile		= 0x80,
};
enum RadiationEmitterType : __int32
{
	Constant		= 0x0,
	Pulsing			= 0x1,
	PulsingAveraged	= 0x2,
	SimplePulse		= 0x3,
	RadialBeams		= 0x4,
	Attractor		= 0x5,
};
enum AddSolidMassSubType : uint8_t
{
	DoVerticalDisplacement = 0x0,
	OnlyIfSameElement      = 0x1,
};
namespace Hashes {
	enum SimMessageHashes : int
	{
		Elements_CreateTable					=  1108437482,
		Elements_CreateInteractions				=  -930289787,
		SetWorldZones							=  -457308393,
		ModifyCellWorldZone						=  -449718014,
		Disease_CreateTable						=   825301935,
		Load									=  -672538170,
		Start									=  -931446686,
		AllocateCells							=  1092408308,
		ClearUnoccupiedCells					= -1836204275,
		DefineWorldOffsets						=  -895846551,
		PrepareGameData							=  1078620451,
		SimData_InitializeFromCells				=  2062421945,
		SimData_ResizeAndInitializeVacuumCells	=  -752676153,
		SimData_FreeCells						= -1167792921,
		SimFrameManager_NewGameFrame			=  -775326397,
		Dig										=   833038498,
		ModifyCell								= -1252920804,
		ModifyCellEnergy						=   818320644,
		SetInsulationValue						=  -898773121,
		SetStrengthValue						=  1593243982,
		SetVisibleCells							=  -563057023,
		ChangeCellProperties					=  -469311643,
		AddBuildingHeatExchange					=  1739021608,
		ModifyBuildingHeatExchange				=  1818001569,
		ModifyBuildingEnergy					= -1348791658,
		RemoveBuildingHeatExchange				=  -456116629,
		AddBuildingToBuildingHeatExchange							= -1338718217,
		AddInContactBuildingToBuildingToBuildingHeatExchange		= -1586724321,
		RemoveBuildingInContactFromBuildingToBuildingHeatExchange	= -1993857213,
		RemoveBuildingToBuildingHeatExchange						=   697100730,
		SetDebugProperties				= -1683118492,
		MassConsumption					=  1727657959,
		MassEmission					=   797274363,
		AddElementConsumer				=  2024405073,
		RemoveElementConsumer			=   894417742,
		SetElementConsumerData			=  1575539738,
		AddElementEmitter				=  -505471181,
		ModifyElementEmitter			=   403589164,
		RemoveElementEmitter			= -1524118282,
		AddElementChunk					=  1445724082,
		RemoveElementChunk				=  -912908555,
		SetElementChunkData				=  -435115907,
		MoveElementChunk				=  -374911358,
		ModifyElementChunkEnergy		=  1020555667,
		ModifyChunkTemperatureAdjuster	= -1387601379,
		AddDiseaseEmitter				=  1486783027,
		ModifyDiseaseEmitter			= -1899123924,
		RemoveDiseaseEmitter			=   468135926,
		AddDiseaseConsumer				=   348345681,
		ModifyDiseaseConsumer			= -1822987624,
		RemoveDiseaseConsumer			=  -781641650,
		ConsumeDisease					= -1019841536,
		CellDiseaseModification			= -1853671274,
		ToggleProfiler					=  -409964931,
		SetSavedOptions					=  1154135737,
		CellRadiationModification		= -1914877797,
		RadiationSickness				=  -727746602,
		AddRadiationEmitter				= -1505895314,
		ModifyRadiationEmitter			=  -503965465,
		RemoveRadiationEmitter			=  -704259919,
		RadiationParamsModification		=   377112707,
		SetLiquidCompression            =  1396974888
	};

	enum ElementInteractionHashes : int
	{
		GasObliteration = 1747383933,
		LiquidObliteration = 764054208,
		LiquidConversion = 1537054354
	};
}
struct ElementState {
	enum State : __int32
	{
		Vacuum	= 0x0,
		Gas		= 0x1,
		Liquid	= 0x2,
		Solid	= 0x3
	};
	enum Masks : __int32
	{
		StateMask	= 0x3,
		FlagMask	= 0x14
	};
};
#endif
