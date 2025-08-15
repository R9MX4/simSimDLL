#pragma once
#ifndef GLOBAL_H
#define GLOBAL_H

#include "pch.h"
//----- declar struct -----
struct GameDataUpdate;

//----- declar class -----
class FrameSync;
class Sim;
class SimData;
class SimEvents;
class Disease;
class CellSOA;

//----- define struct -----
//state -> Element.State: Vacuum=0; Gas=1; Liquid=2; Solid=3;
//                        Unbreakable=4; Unstable=8; TemperatureInsulated=16;
struct ElementLiquidData
{
	uint8_t state;
	float flow;
	float viscosity;
	float minHorizontalFlow;
	float minVerticalFlow;
	float maxMass;
	float compression = 1.01f; // New Add, layer compression
};
struct ElementPressureData
{
	uint8_t	state;
	float	flow;
	float	molarVolume; // New Add
};
struct PhysicsData
{
	float temperature;
	float mass;
	float pressure;
};
struct Element
{
	int id;
	uint16_t lowTempTransitionIdx;
	uint16_t highTempTransitionIdx;
	uint16_t elementsTableIdx;
	uint8_t  state;
	uint8_t  pack0;
	float specificHeatCapacity;
	float thermalConductivity;
	float molarMass;
	float solidSurfaceAreaMultiplier;
	float liquidSurfaceAreaMultiplier;
	float gasSurfaceAreaMultiplier;
	float flow;
	float viscosity;
	float minHorizontalFlow;
	float minVerticalFlow;
	float maxMass;
	float lowTemp;
	float highTemp;
	float strength;
	int   lowTempTransitionOreID;
	float lowTempTransitionOreMassConversion;
	int   highTempTransitionOreID;
	float highTempTransitionOreMassConversion;
	uint16_t sublimateIndex;
	uint16_t convertIndex;
	uint32_t colour;
	int   sublimateFX;
	float sublimateRate;
	float sublimateEfficiency;
	float sublimateProbability;
	float offGasPercentage;
	float lightAbsorptionFactor;
	float radiationAbsorptionFactor;
	float radiationPer1000Mass;
	PhysicsData defaultValues;
};
struct ElementLightAbsorptionData
{
	float factor;
	float massScale;
};
struct ElementRadiationData
{
	float factor;
	float rads_per_1000;
};
struct ElementPostProcessData
{
	uint16_t sublimateIndex;
	uint16_t convertIndex;
	uint8_t  state;
	float sublimateRate;
	float sublimateEfficiency;
	float sublimateProbability;
	float offGasPercentage;
	float molarMass;
	float strength;
	float maxMass;
	float minHorizontalFlow;
	CellProperties sublimateFX; // SpawnFXHashes::Type sublimateFX;
	float maxCompression = 1.35f; // New Add, stacked compression, 1.01x30 layers
	float compression    = 1.01f; // New Add, layer compression
};
struct ElementStateData
{
	uint8_t state;
};
struct ElementPropertyTextureData
{
	uint32_t colour;
	uint8_t state;
};
struct ElementTemperatureData
{
	uint8_t  state;
	uint16_t lowTempTransitionIdx;
	uint16_t highTempTransitionIdx;
	uint16_t lowTempTransitionOreIdx;
	uint16_t highTempTransitionOreIdx;
	float specificHeatCapacity;
	float thermalConductivity;
	float defaultMass;
	float massAreaScale;
	float gasSurfaceAreaMultiplier;
	float liquidSurfaceAreaMultiplier;
	float solidSurfaceAreaMultiplier;
	float lowTemp;
	float highTemp;
	float lowTempTransitionOreMassConversion;
	float highTempTransitionOreMassConversion;
};
struct GasObliteration
{
	uint16_t elemIdx1;
	uint16_t elemIdx2;
	uint16_t elemResultIdx;
	float interactionProbability;
	float minMass;
	float elem1MassDestructionPercent;
	float elem2MassRequiredMultiplier;
	float elemResultMassCreationMultiplier;
};
struct LiquidConversion
{
	uint16_t elemIdx1;
	uint16_t elemIdx2;
	float massRatio;
	float probability;
};
struct LiquidObliteration
{
	uint16_t elemIdx1;
	uint16_t elemIdx2;
	int16_t prefabIdx;
	float massThreshold;
	float probability;
};
struct Buffer
{
	uint32_t mSize;
	char* mGuardData;
	char* mData;

	explicit Buffer(int szie) {
		if (szie) {
			//this->mGuardData = new char[](szie + 8);
			this->mGuardData = (char*)calloc(szie + 8, sizeof(char));
		}

		if (szie && this->mGuardData) {
			this->mSize = szie;
			this->mData = this->mGuardData + 4;
			*(int*)&this->mGuardData[0       ] = -559038737; // 0xdeadbeef
			*(int*)&this->mGuardData[szie + 4] = -559038737; // 0xdeadbeef
		}
		else {
			this->mSize = 0;
			this->mGuardData = 0;
			this->mData = 0;
		}
	};
};

//----- global class -----
class BinaryBufferReader
{
	int   mOffset;
	char* mBufferData;
	int   mBufferSize;
public:
	explicit BinaryBufferReader(int msg_length, char* msg) {
		this->mOffset = 0;
		this->mBufferData = msg;
		this->mBufferSize = msg_length;
	}

	void PRINT_USAGE() {
		LOGGER_PRINT("BinaryBufferReader: data usage %d/%d.\n", this->mOffset, this->mBufferSize);
	}

	// msg is created by C# and can't be release by C++(SimDLL)
	~BinaryBufferReader() {
		//free(this->mBufferData);
		PRINT_USAGE();
		LOGGER_PRINT("BinaryBufferReader Wasted\n\n");
	}

	template <typename T>
	BinaryBufferReader& operator>>(T& _right) {
		_right = *(T*)&this->mBufferData[this->mOffset];
		this->mOffset += sizeof(T);
		if (this->mOffset > this->mBufferSize)
			LOGGER_PRINT("BinaryBufferReader %s Over Size!!! Acquire %d with data size %d.\n", __func__, this->mOffset, this->mBufferSize);
		return *this;
	}

	void ReadBytes(int len, void* dest) {
		memmove(dest, &this->mBufferData[this->mOffset], len);
		this->mOffset += len;
		if (this->mOffset > this->mBufferSize)
			LOGGER_PRINT("BinaryBufferReader %s Over Size!!! Acquire %d with data size %d.\n", __func__, this->mOffset, this->mBufferSize);
	}

	void SkipBytes(int num_bytes) {
		this->mOffset += num_bytes;
		if (this->mOffset > this->mBufferSize)
			LOGGER_PRINT("BinaryBufferReader %s Over Size!!! Acquire %d with data size %d.\n", __func__, this->mOffset, this->mBufferSize);
	}

	void ReadString(std::string* str) {
		int length;
		*this >> length;
		str->assign(this->GetPointer(), length);
		this->SkipBytes(length);
	}

	char* GetPointer() {
		return &this->mBufferData[this->mOffset];
	}
};

class BinaryBufferWriter
{
	Buffer*  mBuffer;
	uint64_t mOffset;
public:
	explicit BinaryBufferWriter(Buffer* buffer) {
		this->mBuffer = buffer;
		this->mOffset = 0;
	}

	~BinaryBufferWriter() {
		// Don't release this->mBuffer. It will be sent to main program for saving.
		// this->mBuffer will be released with SaveBuffer.
		PRINT_USAGE();
		LOGGER_PRINT("BinaryBufferWriter Wasted\n\n");
	}

	template <typename T>
	BinaryBufferWriter& operator<<(T _right) {
		*(T*)&this->mBuffer->mData[this->mOffset] = _right;
		this->mOffset += sizeof(T);
		if (this->mOffset > this->mBuffer->mSize)
			ASSERT_TEXT("BinaryBufferWriter Over Size!!!");
		return *this;
	}

	void WriteBytes(uint64_t len, void* src) {
		memmove(&this->mBuffer->mData[this->mOffset], src, len);
		this->mOffset += len;
		if (this->mOffset > this->mBuffer->mSize)
			ASSERT_TEXT("BinaryBufferWriter Over Size!!!");
	}

	void WriteString(std::string* str) {
		uint64_t length = str->size();
		*this << length;
		this->WriteBytes(length, str);
	}

	void PRINT_USAGE() {
		LOGGER_PRINT("BinaryBufferWriter: data usage %lld/%d.\n", this->mOffset, this->mBuffer->mSize);
	}
};

//----- global variable -----
extern std::vector<Element> gElements;
extern std::vector<ElementLiquidData> gElementLiquidData;
extern std::vector<ElementLightAbsorptionData> gElementLightAbsorptionData;
extern std::vector<PhysicsData> gElementPhysicsData;
extern std::vector<ElementRadiationData> gElementRadiationData;
extern std::vector<ElementPostProcessData> gElementPostProcessData;
extern std::vector<ElementStateData> gElementStateData;
extern std::vector<ElementPressureData> gElementPressureData;
extern std::vector<std::string> gElementNames;
extern std::unordered_map<uint32_t, uint16_t> gElementIndices;
extern std::vector<ElementPropertyTextureData> gElementPropertyTextureData;
extern std::vector<ElementTemperatureData> gElementTemperatureData;
extern std::unique_ptr<SimData> gSimData;
extern std::unique_ptr<Disease> gDisease;
extern std::vector<GasObliteration> gGasObliterations;
extern std::vector<LiquidConversion> gLiquidConversions;
extern std::vector<LiquidObliteration> gLiquidObliterations;
extern std::unique_ptr<Sim> gSim;
extern void* (*gGameMessageHandler)(int, const void*);
extern FrameSync gFrameSync;
extern GameDataUpdate gGameDataUpdate;
extern Buffer* SaveBuffer;

//----- global Function -----
extern uint16_t GetElementIndex(uint32_t hash);
extern void  FindCellbyDepth(std::vector<uint16_t> elemIdx, int width, uint16_t src_x, uint16_t src_y, short remain_depth,
	bool stop_at_solid, bool stop_at_liquid, bool stop_at_gas, std::unordered_map<int, short>* cellDic);
extern bool  cmpDepth(std::pair<int, short> a, std::pair<int, short> b);
extern float CalculateFinalTemperature(float mass1, float temp1, float mass2, float temp2);
extern void  AddMassAndUpdateTemperature(CellSOA* read_cells, CellSOA* write_cells, int cell_idx, float mass_delta, float temperature, uint8_t disease_idx, int disease_count);
extern bool  DisplaceGas   (SimData* simData, SimEvents* simEvents, CellSOA* cells, int cell_idx, uint16_t elem_idx);
extern bool  DisplaceLiquid(SimData* simData, SimEvents* simEvents, CellSOA* cells, int cell_idx, uint16_t elem_idx);
extern bool  DoStateTransition(SimData* simData, SimEvents* simEvents, int cell, const ElementTemperatureData* elem);
extern void  Evaporate(int cell_idx, SimData* simData, SimEvents* simEvents);
extern float DoGasPressureDisplacement(uint16_t elem_idx, const int cell, const int ncell, const int nncell, SimData* simData);
extern float DoLiquidPressureDisplacement(uint16_t elem_idx, const int cell, const int ncell, int nncell, SimData* simData);
extern void  PostProcessCell(SimData* simData, SimEvents* simEvents, int cell_idx);
extern float UpdatePressure(SimData* simData, SimEvents* simEvents, CellSOA* cells, CellSOA* updated_cells,
	int cell, uint16_t elem_idx, uint8_t elem_state, float elem_flow, int ncell);
extern void  UpdateLiquid(SimData* simData, SimEvents* simEvents, int cell_idx);
extern void  UpdateTemperature(SimData* simData, SimEvents* simEvents, const int cell, const int ncell);
#endif
