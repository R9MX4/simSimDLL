#include "ClassDisease.h"
#include "ClassBase.h"

#include <chrono>
#include <iostream>
#include <sstream> 

using namespace std;

// New_add
FILE*  gLogger, * gLogger2;
int    gLogLevel = 0;
double gGasDisplace = 1;

std::vector<Element> gElements;
std::vector<ElementLiquidData> gElementLiquidData;
std::vector<ElementLightAbsorptionData> gElementLightAbsorptionData;
std::vector<PhysicsData> gElementPhysicsData;
std::vector<ElementRadiationData> gElementRadiationData;
std::vector<ElementPostProcessData> gElementPostProcessData;
std::vector<ElementStateData> gElementStateData;
std::vector<ElementPressureData> gElementPressureData;
std::vector<std::string> gElementNames;
std::unordered_map<uint32_t, uint16_t> gElementIndices;
std::vector<ElementPropertyTextureData> gElementPropertyTextureData;
std::vector<ElementTemperatureData> gElementTemperatureData;
std::unique_ptr<SimData> gSimData;
std::unique_ptr<Disease> gDisease;
std::vector<GasObliteration> gGasObliterations;
std::vector<LiquidConversion> gLiquidConversions;
std::vector<LiquidObliteration> gLiquidObliterations;
std::unique_ptr<Sim> gSim;
void* (*gGameMessageHandler)(int, const void*) = NULL;
FrameSync gFrameSync;
GameDataUpdate gGameDataUpdate;
Buffer* SaveBuffer;

std::chrono::steady_clock::time_point  tp_init;
void INIT_TIMER()
{
	tp_init = std::chrono::steady_clock::now();
}

void PRINT_TIME(FILE* _Logger)
{
	std::chrono::steady_clock::time_point tp = std::chrono::steady_clock::now();
	std::ostringstream ostring;
	auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_init);
	auto micros = std::chrono::duration_cast<std::chrono::microseconds>(tp - tp_init);
	ostring << millis.count() << '.' << setw(3) << micros.count() % 1000 << '\t';
	fprintf(_Logger, ostring.str().c_str());
}

__inline float MAX_F(float a, float b)
{
	if (a > b) return a;
	return b;
}

__inline float MIN_F(float a, float b)
{
	if (a < b) return a;
	return b;
}

__inline float CLAMP_F(float val, float max, float min)
{
	if (val >= max) return max;
	if (val <= min) return min;
	return val;
}