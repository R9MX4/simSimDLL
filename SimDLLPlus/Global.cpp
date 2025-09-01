#include "ClassDisease.h"
#include "ClassBase.h"

#include <chrono>
#include <iostream>
#include <sstream> 
#include <stdarg.h>

using namespace std;

// New_add
FILE* gLogger = 0, * gLogger2 = 0;
int   gLogLevel    = 0;
float gGasDisplace = 1;

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

std::chrono::steady_clock::time_point tp_init;
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

void PRINT_LOG(FILE* _Logger, const char* fmt, ...)
{
	char buffer[256] = { 0 };
	va_list args;
	va_start(args, fmt);
	vsnprintf(buffer, 255, fmt, args);
	va_end(args);

	std::chrono::steady_clock::time_point tp = std::chrono::steady_clock::now();
	std::ostringstream ostring;
	auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_init);
	auto micros = std::chrono::duration_cast<std::chrono::microseconds>(tp - tp_init);
	ostring << millis.count() << '.' << setw(3) << micros.count() % 1000 << '\t' << buffer;
	fprintf(_Logger, ostring.str().c_str());

	std::fflush(_Logger);
}

void PRINT_LOG_PARA(FILE* _Logger, const char* fmt, ...)
{
	int thrIdx = omp_get_thread_num();

	char buffer[256] = { 0 };
	va_list args;
	va_start(args, fmt);
	vsnprintf(buffer, 255, fmt, args);
	va_end(args);

	std::chrono::steady_clock::time_point tp = std::chrono::steady_clock::now();
	std::ostringstream ostring;
	auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tp_init);
	auto micros = std::chrono::duration_cast<std::chrono::microseconds>(tp - tp_init);
	ostring << millis.count() << '.' << setw(3) << micros.count() % 1000 << "\tThread-" << thrIdx << ": " << buffer;
	fprintf(_Logger, ostring.str().c_str());

	std::fflush(_Logger);
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