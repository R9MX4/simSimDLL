#pragma once
#ifndef CLASS_DISEASE_H
#define CLASS_DISEASE_H

#include "global.h"
//----- define struct -----
struct DiseaseResult
{
	uint8_t idx;
	int count;
};
struct DiseaseInfo
{
public:
	struct RangeInfo
	{
		float minViable;
		float minGrowth;
		float maxGrowth;
		float maxViable;
	};

	struct ElemGrowthInfo
	{
		std::vector<float> underPopulationDeathRate;
		std::vector<float> populationHalfLife;
		std::vector<float> overPopulationHalfLife;
		std::vector<float> diffusionScale;
		std::vector<float> minCountPerKG;
		std::vector<float> maxCountPerKG;
		std::vector<int>   minDiffusionCount;
		std::vector<uint8_t> minDiffusionInfestationTickCount;
	};

	uint32_t hashID;
	float    strength;
	DiseaseInfo::RangeInfo		temperatureRange;
	DiseaseInfo::RangeInfo		temperatureHalfLives;
	DiseaseInfo::RangeInfo		pressureRange;
	DiseaseInfo::RangeInfo		pressureHalfLives;
	DiseaseInfo::ElemGrowthInfo	elemGrowthInfo;
	float radiationKillRate;
};

//----- define class -----
class Disease
{
	std::vector<std::string> diseaseNames;

public:
	std::vector<DiseaseInfo> diseases;

	explicit Disease            (BinaryBufferReader* reader);
	void	AddDiseaseToCell	(CellSOA* cells, int cell_idx, uint8_t new_disease_idx, int new_disease_count);
	float	GetDiffusionScale	(const SimData* simData, int cell_idx);
	uint8_t	GetDiseaseIndex		(uint32_t disease_hash);
	void	PostProcess			(SimData* simData, SimEvents* simEvents, int x_start, const int x_end, int y_start, int y_end);
	void	UpdateCells			(float dt, SimData* simData, SimEvents* simEvents, int cell, int ncell);
	DiseaseResult* CalculateFinalDiseaseCount(DiseaseResult* result, uint8_t src1_idx, int src1_count, uint8_t src2_idx, int src2_count);
};

#endif