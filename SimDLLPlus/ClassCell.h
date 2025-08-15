#pragma once
#ifndef CLASS_CELL_H
#define CLASS_CELL_H

#include "global.h"
//----- define class -----
class CellSOA
{
public:
		std::vector<uint16_t> elementIdx;
		std::vector<float	> temperature;
		std::vector<float	> mass;
		std::vector<uint8_t	> properties; 
		//Sim.Cell.Properties: GasImpermeable=1; LiquidImpermeable=2; SolidImpermeable=4; 
		//Unbreakable=8; Transparent=16; Opaque=32; NotifyOnMelt=64; ConstructedTile=128;
		std::vector<uint8_t	> insulation;
		std::vector<uint8_t	> strengthInfo;
		// SimMessages.SetStrength: ConstructedTile: strengthMultiplier * 4, Nature Tile 132
		std::vector<uint8_t	> diseaseIdx;
		std::vector<int		> diseaseCount;
		std::vector<uint8_t	> diseaseInfestationTickCount;
		std::vector<float	> diseaseGrowthAccumulatedError;
		std::vector<float	> radiation;

		explicit CellSOA		(uint64_t num_cells);
		void CopyFrom			(const CellSOA * src);
		void SwapCells			(int cell_idx1, int cell_idx2);
		void ClearDisease		(int cell_idx);
		void ModifyDiseaseCount	(int cell_idx, int delta);
};

class CellAccessor
{
public:
	CellSOA* cells;
	int      cellIdx;

	CellAccessor(CellSOA* cells, int cellIdx)
	{
		this->cells = cells;
		this->cellIdx = cellIdx;
	}

	void     ClearDisease();
	void     ModifyDiseaseCount(int delta);
	void     SetDisease(uint8_t disease_idx, int disease_count);
	float    temperature();
	int      diseaseCount();
	uint8_t  diseaseIdx();
	uint16_t elementIdx();
	float    mass();
	void     modifyMass(float delta);
};
#endif
