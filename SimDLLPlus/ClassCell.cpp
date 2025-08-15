#include "ClassCell.h"

CellSOA::CellSOA(uint64_t num_cells) 
{
#ifdef __DEBUGED__
	this->elementIdx	= std::vector<uint16_t>	(num_cells, -1);
#else
	this->elementIdx    = std::vector<uint16_t> (num_cells, 255);
#endif
	this->temperature	= std::vector<float>	(num_cells);
	this->mass			= std::vector<float>	(num_cells);
	this->properties	= std::vector<uint8_t>	(num_cells);
	this->insulation	= std::vector<uint8_t>	(num_cells, -1);
	this->strengthInfo	= std::vector<uint8_t>	(num_cells, 0x84);
	this->diseaseIdx	= std::vector<uint8_t>	(num_cells, -1);
	this->diseaseCount	= std::vector<int>		(num_cells);
	this->diseaseInfestationTickCount	= std::vector<uint8_t>(num_cells);
	this->diseaseGrowthAccumulatedError	= std::vector<float>  (num_cells);
	this->radiation		= std::vector<float>	(num_cells);
}

void CellSOA::CopyFrom(const CellSOA* src)
{
#define CopyVector(_Key, _Type) \
	this->_Key = std::vector<_Type>(src->_Key.size()); \
	std::copy(src->_Key.begin(), src->_Key.end(), this->_Key.begin());

	//LOGGER_PRINT2("CellSOA %s: %lld->%lld\n", __func__, (uint64_t)src, (uint64_t)this);

	CopyVector(elementIdx,   uint16_t);
	CopyVector(temperature,  float);
	CopyVector(mass,         float);
	CopyVector(properties,   uint8_t);
	CopyVector(insulation,   uint8_t);
	CopyVector(strengthInfo, uint8_t);
	CopyVector(diseaseIdx,   uint8_t);
	CopyVector(diseaseCount, int);
	CopyVector(diseaseInfestationTickCount,   uint8_t);
	CopyVector(diseaseGrowthAccumulatedError, float);
	CopyVector(radiation,    float);
}

void CellSOA::SwapCells(int cell_idx1, int cell_idx2)
{
	std::swap(this->elementIdx	[cell_idx1], this->elementIdx	[cell_idx2]);
	std::swap(this->diseaseIdx	[cell_idx1], this->diseaseIdx	[cell_idx2]);
	std::swap(this->temperature	[cell_idx1], this->temperature	[cell_idx2]);
	std::swap(this->mass		[cell_idx1], this->mass			[cell_idx2]);
	std::swap(this->diseaseCount[cell_idx1], this->diseaseCount	[cell_idx2]);
	std::swap(this->diseaseInfestationTickCount  [cell_idx1], this->diseaseInfestationTickCount  [cell_idx2]);
	std::swap(this->diseaseGrowthAccumulatedError[cell_idx1], this->diseaseGrowthAccumulatedError[cell_idx2]);
	// Radiation does not swap with element
}

void CellSOA::ClearDisease(int cell_idx)
{
	this->diseaseIdx					[cell_idx] = -1;
	this->diseaseCount					[cell_idx] = 0;
	this->diseaseInfestationTickCount	[cell_idx] = 0;
	this->diseaseGrowthAccumulatedError	[cell_idx] = 0;
}

void CellSOA::ModifyDiseaseCount(int cell_idx, int delta)
{
	if (this->diseaseCount[cell_idx] + delta > 0) {
		this->diseaseCount[cell_idx] += delta;
	}
	else {
		this->diseaseIdx					[cell_idx] = -1;
		this->diseaseCount					[cell_idx] = 0;
		this->diseaseInfestationTickCount	[cell_idx] = 0;
		this->diseaseGrowthAccumulatedError	[cell_idx] = 0;
	}
}

void CellAccessor::ClearDisease()
{
	this->cells->diseaseIdx                   [this->cellIdx] = -1;
	this->cells->diseaseCount                 [this->cellIdx] = 0;
	this->cells->diseaseInfestationTickCount  [this->cellIdx] = 0;
	this->cells->diseaseGrowthAccumulatedError[this->cellIdx] = 0;
}

void CellAccessor::ModifyDiseaseCount(int delta)
{
	this->cells->diseaseCount[cellIdx] += delta;
	if (cells->diseaseCount[cellIdx] > 0)
		return;
	this->cells->diseaseIdx                   [this->cellIdx] = -1;
	this->cells->diseaseCount                 [this->cellIdx] = 0;
	this->cells->diseaseInfestationTickCount  [this->cellIdx] = 0;
	this->cells->diseaseGrowthAccumulatedError[this->cellIdx] = 0;
}

void CellAccessor::SetDisease(uint8_t disease_idx, int disease_count)
{
	this->cells->diseaseIdx                   [this->cellIdx] = disease_idx;
	this->cells->diseaseCount                 [this->cellIdx] = disease_count;
	this->cells->diseaseInfestationTickCount  [this->cellIdx] = 0;
	this->cells->diseaseGrowthAccumulatedError[this->cellIdx] = 0;
}

float    CellAccessor::temperature()  { return this->cells->temperature [this->cellIdx]; }
int      CellAccessor::diseaseCount() { return this->cells->diseaseCount[this->cellIdx]; }
uint8_t  CellAccessor::diseaseIdx()   { return this->cells->diseaseIdx  [this->cellIdx]; }
uint16_t CellAccessor::elementIdx()   { return this->cells->elementIdx  [this->cellIdx]; }
float    CellAccessor::mass()         { return this->cells->mass        [this->cellIdx]; }
void     CellAccessor::modifyMass(float delta){this->cells->mass[this->cellIdx] += delta; }