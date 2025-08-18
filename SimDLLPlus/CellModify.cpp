#include "global.h"
#include "ClassCell.h"
#include "ClassDisease.h"
#include "ClassSim.h"
#include "GameData.h"

//----- (Local Struct) -----------------------------------------------
struct FloodRemovedInfo
{
    int16_t  remove_elem_idx;
    SimData* simData;
    float*   mass_to_remove;
    ConsumedMassInfo* removed_info;

    bool operator()(int cell);
};

//----- (Variable) ---------------------------------------------------
int displacement_offset = -1;

//----- (New add) ----------------------------------------------------
float CalculateFinalTemperature(float mass1, float temp1, float mass2, float temp2)
{
    if (mass1 + mass2 == 0) 
        return (temp1 + temp2) * 0.5f;

    float temp_min = MIN_F(temp1, temp2);
    float temp_max = MAX_F(temp1, temp2);
    float temp_new = ((mass1 * temp1) + (mass2 * temp2)) / (mass1 + mass2);
    return CLAMP_F(temp_new, temp_max, temp_min);
}

bool FloodRemovedInfo::operator()(int cell)
{
    if (this->simData->updatedCells->elementIdx[cell] == this->remove_elem_idx) {
        CellAccessor cellAccessor(this->simData->updatedCells.get(), cell);
        float massDelta = MIN_F(this->simData->updatedCells->mass[cell], *this->mass_to_remove);
        int   diseDelta = (int)(this->simData->updatedCells->diseaseCount[cell] * massDelta / this->simData->updatedCells->mass[cell]);
        this->simData->updatedCells->diseaseCount[cell] -= diseDelta;
        this->simData->updatedCells->mass[cell]         -= massDelta;
        *this->mass_to_remove                           -= massDelta;

        this->removed_info->temperature = CalculateFinalTemperature(
            massDelta, this->simData->updatedCells->temperature[cell], this->removed_info->mass, this->removed_info->temperature);
        this->removed_info->mass += massDelta;

        if (this->simData->updatedCells->mass[cell] <= FLT_MIN)
            Evaporate(cell, this->simData, this->simData->simEvents.get());
        else if (this->simData->updatedCells->diseaseCount[cell] <= 0)
            cellAccessor.ClearDisease();
        else if (this->simData->updatedCells->diseaseIdx[cell] != 0xFF) {
            DiseaseResult final_disease;
            gDisease->CalculateFinalDiseaseCount(&final_disease, this->removed_info->diseaseIdx,
                this->removed_info->diseaseCount, this->simData->updatedCells->diseaseIdx[cell], diseDelta);
            this->removed_info->diseaseIdx   = final_disease.idx;
            this->removed_info->diseaseCount = final_disease.count;
        }
    }
    return (*this->mass_to_remove) <= 0.0f;
}

int Flood(SimData* simData, uint16_t src_x, uint16_t src_y, int max_depth,
    bool stop_at_solid, bool stop_at_liquid, bool stop_at_gas, std::vector<int>* visited, std::queue<FloodFillInfo>* next, FloodRemovedInfo& fn)
{
    std::unordered_map<int, short> cellDic;
    for (; !next->empty(); next->pop()) {
        int cell = next->front().x + next->front().y * simData->width;
        cellDic[cell] = max_depth - next->front().depth;
    }
    FindCellbyDepth(simData->updatedCells->elementIdx, simData->width, src_x, src_y, max_depth, stop_at_solid, stop_at_liquid, stop_at_gas, &cellDic);
#if 1 // rewrite
#if 0 // unsorted
    for (const auto& pair : cellDic) {
#else // sorted
    std::vector<std::pair<int, short>> cellSorted(cellDic.begin(), cellDic.end());
    std::sort(cellSorted.begin(), cellSorted.end(), cmpDepth);

    for (const auto& pair : cellSorted) {
#endif
        visited->push_back(pair.first);
        if (fn(pair.first))
            return pair.first;
    }
    return -1;
#else // origin
    FloodFillInfo floodFillInfo = { .x = src_x, .y = src_y, .depth = 0 };
    next->push(floodFillInfo);
    if (!next->size())
        return -1;
    uint16_t locX, locY, depth, element;
    uint32_t cell;
    uint8_t  state;

    for (; !next->empty(); next->pop()) {
        locX = next->front().x;
        locY = next->front().y;
        depth = next->front().depth;
        cell = locX + locY * simData->width;
        if (depth >= max_depth || locX == 0 || locX >= simData->width || locY == 0 || locY >= simData->height) {
            if (!next->size())
                return -1;
            continue;
        }

        element = simData->updatedCells->elementIdx[cell];
        state = gElementPostProcessData[element].state & 3;
        if ((state != 3 || !stop_at_solid) &&
            (state != 2 || !stop_at_liquid) &&
            (state != 1 || !stop_at_gas)) {
            break;
        }
    label:;
    }

    int idx = 0;
    if (visited->size()) {
        for (; idx < visited->size(); idx++) {
            if (cell == (*visited)[idx])
                break;
        }
    }
    if (idx != visited->size()) {
        if (!next->size())
            return -1;
        goto label;
    }
    visited->push_back(cell);

    if (!fn(cell)) {
        FloodFillInfo floodFillInfo1 = { .x = locX,      .y = locY - 1u, .depth = depth + 1u };
        FloodFillInfo floodFillInfo2 = { .x = locX - 1u, .y = locY,      .depth = depth + 1u };
        FloodFillInfo floodFillInfo3 = { .x = locX,      .y = locY + 1u, .depth = depth + 1u };
        FloodFillInfo floodFillInfo4 = { .x = locX + 1u, .y = locY,      .depth = depth + 1u };
        next->push(floodFillInfo1);
        next->push(floodFillInfo2);
        next->push(floodFillInfo3);
        next->push(floodFillInfo4);
        if (!next->size())
            return -1;
        goto label;
    }

    return cell;
#endif
}

//----- (Decompiled Basic) -------------------------------------------
void AddMassAndUpdateTemperature(CellSOA* read_cells, CellSOA* write_cells, int cell_idx, float mass_delta, float temperature, uint8_t disease_idx, int disease_count)
{
    if ((read_cells->mass[cell_idx] + mass_delta) <= 0.0) {
        write_cells->temperature[cell_idx] = 0.0;
        CellAccessor write_cell(write_cells, cell_idx);
        write_cell.ClearDisease();
    }
    else {
        write_cells->temperature[cell_idx] = CalculateFinalTemperature(read_cells->mass[cell_idx], read_cells->temperature[cell_idx], mass_delta, temperature);
        write_cells->mass       [cell_idx] = read_cells->mass[cell_idx] + mass_delta;
        gDisease->AddDiseaseToCell(write_cells, cell_idx, disease_idx, disease_count);
    }

    ASSERT_TEMP(write_cells->mass[cell_idx], write_cells->temperature[cell_idx]);
}

void DoDisplacement(SimData* simData, SimEvents* simEvents, CellSOA* cells, int src_cell_idx, int dest_cell_idx)
{
    AddMassAndUpdateTemperature(cells, cells, dest_cell_idx,
        cells->mass      [src_cell_idx], cells->temperature [src_cell_idx],
        cells->diseaseIdx[src_cell_idx], cells->diseaseCount[src_cell_idx]);
    cells->elementIdx[dest_cell_idx] = cells->elementIdx[src_cell_idx];

    ASSERT_TEMP(cells->mass[dest_cell_idx], cells->temperature[dest_cell_idx]);
    CellAccessor src_cell(cells, src_cell_idx);
    simData->ClearCell(&src_cell);
    simEvents->ChangeSubstance(simData, src_cell_idx);
    simEvents->ChangeSubstance(simData, dest_cell_idx);
    return;
}

bool DisplaceGas(SimData* simData, SimEvents* simEvents, CellSOA* cells, int cell_idx, uint16_t elem_idx)
{
    if (cells->mass[cell_idx] <= 0.0)
        return false;
    if ((gElementPostProcessData[cells->elementIdx[cell_idx]].state & 3) != 1)
        return false;

    int width              = simData->width;
    uint16_t tickCount     = simData->tickCount;
    int candidate_cells[4] = { cell_idx + width, cell_idx - 1, cell_idx + 1, cell_idx - width };
    for (int direct = 0; direct < 4; direct++) {
        int dest_cell_idx = candidate_cells[(direct + tickCount) & 3];
        uint16_t elem_d = cells->elementIdx[dest_cell_idx];
        if ((elem_d == elem_idx || elem_d == simData->vacuumElementIdx) && ((cells->properties[dest_cell_idx] & 1) == 0)) {
            DoDisplacement(simData, simEvents, cells, cell_idx, dest_cell_idx);
            return true;
        }
    }

    int neighbor_cells[2] = { cell_idx - 1, cell_idx + 1 };
    int diagonal_candidates[2] = { cell_idx + width - 1, cell_idx + width + 1 };
    for (int direct = 0; direct < 2; direct++) {
        int dest_cell_idx = diagonal_candidates[(direct + tickCount) & 1];

        if (cells->elementIdx[dest_cell_idx] == elem_idx) {
            int neighbor_cell_idx = neighbor_cells[(direct + tickCount) & 1];
            uint16_t elem_n = simData->updatedCells->elementIdx[neighbor_cell_idx];

            if (((gElementLiquidData[elem_n].state & 3) != 3) && ((cells->properties[dest_cell_idx] & 1) == 0)) {
                DoDisplacement(simData, simEvents, cells, cell_idx, dest_cell_idx);
                return true;
            }
        }
    }
    return false;
}

bool DisplaceLiquidSimple(SimData* simData, SimEvents* simEvents, int cell_idx, uint16_t elem_idx, int num_candidate_cells, int* candidate_cells)
{
    if (num_candidate_cells > 4)
        ASSERT_TEXT("Assert failed: num_candidate_cells <= kMaxNeighbors");
    if (simData->updatedCells->mass[cell_idx] < 0.01 || num_candidate_cells <= 0)
        return false;

    int available_cells[4];
    int available_count = 0;
    for (int i = 0; i < num_candidate_cells; i++) {
        int      cell     = candidate_cells[i];
        uint16_t elemCand = simData->updatedCells->elementIdx[cell];
        if (elemCand == elem_idx || elemCand == simData->vacuumElementIdx)
            if (!(simData->updatedCells->properties[cell] & 2))
                available_cells[available_count++] = cell;
    }
    if (available_count == 0) return false;

    //Round Up
    int disease_count = (simData->updatedCells->diseaseCount[cell_idx] + available_count - 1) / available_count;
    for (int i = 0; i < available_count; i++) {
        int cellDest = available_cells[i];
        AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), cellDest,
            simData->updatedCells->mass[cell_idx] / available_count,
            simData->updatedCells->temperature[cell_idx],
            simData->updatedCells->diseaseIdx[cell_idx], disease_count);
        simData->updatedCells->elementIdx[cellDest] = elem_idx;
        int cellInner = CELL_S2G(cellDest, simData->width);
        if (CELL_AVAILABLE(cellInner, simData)) {
            SubstanceChangeInfo info{ cellInner, (uint16_t)-1, (uint16_t)-1 };
            simEvents->substanceChangeInfo.push_back(info);
        }
        simData->timers[cellDest].stableCellTicks |= 0x1F;
    }
    CellAccessor src_cell(simData->updatedCells.get(), cell_idx);
    simData->ClearCell(&src_cell);
    simEvents->ChangeSubstance(simData, cell_idx);
    return true;
}

bool DisplaceLiquid(SimData* simData, SimEvents* simEvents, CellSOA* cells, int cell_idx, uint16_t elem_idx)
{
    if (cells->mass[cell_idx] <= 0.0)
        return false;
    if ((gElementPostProcessData[cells->elementIdx[cell_idx]].state & 3) != 2)
        return false;

    int candidate_cells[4] = { cell_idx + 1, cell_idx - 1, cell_idx + simData->width, cell_idx - simData->width };
    if (DisplaceLiquidSimple(simData, simEvents, cell_idx, elem_idx, 4, candidate_cells))
        return true;

    for (int direct = 0; direct < 4; direct++) {
        int cell = candidate_cells[(direct + simData->tickCount) & 3];
        if (GET_STATE(cells, cell) == 1 && DisplaceGas(simData, simEvents, cells, cell, cells->elementIdx[cell]))
            return DisplaceLiquidSimple(simData, simEvents, cell_idx, elem_idx, 4, candidate_cells);
    }
    return false;
}

bool DoStateTransition(SimData* simData, SimEvents* simEvents, int cell, const ElementTemperatureData* elem)
{
	if (simData->updatedCells->mass[cell] <= 0.0) return false;

	//Low Temp Transition
#ifdef __SIMDLL_PLUS__ // Same as DoLoadTimeStateTransition
    if (simData->updatedCells->mass[cell] > elem->compressedLiquifyMass ||
        simData->updatedCells->temperature[cell] < (elem->lowTemp - 3.0f) && elem->lowTempTransitionIdx != 0xFFFF) {
#else
	if (simData->updatedCells->temperature[cell] < (elem->lowTemp - 3.0f) && elem->lowTempTransitionIdx != 0xFFFF) {
#endif
		//Update temperature
		//simData->updatedCells->temperature[cell] = MAX_F(simData->updatedCells->temperature[cell] + 1.5f, 0);
		simData->updatedCells->temperature[cell] = CLAMP_F(simData->updatedCells->temperature[cell] + 1.5f, SIM_MAX_TEMPERATURE, 1);
		//Sub Transition Target
		uint16_t subTargetIdx = elem->lowTempTransitionOreIdx;
		if (subTargetIdx != 0xFFFF) {
			float subTargetMass = elem->lowTempTransitionOreMassConversion * simData->updatedCells->mass[cell];
			if (subTargetMass > 0.001) {
				int disease_count = (int)(simData->updatedCells->diseaseCount[cell] * elem->lowTempTransitionOreMassConversion);
				if (simEvents->SpawnOre(simData, cell, subTargetIdx, subTargetMass, 
					simData->updatedCells->temperature[cell], simData->updatedCells->diseaseIdx[cell], 
					disease_count, true))
                {
					simData->updatedCells->mass[cell] -= subTargetMass;
					CellAccessor sim_cell(simData->updatedCells.get(), cell);
					sim_cell.ModifyDiseaseCount(-disease_count);
				}
			}
		}
		//Transition Target
		uint16_t mainTargetIdx = elem->lowTempTransitionIdx;
		if ((gElementTemperatureData[mainTargetIdx].state & 3) == 3) {
			//Target is Solid
			//Small mass, spwan debris
			if ((simData->updatedCells->mass[cell] / gElementTemperatureData[mainTargetIdx].defaultMass) <= 0.8
				&& simEvents->SpawnOre(simData, cell, mainTargetIdx, 
					simData->updatedCells->mass[cell], 
					simData->updatedCells->temperature[cell], 
					simData->updatedCells->diseaseIdx[cell], 
					simData->updatedCells->diseaseCount[cell], 
					simData->debugProperties.isDebugEditing)) 
			{
				CellAccessor sim_cell(simData->updatedCells.get(), cell);
				simData->ClearCell(&sim_cell);
				simEvents->ChangeSubstance(simData, cell);
				return true;
			}
			//Large mass, spwan cell
			simData->updatedCells->elementIdx[cell] = mainTargetIdx;
			simEvents->ChangeSubstance(simData, cell);
			return true;
		}
		else if (!simData->headless && (simData->cells->properties[cell - simData->width] & 2) == 0) {
			//Target is Liquid. Bottom cell is Liquid-permeable
			//Remove mass by spawning falling liquid
			int innerCell = CELL_S2G(cell, simData->width);
			if (CELL_AVAILABLE(innerCell, simData) && CELL_VISIABLE(innerCell, simData)
#ifdef __SIMDLL_PLUS__
                && simData->updatedCells->mass[cell] > elem->compressedLiquifyMass
                || (   simData->updatedCells->temperature[cell] >= gElements[mainTargetIdx].lowTemp  - 3.0f
                    && simData->updatedCells->temperature[cell] <= gElements[mainTargetIdx].highTemp + 3.0f))
#else
                && simData->updatedCells->temperature[cell] >= gElements[mainTargetIdx].lowTemp  - 3.0f
                && simData->updatedCells->temperature[cell] <= gElements[mainTargetIdx].highTemp + 3.0f)
#endif
			{
				SpawnFallingLiquidInfo spawnFallingLiquidInfo = {
					.cellIdx      = innerCell, 
					.elementIdx   = mainTargetIdx,
					.diseaseIdx	  = simData->updatedCells->diseaseIdx  [cell],
					.mass         = simData->updatedCells->mass        [cell],
					.temperature  = simData->updatedCells->temperature [cell],
					.diseaseCount = simData->updatedCells->diseaseCount[cell]
				};
				simEvents->spawnLiquidInfo.push_back(spawnFallingLiquidInfo);
				//Set original cell to vacuum 
				CellAccessor sim_cell(simData->updatedCells.get(), cell);
				simData->ClearCell(&sim_cell);
				simEvents->ChangeSubstance(simData, cell);
				return true;
			}
		}
		//Update element to Transition Target
		simData->updatedCells->elementIdx[cell] = mainTargetIdx;
		simEvents->ChangeSubstance(simData, cell);
		return true;
	}

	//High Temp Transition
#ifdef __SIMDLL_PLUS__
    if (simData->updatedCells->mass[cell] < elem->uncompressedgasifyMass &&
        simData->updatedCells->temperature[cell] > (elem->highTemp + 3.0f) && elem->highTempTransitionIdx != 0xFFFF) {
#else
	if (simData->updatedCells->temperature[cell] > (elem->highTemp + 3.0f) && elem->highTempTransitionIdx != 0xFFFF) {
#endif
		//Update temperature
		//simData->updatedCells->temperature[cell] = MAX_F(simData->updatedCells->temperature[cell] - 1.5f, 0);
		simData->updatedCells->temperature[cell] = CLAMP_F(simData->updatedCells->temperature[cell] - 1.5f, SIM_MAX_TEMPERATURE, 1);
		//Update element to Transition Target
		simData->updatedCells->elementIdx [cell] = elem->highTempTransitionIdx;

		//Sub Transition Target
		uint16_t subTargetIdx = elem->highTempTransitionOreIdx;
		if (subTargetIdx != 0xFFFF) {
			float subTargetMass = elem->highTempTransitionOreMassConversion * simData->updatedCells->mass[cell];
			if (subTargetMass > 0.001) {
				int disease_count = (int)(simData->updatedCells->diseaseCount[cell] * elem->highTempTransitionOreMassConversion);
				if (simEvents->SpawnOre(simData, cell, subTargetIdx, subTargetMass,
					simData->updatedCells->temperature[cell], simData->updatedCells->diseaseIdx[cell],
					disease_count, true)) 
                {
					simData->updatedCells->mass[cell] -= subTargetMass;
					CellAccessor sim_cell(simData->updatedCells.get(), cell);
					sim_cell.ModifyDiseaseCount(-disease_count);
				}
			}
		}
		//Property: NotifyOnMelt
		if ((simData->updatedCells->properties[cell] & 0x40) != 0) {
			int innerCell = CELL_S2G(cell, simData->width);
			if (CELL_AVAILABLE(innerCell, simData)) {
				struct CellMeltedInfo cellMeltedInfo = { .gameCell = (uint32_t)innerCell };
				simEvents->cellMeltedInfo.push_back(cellMeltedInfo);
			}
			simData->timers[cell].stableCellTicks |= 0x1F; // *(unsigned char*)&simData.timers[cell] |= 0x1Fu;
		}
		simEvents->ChangeSubstance(simData, cell);
		return true;
	}

	return false;
}

void Evaporate(int cell_idx, SimData* simData, SimEvents* simEvents)
{
    CellAccessor sim_cell(simData->updatedCells.get(), cell_idx);
    simData->ClearCell(&sim_cell);
    simEvents->ChangeSubstance(simData, cell_idx);
}

//----- (Decompiled) -------------------------------------------------
bool DisplaceLiquidDirectional(SimData* simData, SimEvents* simEvents, int src_cell_idx, const int dest_cell_idx)
{
    CellAccessor src_cell(simData->updatedCells.get(), src_cell_idx);
    uint16_t elem_d = simData->updatedCells->elementIdx[dest_cell_idx];
    uint8_t  state  = gElementLiquidData[elem_d].state & 3;

    if ((simData->updatedCells->mass[src_cell_idx] > 0.0) && ((simData->updatedCells->properties[dest_cell_idx] & 2) == 0)) {
        if ((state == 2) && (simData->updatedCells->elementIdx[src_cell_idx] == elem_d)) {
            AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), dest_cell_idx,
                simData->updatedCells->mass[src_cell_idx], simData->updatedCells->temperature[src_cell_idx],
                simData->updatedCells->diseaseIdx[src_cell_idx], simData->updatedCells->diseaseCount[src_cell_idx]);
            simData->ClearCell(&src_cell);
            simEvents->ChangeSubstance(simData, src_cell_idx);
            return true;
        }
        else if ((state <= 1) && DisplaceGas(simData, simEvents, simData->updatedCells.get(), dest_cell_idx, elem_d)) {
            if (simData->updatedCells->mass[dest_cell_idx] != 0.0)
                printf("dest_cell.mass() == 0.0f. File:%s Func:%s Line:%d\n", __FILE__, __func__, __LINE__);
            simData->updatedCells->SwapCells(src_cell_idx, dest_cell_idx);
            simEvents->ChangeSubstance(simData, src_cell_idx);
            simEvents->ChangeSubstance(simData, dest_cell_idx);
            return true;
        }
    }
    return false;
}

bool DoSublimation(SimData* simData, SimEvents* simEvents, int cell_idx, const ElementPostProcessData* elem)
{
    if (simData->updatedCells->mass[cell_idx] >= 1.8)
        return false;

    bool result = false;
    int neighbours[4] = { cell_idx + simData->width , cell_idx + 1 ,cell_idx - 1,cell_idx - simData->width };
    for (int celln : neighbours) {
        ElementPostProcessData* data = &gElementPostProcessData[simData->updatedCells->elementIdx[celln]];
        if ((data->state & 3) != 3) // Not Solid
            continue;
        if (data->sublimateIndex == 0xFFFF) // Not sublimatable
            continue;
        if (RAND_FLOAT_01(simData->randomSeed) > data->sublimateProbability) // Probability
            continue;

        float mass    = simData->updatedCells->mass[celln];
        float massSrc = MIN_F(data->sublimateRate * 0.2f, mass);
        // Fully sublimate. Condition: mass < 2 * mass_sublimate 
        if (mass - massSrc < massSrc) {
            simData->updatedCells->elementIdx[celln] = data->sublimateIndex;
            simData->updatedCells->mass      [celln] = mass * data->sublimateEfficiency;
            simEvents->ChangeSubstance(simData, celln);
            simEvents->SpawnFX(simData, celln, data->sublimateFX, 0);
            result = true;
            continue;
        }

        int disease_count = (int)(massSrc / simData->updatedCells->mass[celln] * simData->updatedCells->diseaseCount[celln]);
        // Same element, merge mass
        if (simData->updatedCells->elementIdx[cell_idx] == data->sublimateIndex) {
            float massDest = fmin(massSrc * data->sublimateEfficiency, 1.8f - simData->updatedCells->mass[cell_idx]);
            massSrc        = massDest / data->sublimateEfficiency;
            disease_count  = (int)(massSrc / simData->updatedCells->mass[celln] * simData->updatedCells->diseaseCount[celln]);
            AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), cell_idx, massDest,
                simData->updatedCells->temperature[celln], simData->updatedCells->diseaseIdx[celln], disease_count);
        }
        else if (!(elem->state & 3) 
            || DisplaceGas(simData, simEvents, simData->updatedCells.get(), cell_idx, simData->updatedCells->elementIdx[cell_idx])) 
        { // Vacuum or Remove Gas
            simData->updatedCells->elementIdx  [cell_idx] = data->sublimateIndex;
            simData->updatedCells->mass        [cell_idx]+= massSrc * data->sublimateEfficiency;
            simData->updatedCells->temperature [cell_idx] = simData->updatedCells->temperature[celln];
            simData->updatedCells->diseaseIdx  [cell_idx] = simData->updatedCells->diseaseIdx [celln];
            simData->updatedCells->diseaseCount[cell_idx] = disease_count;
            simEvents->ChangeSubstance(simData, cell_idx);
        }
#ifdef __DEBUGED__
        else continue;
#else
        else if (!result) continue;
#endif

        result = true;
        simData->updatedCells->mass[celln] -= massSrc;
        CellAccessor nsim_cell(simData->updatedCells.get(), celln);
        nsim_cell.ModifyDiseaseCount(-disease_count);
        simData->accumulatedFlow[celln] += massSrc;
        float rotation = 0;
        if      (celln - cell_idx == -1            ) rotation = -90;
        else if (celln - cell_idx == simData->width) rotation = 180;
        else if (celln - cell_idx == 1             ) rotation =  90;
        simEvents->SpawnFX(simData, celln, data->sublimateFX, rotation);
        if (simData->updatedCells->mass[celln] <= FLT_MIN)
            Evaporate(celln, simData, simEvents);
        continue;
    }
    return result;
}

bool DoDensityDisplacement(SimData* simData, SimEvents* simEvents, const int cell_idx, const ElementPostProcessData* elem, float required_probability)
{
    int cellDown = cell_idx - simData->width;
    Element* elemDown = &gElements[simData->updatedCells->elementIdx[cellDown]];
#ifndef __DEBUGED__ //123 Self Write
    if (elem->state != elemDown->state) return false;

    // Different element
    if (simData->updatedCells->elementIdx[cellDown] != simData->updatedCells->elementIdx[cell_idx]) {
        if (RAND_FLOAT_01(simData->randomSeed) <= required_probability)
            return false;
        if ((elem->state & 3) == 1) { // Gas
            if (elem->molarMass < elemDown->molarMass) required_probability = 1 - elem->molarMass / elemDown->molarMass * 0.5f;
            else                                       required_probability =     elemDown->molarMass / elem->molarMass * 0.5f;
            if (RAND_FLOAT_01(simData->randomSeed) <= required_probability)
                return false;
        }
        else if (elem->molarMass < elemDown->molarMass)
            false;
        // Do Random density replacement
        simData->updatedCells->SwapCells(cell_idx, cellDown);
        simEvents->ChangeSubstance(simData, cell_idx);
        simEvents->ChangeSubstance(simData, cellDown);
        return true;
    }
    else {
        // Same element, do thermal convection (Swap)
        float tempDown = simData->updatedCells->temperature[cellDown];
        float tempCent = simData->updatedCells->temperature[cell_idx];
        if (tempDown > tempCent)
            simData->updatedCells->SwapCells(cell_idx, cellDown);
    }
    return false;
#else
    if (elem->state == elemDown->state) {
        if (RAND_FLOAT_01(simData->randomSeed) > required_probability) {
            // Do density replacement
            if (elem->molarMass > elemDown->molarMass) {
                simData->updatedCells->SwapCells(cell_idx, cellDown);
                simEvents->ChangeSubstance(simData, cell_idx);
                simEvents->ChangeSubstance(simData, cellDown);
                return true;
            }
            // Same element, execute thermal convection
            if (simData->updatedCells->elementIdx[cell_idx] == simData->updatedCells->elementIdx[cellDown]) {
                float tempDown = simData->updatedCells->temperature[cellDown];
                float tempCent = simData->updatedCells->temperature[cell_idx];
                if (tempDown > tempCent) {
                    // Gas. Swap
                    if ((elem->state & 3) != 2) {
                        simData->updatedCells->SwapCells(cell_idx, cellDown);
                        return true;
                    }
                    // Liquid. Average temperature
                    float tempFin = CalculateFinalTemperature(simData->updatedCells->mass[cellDown], tempDown, simData->updatedCells->mass[cell_idx], tempCent);
                    simData->updatedCells->temperature[cell_idx] = tempFin;
                    simData->updatedCells->temperature[cellDown] = tempFin;
                    return true;
                }
            }
        }
    }
    return false;
#endif
}

bool DoPartialMelt(SimData* simData, SimEvents* simEvents, int cell_idx, int ncell_idx)
{
    // Situation: [Heat source: Gas] and [Target cell: Solid] and [Trans Target: Gas]
    uint16_t elemIdxC = simData->updatedCells->elementIdx[cell_idx];  // Heat source
    uint16_t elemIdxN = simData->updatedCells->elementIdx[ncell_idx]; // Melted cell
    Element* elemC    = &gElements[elemIdxC];
    Element* elemN    = &gElements[elemIdxN];
    float tempC       = simData->updatedCells->temperature[cell_idx];
    float tempN       = simData->updatedCells->temperature[ncell_idx];
    float tempTransN  = elemN->highTemp + 3.0f;
    uint16_t highTempTransitionIdx = elemN->highTempTransitionIdx;

    if ((elemC->state & 3) >  1)                             return false; // Heat source is not gas
    if ((elemN->state & 3) != 3)                             return false; // Melted cell is not solid
    if (simData->updatedCells->mass[ncell_idx] <= 5.0f)      return false; // Melted cell doesn't have enough mass
    if (simData->updatedCells->properties[ncell_idx] & 0x48) return false; // NotifyOnMelt or Unbreakable
    if (tempC <= tempTransN)                                 return false; // Heat source is not hot enough
    if (tempN >= elemN->highTemp - 3.0)                      return false; // Melted cell is too hot
    if (highTempTransitionIdx == 0xFFFF)                     return false; // Melted cell is not meltable

    float heatTrans = ((tempTransN - tempN) * gElements[highTempTransitionIdx].specificHeatCapacity) * 5.0f;
    float SHC_C     = simData->updatedCells->mass[cell_idx] * elemN->specificHeatCapacity;
    if ((tempC - elemN->lowTemp - 6.0f) * SHC_C < heatTrans) return false; // Heat source can't provide enough heat
    if ((gElements[highTempTransitionIdx].state & 3) != 2)   return false; // Melt target is not liquid

    // Melt 5kg, Spawn Falling Liquid
    int disease_count = (int)(5.0f / simData->updatedCells->mass[ncell_idx] * simData->updatedCells->diseaseCount[ncell_idx]);
    if (simEvents->SpawnFallingLiquid(simData, cell_idx, elemN->highTempTransitionIdx, 5.0f, tempTransN, 
        simData->updatedCells->diseaseIdx[ncell_idx], disease_count, simData->debugProperties.isDebugEditing)) 
    {
        simData->updatedCells->temperature [cell_idx]  = (SHC_C * tempC - heatTrans) / SHC_C;
        simData->updatedCells->mass        [ncell_idx] -= 5.0f;
        simData->updatedCells->diseaseCount[ncell_idx] -= disease_count;
        if (simData->updatedCells->mass[ncell_idx] <= FLT_MIN)
            Evaporate(ncell_idx, simData, simEvents);
        return true;
    }
    // Melt 5kg, Remove gas in Heat source cell
#ifdef __DEBUGED__
    simData->updatedCells->temperature[ncell_idx] = (SHC_C * tempC - heatTrans) / SHC_C;
    if (DisplaceGas(simData, simEvents, simData->updatedCells.get(), cell_idx, elemIdxC)) {
        simData->updatedCells->mass        [cell_idx]  = 5.0f;
        simData->updatedCells->temperature [cell_idx]  = tempTransN;
        simData->updatedCells->elementIdx  [cell_idx]  = elemN->highTempTransitionIdx;
        simData->updatedCells->diseaseIdx  [cell_idx]  = simData->updatedCells->diseaseIdx[ncell_idx];
        simData->updatedCells->diseaseCount[cell_idx]  = disease_count;
        simData->updatedCells->mass        [ncell_idx] -= 5.0f;
        simData->updatedCells->diseaseCount[ncell_idx] -= disease_count;
        if (simData->updatedCells->mass[ncell_idx] <= FLT_MIN)
            Evaporate(ncell_idx, simData, simEvents);
        simEvents->ChangeSubstance(simData, cell_idx);
        return true;
    }
    // Fail to DisplaceGas, recover temperature
    simData->updatedCells->temperature[cell_idx] = tempC;
#else
    if (DisplaceGas(simData, simEvents, simData->updatedCells.get(), cell_idx, elemIdxC)) {
        simData->updatedCells->temperature[cell_idx]  = (SHC_C * tempC - heatTrans) / SHC_C;
        simData->updatedCells->mass       [cell_idx]  = 5.0f;
        simData->updatedCells->temperature[cell_idx]  = tempTransN;
        simData->updatedCells->elementIdx [cell_idx]  = elemN->highTempTransitionIdx;
        simData->updatedCells->mass       [ncell_idx] -= 5.0f;
        if (simData->updatedCells->mass[ncell_idx] <= FLT_MIN)
            Evaporate(ncell_idx, simData, simEvents);
        simEvents->ChangeSubstance(simData, cell_idx);
        return true;
    }
#endif
    return false;
}

void HeadlessUnstableFallSimpleSwap(SimData* simData, SimEvents* simEvents, CellAccessor* source_cells)
{
    for (int cellDown = source_cells->cellIdx - simData->width; cellDown < 0; cellDown -= simData->width) {
        int      cellTop  = simData->width + cellDown;
        uint16_t elemDown = simData->updatedCells->elementIdx[cellDown];

        // Bottom cell is solid, stop drop
        if ((gElements[elemDown].state & 3) == 3) {
            simEvents->ChangeSubstance(simData, cellTop);
            return;
        }
        // Keep droping. Swap with lower cells
        simData->updatedCells->SwapCells(cellDown, cellTop);
        int cellInner = CELL_S2G(cellTop, simData->width);
        if (CELL_AVAILABLE(cellInner, simData)) {
            SubstanceChangeInfo info{ cellInner, (uint16_t)-1, (uint16_t)-1 };
            simEvents->substanceChangeInfo.push_back(info);
        }
        simData->timers[cellTop].stableCellTicks |= 0x1F;
    }
    return;
}

void DoUnstableCheckBasic(SimData* simData, SimEvents* simEvents, CellAccessor* sim_cell)
{
    int      cellDown = sim_cell->cellIdx - simData->width;
    uint16_t elemDown = simData->updatedCells->elementIdx[cellDown];

    // Bottom cell is solid/Void/SolidImpermeable
    if ((gElements[elemDown].state & 3) == 3
        || gElements[elemDown].id == simData->voidElementIdx
        || simData->updatedCells->properties[cellDown] & 4)
    {
        simData->timers[sim_cell->cellIdx].stableCellTicks |= 0x1F;
        return;
    }
    // Hold 3~5 ticks
    if (simData->GetStableTicksRemaining(sim_cell->cellIdx))
        return;
    // Drop
    if (simData->headless) {
        HeadlessUnstableFallSimpleSwap(simData, simEvents, sim_cell);
    }
    else {
        UnstableCellInfo info{
            .cellIdx      = CELL_S2G(sim_cell->cellIdx, simData->width),
            .elemIdx      = sim_cell->elementIdx(),
            .fallingInfo  = 0,
            .diseaseIdx   = sim_cell->diseaseIdx(),
            .mass         = sim_cell->mass(),
            .temperature  = sim_cell->temperature(),
            .diseaseCount = sim_cell->diseaseCount(),
        };
        simEvents->unstableCellInfo.push_back(info);
        Evaporate(sim_cell->cellIdx, simData, simEvents);
    }
}

void DoUnstableCheckWithDiagonals(SimData* simData, SimEvents* simEvents, CellAccessor* sim_cell)
{
    UnstableOffsetInfo infos[3] = {
        //offset            , massPercent, srcMinMassScale, checkAbove
         -simData->width    , 1.0f       , 0.0f           , false,
         -simData->width - 1, 0.5f       , 1.0f           , true,
         -simData->width + 1, 0.5f       , 1.0f           , true,
    };
    uint16_t elemIdx    = sim_cell->cells->elementIdx[sim_cell->cellIdx];
    uint8_t  remainTick = 255;
    bool     resetFlag  = true;
    for (int i = 0; i < 3; i++) {
        UnstableOffsetInfo& info = infos[i];
        if (sim_cell->cells->mass[sim_cell->cellIdx] <= info.srcMinMassScale * gElementPostProcessData[elemIdx].minHorizontalFlow * info.massPercent)
            continue;
        int cellOfs  = sim_cell->cellIdx + info.offset;
        uint8_t state = gElementPostProcessData[simData->updatedCells->elementIdx[cellOfs]].state & 3;
        // Offset cell is Solid/SolidImpermeable/Void
        if (state == 3) continue;
        if (simData->updatedCells->elementIdx[cellOfs] == simData->voidElementIdx) continue;
        if (simData->updatedCells->properties[cellOfs] & 4) continue;
        if (info.checkAbove) {
            int cellAbv = cellOfs + simData->width;
            if ((gElementPostProcessData[simData->updatedCells->elementIdx[cellAbv]].state & 3) == 3) continue;
            if (simData->updatedCells->properties[cellAbv] & 4) continue;
        }
        resetFlag = false;

        if (!i) remainTick = simData->GetStableTicksRemaining(sim_cell->cellIdx);
        if (remainTick) return;

        simData->timers[sim_cell->cellIdx].stableCellTicks |= 0x1F;
        if (simData->headless) {
            CellAccessor target_cell{ simData->updatedCells.get(),cellOfs };
            if (state == 1)
                DisplaceGas(simData, simEvents, simData->updatedCells.get(), cellOfs, simData->updatedCells->elementIdx[cellOfs]);
            else if (state == 2)
                DisplaceLiquid(simData, simEvents, simData->updatedCells.get(), cellOfs, simData->updatedCells->elementIdx[cellOfs]);

            simData->updatedCells->elementIdx  [cellOfs] = sim_cell->elementIdx();
            simData->updatedCells->mass        [cellOfs] = sim_cell->mass() * info.massPercent;
            simData->updatedCells->temperature [cellOfs] = sim_cell->temperature();
            simData->updatedCells->diseaseIdx  [cellOfs] = sim_cell->diseaseIdx();
            simData->updatedCells->diseaseCount[cellOfs] = int(sim_cell->diseaseCount() * info.massPercent);
            simEvents->ChangeSubstance(simData, cellOfs);
            sim_cell->modifyMass(fabs(sim_cell->mass() * info.massPercent));
            sim_cell->ModifyDiseaseCount((int)(-sim_cell->CellAccessor::diseaseCount() * info.massPercent));
            simEvents->ChangeSubstance(simData, sim_cell->cellIdx);
            if (sim_cell->mass() <= FLT_MIN)
                Evaporate(sim_cell->cellIdx, simData, simEvents);
            HeadlessUnstableFallSimpleSwap(simData, simEvents, &target_cell);
        }
        else {
            UnstableCellInfo infoCell{
                .cellIdx      = CELL_S2G(cellOfs, simData->width),
                .elemIdx      = sim_cell->elementIdx(),
                .fallingInfo  = 0,
                .diseaseIdx   = sim_cell->diseaseIdx(),
                .mass         = sim_cell->cells->mass[sim_cell->cellIdx] * info.massPercent,
                .temperature  = sim_cell->temperature(),
                .diseaseCount = int(sim_cell->cells->diseaseCount[sim_cell->cellIdx] * info.massPercent),
            };
            simEvents->unstableCellInfo.push_back(infoCell);

            sim_cell->modifyMass(fabs(infoCell.mass));
            sim_cell->ModifyDiseaseCount(-infoCell.diseaseCount);
            simEvents->ChangeSubstance(simData, sim_cell->cellIdx);
            if (sim_cell->mass() <= 0.0)
                Evaporate(sim_cell->cellIdx, simData, simEvents);
        }
    }
    if (resetFlag)
        simData->timers[sim_cell->cellIdx].stableCellTicks |= 0x1F;
}

bool DoPressureBreak(SimData* simData, SimEvents* simEvents, int cell_idx, float pressure, const int ncell_offset)
{
    float pressureResist = 1.0f;
    int layer = 0;
    for (; layer < 3; layer++) {
        int cell = cell_idx + ncell_offset * (layer + 1);
        ElementPostProcessData* data = &gElementPostProcessData[simData->updatedCells->elementIdx[cell]];

        if ((data->state & 3) != 3)                             break;        // Not Solid
        if ((data->state & 4) != 0)                             return false; // Unbreakable
        if ((simData->updatedCells->properties[cell] & 8) != 0) return false; // Unbreakable
        uint8_t strengthInfo = simData->updatedCells->strengthInfo[cell];     // ConstructedTile: 0~127, Nature Tile 132
        float   massFactor   = strengthInfo > 127 ? simData->updatedCells->mass[cell] / data->maxMass : 1.0f;
        pressureResist      += massFactor * data->strength * (strengthInfo & 0x7F) * 0.25f;
        if (pressureResist > pressure)                          return false; // pressure = mass / maxmass
    }
    if (layer == 0 || layer == 3 || pressureResist > pressure) return false;

    for (int i = 0; i < layer; i++) {
        int cell          = cell_idx + ncell_offset * (i + 1);
        int cellInnerDest = CELL_S2G(cell    , simData->width);
        int cellInnerSrc  = CELL_S2G(cell_idx, simData->width);
        if (CELL_AVAILABLE(cellInnerDest, simData) && CELL_AVAILABLE(cellInnerSrc, simData)) {
            WorldDamageInfo info{ .cellIdx = cellInnerDest, .damageSourceCellIdx = cellInnerSrc };
            simEvents->worldDamageInfo.push_back(info);
        }
    }
    return true;
}

bool DoPartialHeatTransition(SimData* simData, SimEvents* simEvents, int transition_cell_idx, const int heat_src_cell_idx)
{
    // Entering Condition: [Trans cell: Liquid]
    // Situation1: [Trans Target: Gas] and {[Top cell: Gas] or [Target cell: Liquid]}
    // Situation2: [Trans Target: Liquid] and [Target cell: Liquid]
    if (simData->updatedCells->mass[transition_cell_idx] < 5.0f) return false;
    uint16_t elemSrc   = simData->updatedCells->elementIdx[heat_src_cell_idx];
    uint16_t elemTarg  = simData->updatedCells->elementIdx[transition_cell_idx];
    uint16_t elemTrans = gElements[elemTarg].highTempTransitionIdx;
    float    tempSrc   = simData->updatedCells->temperature[heat_src_cell_idx];
    float    tempTarg  = simData->updatedCells->temperature[transition_cell_idx];
    float    tempTrans = gElements[elemTarg].highTemp + 3.0f;

    if (elemTrans == 0xFFFF)                               return false; // No high transition 
    if (tempSrc   <= tempTrans)                            return false; // Heat source is not hot enough
    if (tempSrc   <  gElements[elemSrc ].lowTemp  + 13.0f) return false; // Heat source is not hot enough
    if (tempTarg  >= gElements[elemTarg].highTemp -  3.0f) return false; // Transition cell is too hot

    float heatRequire  = (tempTrans - tempTarg) * gElements[elemTarg].specificHeatCapacity * 5.0f;
    float heatCapacity = simData->updatedCells->mass[heat_src_cell_idx] * gElements[elemSrc].specificHeatCapacity;
    if (heatCapacity * 10.0f < heatRequire)                return false; // Heat source can't provide enough heat

    CellAccessor transition_cell(simData->updatedCells.get(), transition_cell_idx);
    int   diseaseMove = (int)(simData->updatedCells->diseaseCount[transition_cell_idx] * 5.0f / simData->updatedCells->mass[transition_cell_idx]);
    float tempFinal   = (heatCapacity * tempSrc - heatRequire) / heatCapacity;

    switch (gElements[elemTrans].state & 3) {
    case 1: { // Transition target is gas
        int      cellTop = simData->width + transition_cell_idx;
        uint16_t elemTop = simData->updatedCells->elementIdx[cellTop];

        // Remove mass from transition cell
        simData->updatedCells->mass[transition_cell_idx] -= 5.0f;
        transition_cell.ModifyDiseaseCount(-diseaseMove);
        simData->updatedCells->temperature[heat_src_cell_idx] = tempFinal;

        // If top cell is gas, try to remove gas on top cell
        if ((gElements[elemTop].state & 3) == 1 && DisplaceGas(simData, simEvents, simData->updatedCells.get(), cellTop, elemTop)) {
            simData->updatedCells->mass       [cellTop] = 5.0f;
            simData->updatedCells->temperature[cellTop] = tempTrans;
            simData->updatedCells->elementIdx [cellTop] = gElements[elemTarg].highTempTransitionIdx;
            simEvents->ChangeSubstance(simData, cellTop);
            return true;
        }
        // Try to remove liquid on transition cell
        else if (DisplaceLiquid(simData, simEvents, simData->updatedCells.get(), transition_cell_idx, elemTarg)) {
            simData->updatedCells->mass       [transition_cell_idx] = 5.0;
            simData->updatedCells->temperature[transition_cell_idx] = tempTrans;
            simData->updatedCells->elementIdx [transition_cell_idx] = gElements[elemTarg].highTempTransitionIdx;
            simEvents->ChangeSubstance(simData, transition_cell_idx);
            return true;
        }
        else {
            // Transition fail, Recover mass.
            simData->updatedCells->mass[transition_cell_idx] += 5.0f;
            transition_cell.ModifyDiseaseCount(diseaseMove);
            simData->updatedCells->temperature[heat_src_cell_idx] = tempSrc;
            return false;
        }
    }
    case 2: // Transition target is liquid
        // Remove mass from transition cell
        simData->updatedCells->mass[transition_cell_idx] -= 5.0f;
        transition_cell.ModifyDiseaseCount(-diseaseMove);
        simData->updatedCells->temperature[heat_src_cell_idx] = tempFinal;

        // Following process only triggered at [Liquid -> Liquid]
        if (DisplaceLiquid(simData, simEvents, simData->updatedCells.get(), transition_cell_idx, elemTarg)) {
            if (simEvents->SpawnFallingLiquid(simData, transition_cell_idx, gElements[elemTarg].highTempTransitionIdx, 5.0f,
                tempTrans, simData->updatedCells->diseaseIdx[transition_cell_idx], diseaseMove, simData->debugProperties.isDebugEditing)) {
                return true;
            }
            simData->updatedCells->mass       [transition_cell_idx] = 5.0f;
            simData->updatedCells->temperature[transition_cell_idx] = tempTrans;
            simData->updatedCells->elementIdx [transition_cell_idx] = gElements[elemTarg].highTempTransitionIdx;
            simEvents->ChangeSubstance(simData, transition_cell_idx);
            return true;
        }
        else {
            // Transition fail, Recover mass.
            simData->updatedCells->mass[transition_cell_idx] += 5.0f;
            transition_cell.ModifyDiseaseCount(diseaseMove);
            simData->updatedCells->temperature[heat_src_cell_idx] = tempSrc;
            return false;
        }
    default:
        return false;
    };
}

// Add Liquid(mass = flow, elem = elemIdx, temp <- cell) to ncell
// If ncell is solid / different liquid / unmoved gas, action will fail
bool UpdateNeighbourLiquidMass(SimData* simData, SimEvents* simEvents, int cell, int elemIdx, int ncell, int nelemIdx, float flow, uint8_t disease_idx, int disease_count)
{
    LOGGER_LEVEL(1, "%s-cell:%d->%d, elem:%d, temp:%f, mass:%f\n", __func__, cell, ncell, elemIdx, simData->cells->temperature[cell], flow);

	float temperature = simData->cells->temperature[cell];
	CellSOA* p_cells = simData->updatedCells.get();
	// Source cell and Target cell have same element
	if (elemIdx == simData->updatedCells->elementIdx[ncell]) {
		AddMassAndUpdateTemperature(p_cells, p_cells, ncell, flow, temperature, disease_idx, disease_count);
		return true;
	}
	// Target cell is void. Moved mass will be vanished
	if (simData->updatedCells->elementIdx[ncell] == simData->voidElementIdx) {
		simData->updatedCells->mass[ncell] = 0.0;
		return true;
	}
	// Target cell is vacuum or (gas and displaced)
    uint8_t state = gElementLiquidData[simData->updatedCells->elementIdx[ncell]].state & 3;
#ifndef __DEBUGED__
	if (state != 2 && (state != 1 || DisplaceGas(simData, simEvents, p_cells, ncell, simData->updatedCells->elementIdx[ncell]))) {
#else
	if (state == 0 || (state == 1 && DisplaceGas(simData, simEvents, p_cells, ncell, simData->updatedCells->elementIdx[ncell]))) {
#endif
		simData->updatedCells->elementIdx[ncell] = elemIdx;
		AddMassAndUpdateTemperature(p_cells, p_cells, ncell, flow, temperature, disease_idx, disease_count);
		simEvents->ChangeSubstance(simData, ncell);
		return true;
	}
	return false;
}

bool IsLiquidPermeable(const CellSOA* cells, uint64_t cell_idx)
{
	if (cells->properties[cell_idx] & 2) //LiquidImpermeable
		return false;

	return (gElementLiquidData[cells->elementIdx[cell_idx]].state & 3) < 2; // vacuum or gas
}

bool IsSolid(const CellSOA* cells, uint64_t cell_idx)
{
	if ((gElementLiquidData[cells->elementIdx[cell_idx]].state & 3) == 3) // Solid
		return true;
	return !(cells->properties[cell_idx] & 2); // LiquidImpermeable
}

//----- (Decompiled Main) --------------------------------------------
float DoGasPressureDisplacement(uint16_t elem_idx, const int cell, const int ncell, const int nncell, SimData* simData)
{
    uint16_t nelem_idx = simData->cells->elementIdx[ncell];
    if (  ((gElementPressureData[nelem_idx].state & 3) == 1)
        && ( simData->updatedCells->elementIdx[cell  ] == elem_idx)
        && ( simData->updatedCells->elementIdx[ncell ] == nelem_idx)
        && ( simData->updatedCells->mass      [ncell ] >  0)
        && ((simData->updatedCells->properties[nncell] &  1) == 0)
#ifdef __SIMDLL_PLUS__ // Introduce gas molar volume
        && (simData->cells->mass[cell] * gElementPressureData[elem_idx].molarVolume > simData->cells->mass[ncell] * gElementPressureData[nelem_idx].molarVolume * gGasDisplace))
#else
        && (simData->cells->mass[cell]                                              > simData->cells->mass[ncell]))
#endif
    {
        if (   (simData->updatedCells->elementIdx[nncell] == nelem_idx)
            || (simData->updatedCells->elementIdx[nncell] == simData->vacuumElementIdx))
        {
            DoDisplacement(simData, simData->simEvents.get(), simData->updatedCells.get(), ncell, nncell);
#ifdef __SIMDLL_PLUS__ // [DoGasPressureDisplacement only] happen on 4 direction. Increase gas flow ratio to avoid element change.
            float mass_delta = MIN_F(simData->updatedCells->mass[cell], simData->cells->mass[cell] * 0.2f);
#else
            float mass_delta = MIN_F(simData->updatedCells->mass[cell], simData->cells->mass[cell] * 0.125f);
#endif
            simData->updatedCells->mass       [ncell] += mass_delta;
            simData->updatedCells->temperature[ncell]  = simData->cells->temperature[cell];
            simData->updatedCells->elementIdx [ncell]  = elem_idx;
            simData->simEvents->ChangeSubstance(simData, ncell);
            float mass_tmp = simData->updatedCells->mass[cell] - mass_delta;
            simData->updatedCells->mass[cell] = MAX_F(mass_tmp, 0.0f);
            int disease_delta = simData->cells->diseaseCount[cell] >> 3; // (int)(simData->cells->diseaseCount[cell] * 0.125);
            gDisease->AddDiseaseToCell(simData->updatedCells.get(), ncell, simData->cells->diseaseIdx[cell], disease_delta);
            simData->updatedCells->ModifyDiseaseCount(cell, -disease_delta);
            return mass_delta;
        }
    }
    return 0.0f;
}

float DoLiquidPressureDisplacement(uint16_t elem_idx, const int cell, const int ncell, const int nncell, SimData* simData)
{
    if (  ((simData->updatedCells->properties[ncell] & 2) == 0)
        && (simData->cells       ->elementIdx[ncell] != elem_idx)
        && (simData->updatedCells->elementIdx[ cell] == elem_idx)
        && (simData->updatedCells->elementIdx[ncell] == simData->cells->elementIdx[ncell])
        && ((gElementPressureData[simData->cells->elementIdx[ncell]].state & 3) == 2))
    {
        float mass_org = simData->cells->mass[cell] - (simData->flow[cell].y + simData->flow[cell].x + simData->flow[cell].z + simData->flow[cell].w);
#ifdef __DEBUGED__
        mass_org = MIN_F(simData->updatedCells->mass[cell] * 8, mass_org);
#endif
        if (mass_org > (simData->cells->mass[nncell] + simData->cells->mass[ncell]) && DisplaceLiquidDirectional(simData, simData->simEvents.get(), ncell, nncell)) {
            float mass_delta = mass_org * 0.125f;
            simData->updatedCells->mass       [ncell] += mass_delta;
            simData->updatedCells->temperature[ncell]  = simData->cells->temperature[cell];
            simData->updatedCells->elementIdx [ncell]  = elem_idx;
            simData->simEvents->ChangeSubstance(simData, ncell);
            float mass_tmp = simData->updatedCells->mass[cell] - mass_delta;
            simData->updatedCells->mass[cell] = MAX_F(mass_tmp, 0.0f);
            int disease_delta = simData->cells->diseaseCount[cell] >> 3; // (int)(simData->cells->diseaseCount[cell] * 0.125);
            gDisease->AddDiseaseToCell(simData->updatedCells.get(), ncell, simData->cells->diseaseIdx[cell], disease_delta);
            simData->updatedCells->ModifyDiseaseCount(cell, -disease_delta);
            return mass_delta;
        }
    }

    return 0;
}

void PostProcessCell(SimData* simData, SimEvents* simEvents, int cell_idx)
{
#define BIG_GAS(_cell) (simData->updatedCells->mass[_cell] >= 1.0f && (gElementPostProcessData[simData->updatedCells->elementIdx[_cell]].state & 3) == 1)

    ElementPostProcessData* data = &gElementPostProcessData[simData->updatedCells->elementIdx[cell_idx]];
    if ((data->state & 3) == 0) { // Vacuum
        if ((simData->updatedCells->properties[cell_idx] & 1) == 0) // Gas Permeable
            DoSublimation(simData, simEvents, cell_idx, data);
        return;
    }
    if ((data->state & 3) == 3) { // Solid
        if ((data->state & 0xB) == 0xB) { // Unstable Solid
            CellAccessor neighbours(simData->updatedCells.get(), cell_idx);
            if (simData->savedOptions & 1)
                DoUnstableCheckWithDiagonals(simData, simEvents, &neighbours);
            else
                DoUnstableCheckBasic(simData, simEvents, &neighbours);
        }
        return;
    }

    // Gas
    if ((data->state & 3) == 1) {
        displacement_offset = -displacement_offset;
        int cell_idx_N = cell_idx - displacement_offset;
        int cell_idx_P = cell_idx + displacement_offset;
        int cell_idx_D = cell_idx - simData->width;
        int cell_idx_T = cell_idx + simData->width;

        // Gas elimination condition:
        // Mass < 1ug
        // Mass < 1g  And (one of surrounding cell is Gas And Mass >= 1)
        if (simData->updatedCells->mass[cell_idx] < 1e-9 || simData->updatedCells->mass[cell_idx] < 0.001
            && (BIG_GAS(cell_idx_N) || BIG_GAS(cell_idx_P) || BIG_GAS(cell_idx_D) || BIG_GAS(cell_idx_T)))
        {
            simData->flow[cell_idx].w = 0.0;
            simData->flow[cell_idx].z = 0.0;
            simData->flow[cell_idx].y = 0.0;
            simData->flow[cell_idx].x = 0.0;
            Evaporate(cell_idx, simData, simEvents);
        }
        if (!DoDensityDisplacement(simData, simEvents, cell_idx, data, 0.99f)) {
            if ((simData->updatedCells->properties[cell_idx] & 2) // LiquidImpermeable
                || !DoPartialMelt(simData, simEvents, cell_idx, cell_idx_D)
                && !DoPartialMelt(simData, simEvents, cell_idx, cell_idx - 1)
                && !DoPartialMelt(simData, simEvents, cell_idx, cell_idx + 1)
                && !DoPartialMelt(simData, simEvents, cell_idx, cell_idx_T))
            {
                // Random swap with left/right/bottom cell
                // Hot cell has more chance
                // Bottom cell must lighter than target cell
                if (RAND_FLOAT_01(simData->randomSeed) > 0.9) {
                    int   cell_swap = -1;
                    float temp_max = FLT_MAX;
                    int   cell_List[3] = { cell_idx_P, cell_idx_N, cell_idx_D };
                    bool  check_Molar[3] = { false, false, true };
                    for (int i = 0; i < 3; i++) {
                        int cell = cell_List[i];
                        if (simData->updatedCells->temperature[cell] <= temp_max && RAND_FLOAT_01(simData->randomSeed) > 0.5) {
                            Element* elemenet = &gElements[simData->updatedCells->elementIdx[cell]];
                            if ((elemenet->state & 3) != 1)                               continue;
                            if (check_Molar[i] && data->molarMass <= elemenet->molarMass) continue;
                            cell_swap = cell;
                            temp_max  = simData->updatedCells->temperature[cell];
                        }
                    }
                    if (cell_swap != -1) {
                        simData->updatedCells->SwapCells(cell_swap, cell_idx);
                        simEvents->ChangeSubstance(simData, cell_swap);
                        simEvents->ChangeSubstance(simData, cell_idx);
                    }
                }
            }
        }
        DoSublimation(simData, simEvents, cell_idx, data);
        return;
    }

    // Liquid
    if (simData->updatedCells->mass[cell_idx] < 0.01) {
        simData->flow[cell_idx].w = 0.0;
        simData->flow[cell_idx].z = 0.0;
        simData->flow[cell_idx].y = 0.0;
        simData->flow[cell_idx].x = 0.0;
        Evaporate(cell_idx, simData, simEvents);
        return;
    }

    if (DoDensityDisplacement(simData, simEvents, cell_idx, data, 0.3f))
        return;

    // High pressure liquid
    if (simData->updatedCells->mass[cell_idx] > data->maxMass) {
        // Break tiles
        float pressure = simData->updatedCells->mass[cell_idx] / data->maxMass;
        if (DoPressureBreak(simData, simEvents, cell_idx, pressure, -simData->width)) return;
        if (DoPressureBreak(simData, simEvents, cell_idx, pressure, -1             )) return;
        if (DoPressureBreak(simData, simEvents, cell_idx, pressure,  1             )) return;
        if (DoPressureBreak(simData, simEvents, cell_idx, pressure,  simData->width)) return;
#ifdef __SIMDLL_PLUS__ // Stacked liquid mass: Fixed 1.01 -> Variable
        if (simData->updatedCells->mass[cell_idx] > data->maxMass * data->maxCompression) {
            // Try to move liquid in top cell
            int      cell_Top = cell_idx + simData->width;
            float    mass_Top = simData->updatedCells->mass      [cell_Top];
            uint16_t elem_Top = simData->updatedCells->elementIdx[cell_Top];

            if (simData->updatedCells->mass[cell_idx] > mass_Top 
                && (gElementLiquidData[elem_Top].state & 3) == 2
                && simData->updatedCells->mass[cell_idx] > MAX_F(mass_Top * data->compression, data->maxMass) + data->minHorizontalFlow
                && DisplaceLiquid(simData, simEvents, simData->updatedCells.get(), cell_Top, elem_Top)) 
            {
                float mass_Final    =       simData->updatedCells->mass        [cell_idx] / (1 + data->compression);
                int   disease_Final = (int)(simData->updatedCells->diseaseCount[cell_idx] / (1 + data->compression));
#else
        if (simData->updatedCells->mass[cell_idx] > data->maxMass * 1.5) {
            // Try to move liquid in top cell
            int      cell_Top = cell_idx + simData->width;
            float    mass_Top = simData->updatedCells->mass      [cell_Top];
            uint16_t elem_Top = simData->updatedCells->elementIdx[cell_Top];

            if (simData->updatedCells->mass[cell_idx] > mass_Top 
                && (gElementLiquidData[elem_Top].state & 3) == 2
#ifndef __DEBUGED__
                && simData->updatedCells->mass[cell_idx] > MAX_F(mass_Top * 1.01f, data->maxMass) + data->minHorizontalFlow
#else
                && simData->updatedCells->mass[cell_idx] > MAX_F(mass_Top * 1.01f, data->maxMass)
#endif
                && DisplaceLiquid(simData, simEvents, simData->updatedCells.get(), cell_Top, elem_Top)) 
            {
                float mass_Final    =       simData->updatedCells->mass        [cell_idx] / 2.01f;
                int   disease_Final = (int)(simData->updatedCells->diseaseCount[cell_idx] / 2.01f);
#endif
                simData->updatedCells->mass        [cell_Top] = mass_Final;
                simData->updatedCells->temperature [cell_Top] = simData->updatedCells->temperature[cell_idx];
                simData->updatedCells->elementIdx  [cell_Top] = simData->updatedCells->elementIdx [cell_idx];
                simData->updatedCells->diseaseIdx  [cell_Top] = simData->updatedCells->diseaseIdx [cell_idx];
                simData->updatedCells->diseaseCount[cell_Top] = disease_Final;
                simData->updatedCells->diseaseInfestationTickCount  [cell_Top] = 0;
                simData->updatedCells->diseaseGrowthAccumulatedError[cell_Top] = 0;
                simEvents->ChangeSubstance(simData, cell_Top);
                simData->updatedCells->mass[cell_idx] -= mass_Final;

                CellAccessor neighbours(simData->updatedCells.get(), cell_idx);
                neighbours.ModifyDiseaseCount(-disease_Final);

                return;
            }
        }
    }

    if (DoPartialHeatTransition(simData, simEvents, cell_idx, cell_idx - simData->width)) return;
    if (DoPartialHeatTransition(simData, simEvents, cell_idx, cell_idx - 1             )) return;
    if (DoPartialHeatTransition(simData, simEvents, cell_idx, cell_idx + 1             )) return;
    if (DoPartialHeatTransition(simData, simEvents, cell_idx, cell_idx + simData->width)) return;

    int cell_Top = cell_idx + simData->width;
    if (simData->updatedCells->mass[cell_Top] >= 1.8f)                    return;
    if (data->sublimateIndex == 0xFFFF)                                   return;
    if ((gElementPostProcessData[simData->updatedCells->elementIdx[cell_Top]].state & 3) > 1) return;
    if (simData->updatedCells->properties[cell_Top] & 1)                  return;
    if (RAND_FLOAT_01(simData->randomSeed) >= data->sublimateProbability) return;

    int      cell_Bottom = cell_idx - simData->width;
    uint8_t  diseIdx_C   = simData->updatedCells->diseaseIdx  [cell_idx];
    uint8_t  diseCnt_C   = simData->updatedCells->diseaseCount[cell_idx];
    uint16_t elem_Top    = simData->updatedCells->elementIdx  [cell_Top];
    float    mass_Bottom = simData->updatedCells->mass        [cell_Bottom];
    float    mass_Center = simData->updatedCells->mass        [cell_idx];
    float    mass_Source = mass_Center;

    if (simData->updatedCells->elementIdx[cell_Bottom] == simData->updatedCells->elementIdx[cell_idx] && mass_Bottom > mass_Center) {
        mass_Source = mass_Bottom;
        if (simData->updatedCells->diseaseIdx[cell_Bottom] == diseIdx_C)
            diseCnt_C += (int)(simData->updatedCells->diseaseCount[cell_Bottom] * mass_Center / mass_Bottom);
    }

    // Evaporation
    float mass_Evap = MIN_F(mass_Source * data->offGasPercentage, 1.0);
    if (mass_Center - mass_Evap >= 0.01) {
        if (mass_Evap <= 0.0) return;

        if (elem_Top == data->sublimateIndex) {
            float mass_Dest = MIN_F(mass_Evap * data->sublimateEfficiency, 1.8F - simData->updatedCells->mass[cell_Top]);
            mass_Evap       = mass_Dest / data->sublimateEfficiency;

            AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), cell_Top, mass_Dest,
                simData->updatedCells->temperature[cell_idx], diseIdx_C, (int)(mass_Evap / mass_Center * diseCnt_C));
            simData->updatedCells->mass[cell_idx] = MAX_F(mass_Center - mass_Evap, 0.0f);
            simData->accumulatedFlow   [cell_idx] = mass_Evap + simData->accumulatedFlow[cell_idx];
            simEvents->SpawnFX(simData, cell_idx, data->sublimateFX, 0.0);
        }
        else {
            if (elem_Top != simData->vacuumElementIdx && !DisplaceGas(simData, simEvents, simData->updatedCells.get(), cell_Top, elem_Top))
                return;

            simData->updatedCells->elementIdx  [cell_Top] = data->sublimateIndex;
            simData->updatedCells->mass        [cell_Top] = mass_Evap * data->sublimateEfficiency;
            simData->updatedCells->temperature [cell_Top] = simData->updatedCells->temperature[cell_idx];
            simData->updatedCells->diseaseIdx  [cell_Top] = simData->updatedCells->diseaseIdx [cell_idx];
            simData->updatedCells->diseaseCount[cell_Top] = (int)(mass_Evap / mass_Center * diseCnt_C);
            simEvents->ChangeSubstance(simData, cell_Top);

            simEvents->SpawnFX(simData, cell_idx, data->sublimateFX, 0.0);
            simData->updatedCells->mass[cell_idx] = MAX_F(mass_Center - mass_Evap, 0.0f);
            simData->accumulatedFlow   [cell_idx] = mass_Evap + simData->accumulatedFlow[cell_idx];
        }
    }
    else {
        simData->updatedCells->elementIdx  [cell_idx] = data->sublimateIndex;
        simData->updatedCells->mass        [cell_idx] = mass_Evap * data->sublimateEfficiency;
        simData->updatedCells->temperature [cell_idx] = simData->updatedCells->temperature[cell_idx];
        simData->updatedCells->diseaseCount[cell_idx] = (int)(diseCnt_C * MAX_F(mass_Evap / mass_Center, 1.0f));

        if (simData->updatedCells->elementIdx[cell_Bottom] != simData->vacuumElementIdx) {
            if (mass_Evap > mass_Center) {
                int diseCnt_Move = (int)((mass_Evap - mass_Center) / mass_Bottom * simData->updatedCells->diseaseCount[cell_Bottom]);
                gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell_idx, simData->updatedCells->diseaseIdx[cell_Bottom], diseCnt_Move);
                CellAccessor neighbours(simData->updatedCells.get(), cell_idx - simData->width);
                neighbours.ModifyDiseaseCount(-diseCnt_Move);
            }
            if (mass_Evap != mass_Center)
                simData->updatedCells->temperature[cell_idx] = CalculateFinalTemperature(
                    mass_Center, simData->updatedCells->temperature[cell_idx], mass_Evap - mass_Center, simData->updatedCells->temperature[cell_Bottom]);
        }
        simEvents->ChangeSubstance(simData, cell_idx);
        simEvents->SpawnFX(simData, cell_idx, data->sublimateFX, 0.0);
    }

#ifndef __DEBUGED__
    if (simData->updatedCells->elementIdx[cell_Bottom] == simData->updatedCells->elementIdx[cell_idx] && mass_Center < mass_Evap)
#else
    if (simData->updatedCells->elementIdx[cell_Bottom] == simData->cells->elementIdx[cell_idx] && mass_Center < mass_Evap)
#endif
        simData->updatedCells->mass      [cell_Bottom] -= mass_Evap - mass_Center;

    return;
}

float UpdatePressure(SimData* simData, SimEvents* simEvents, CellSOA* cells, CellSOA* updated_cells,
	int cell, uint16_t elem_idx, uint8_t elem_state, float elem_flow, int ncell)
{
//    LOGGER_PRINT2("%s: %d->%d\n", __func__, cell, ncell);
	float flow = elem_state == 1 ? elem_flow : gElementPressureData[cells->elementIdx[ncell]].flow;
	float massC2N = CLAMP_F(flow * (cells->mass[cell] - cells->mass[ncell]), cells->mass[cell] * 0.125f, cells->mass[ncell] * -0.125f);

	int      cellSrc, cellDest;
	uint16_t elemSrc, elemDest;
    uint8_t  stateDest;
    if (massC2N > 1e-10) {
		//Flow in
		elemSrc   = elem_idx;
		elemDest  = cells->elementIdx[ncell];
		stateDest = gElementPressureData[cells->elementIdx[ncell]].state & 3;
		cellSrc   = cell;
		cellDest  = ncell;
	}
	else if (massC2N < -1e-10) {
		//Flow out
		elemSrc   = cells->elementIdx[ncell];
		elemDest  = elem_idx;
		stateDest = elem_state;
		cellSrc   = ncell;
		cellDest  = cell;
		cells->elementIdx[ncell] = elem_idx;
	}
    else return 0; // Improve performance, Gas will be delete at 1e-9
	float massMove    = MIN_F(cells->mass[cellSrc], fabsf(massC2N));
	int   diseaseMove = (int)(CLAMP_F(massMove / cells->mass[cellSrc], 1, 0) * cells->diseaseCount[cellSrc]);

	//Update destination cell
	updated_cells->temperature[cellDest] = CalculateFinalTemperature(massMove, updated_cells->temperature[cellSrc],
		updated_cells->mass[cellDest], updated_cells->temperature[cellDest]);

	updated_cells->mass[cellDest] += massMove;
    gDisease->AddDiseaseToCell(updated_cells, cellDest, updated_cells->diseaseIdx[cellSrc], diseaseMove);
	
    //Update source cell
	updated_cells->mass[cellSrc] = MAX_F(updated_cells->mass[cellSrc] - massMove, 0);
	updated_cells->diseaseCount[cellSrc] -= diseaseMove;
	if (updated_cells->diseaseCount[cellSrc] <= 0)
		updated_cells->ClearDisease(cellSrc);
	if (stateDest) return massC2N;

	if (elemDest == simData->voidElementIdx) {
		//Destination cell is Void. Material will be eliminated by Void 
		updated_cells->mass       [cellDest] = 0;
		updated_cells->temperature[cellDest] = 0;
		updated_cells->ClearDisease(cellDest);
	}
	else {
		updated_cells->elementIdx[cellDest] = elemSrc;
		simEvents->ChangeSubstance(simData, cellDest);
	}

//    LOGGER_PRINT2("%s done\n", __func__);
	return massC2N;
}

void UpdateLiquid(SimData* simData, SimEvents* simEvents, int cell_idx)
{
	float massOrg_C = simData->cells->mass[cell_idx];
	int   elemIdx   = simData->cells->elementIdx[cell_idx];
    LOGGER_LEVEL(1, "%s, cell %d, elem %s, mass %f\n", __func__, cell_idx, gElementNames[elemIdx].c_str(), massOrg_C);

	// Bottom cell
	int cellBottom = cell_idx - simData->width;
	int elemIdxB   = simData->cells->elementIdx[cellBottom];
	// Not Solid, not impermeable
	if (((gElementLiquidData[simData->updatedCells->elementIdx[cellBottom]].state & 3) != 3) && ((simData->cells->properties[cellBottom] & 2) == 0)) {
		// [viscosity] stands for [speed] in element table
		// $maxMass_Stack = max(maxMass, compression * mass_upLayer)
		// $mass_Flow = ($maxMass_Stack - mass_Center) / 2; $mass_Flow > minVerticalFlow
		if (elemIdx == elemIdxB) {
			float massOrg_B = 0;
			if ((gElementLiquidData[elemIdxB].state & 3) == 2)
				massOrg_B = simData->cells->mass[cellBottom];
#ifdef __SIMDLL_PLUS__ // Stacked liquid mass: Fixed 1.01 -> Variable
			float massMax   = MAX_F(massOrg_C * gElementLiquidData[elemIdx].compression, gElementLiquidData[elemIdx].maxMass);
#else
            float massMax   = MAX_F(massOrg_C * 1.01f, gElementLiquidData[elemIdx].maxMass);
#endif
			float viscosity = MIN_F(gElementLiquidData[elemIdx].viscosity, massOrg_C);
			float massFlow  = CLAMP_F((massMax - massOrg_B) * 0.5f, viscosity, 0);
			float minVerticalFlow = gElementLiquidData[elemIdx].minVerticalFlow;
			if ((massOrg_C <= minVerticalFlow || massFlow >= minVerticalFlow) && massFlow > 0.0) {
				int diseaseMove = (int)(massFlow / massOrg_C * simData->cells->diseaseCount[cell_idx]);

				if (UpdateNeighbourLiquidMass(simData, simEvents, cell_idx, elemIdx, cellBottom, elemIdxB, massFlow, simData->cells->diseaseIdx[cell_idx], diseaseMove)) {
					simData->cells->diseaseCount[cell_idx] -= diseaseMove;
					if (simData->updatedCells->diseaseIdx[cell_idx] == simData->cells->diseaseIdx[cell_idx])
						simData->updatedCells->ModifyDiseaseCount(cell_idx, -diseaseMove);

					massOrg_C                             -= massFlow;
					simData->updatedCells->mass[cell_idx] -= massFlow;
					simData->flow[cell_idx].z             += massFlow;
				}
			}
		}
		// Bottom cell is Gas / Vacuum
		else if ((gElementLiquidData[elemIdxB].state & 3) <= 1) {
			// Bottom cell is void: Clear center cell
			if (elemIdxB == simData->voidElementIdx) {
				simData->updatedCells->elementIdx [cell_idx] = simData->vacuumElementIdx;
				simData->updatedCells->temperature[cell_idx] = 0.0;
				simData->updatedCells->mass       [cell_idx] = 0.0;
				simData->updatedCells->ClearDisease(cell_idx);
				simEvents->ChangeSubstance(simData, cell_idx);
				return;
			}
			// Condition: L is Permeable and LB is not /  R is Permeable and RB is not 
			else if ((IsLiquidPermeable(simData->updatedCells.get(), cell_idx - 1) && !IsLiquidPermeable(simData->updatedCells.get(), cellBottom - 1))
				  || (IsLiquidPermeable(simData->updatedCells.get(), cell_idx + 1) && !IsLiquidPermeable(simData->updatedCells.get(), cellBottom + 1))) {
				// SpawnFallingLiquid Fail, do nothing
				if (!simEvents->SpawnFallingLiquid(simData, cell_idx, elemIdx, massOrg_C, simData->cells->temperature[cell_idx],
					simData->cells->diseaseIdx[cell_idx], simData->cells->diseaseCount[cell_idx], simData->debugProperties.isDebugEditing))
					return;
				// SpawnFallingLiquid success, clear center cell
				simData->updatedCells->elementIdx [cell_idx] = simData->vacuumElementIdx;
				simData->updatedCells->temperature[cell_idx] = 0.0;
				simData->updatedCells->mass       [cell_idx] = 0.0;
				simData->updatedCells->ClearDisease(cell_idx);
				simEvents->ChangeSubstance(simData, cell_idx);
				return;
			}
			else {
				// Swap with bottom cell
				simData->updatedCells->SwapCells(cell_idx, cellBottom);
				simEvents->ChangeSubstance(simData, cell_idx);
				simEvents->ChangeSubstance(simData, cellBottom);
				return;
			}
		}
	}
	if (massOrg_C <= 0.0) return;

	// Left cell
	int cellLeft = cell_idx - 1;
	int elemIdxL = simData->cells->elementIdx[cellLeft];
	// Same liquid or gas or vacuum
	if ((elemIdx == elemIdxL || (gElementLiquidData[elemIdxL].state & 3) <= 1)
		&& (gElementLiquidData[simData->updatedCells->elementIdx[cellLeft]].state & 3) != 3
		&& (simData->cells->properties[cellLeft] & 2) == 0)
	{
		float massOrg_L = 0;
		if ((gElementLiquidData[elemIdxL].state & 3) == 2)
			massOrg_L = simData->cells->mass[cellLeft];

		// mass_Flow = (mass_Center - mass_Left) / 4; mass_Flow > minHorizontalFlow
		float viscosity = MIN_F(gElementLiquidData[elemIdx].viscosity, massOrg_C);
		float massFlow  = MIN_F(viscosity, (massOrg_C - massOrg_L) * 0.25f);
		if (massFlow >= gElementLiquidData[elemIdx].minHorizontalFlow && massFlow > 0) {
			int diseaseMove = (int)(massFlow / massOrg_C * simData->cells->diseaseCount[cell_idx]);

			bool updatFlag = false;
			// SpawnFallingLiquid Condition: Bottom is Solid, Left Bottom is Vacuum/Gas
			if (IsSolid(simData->cells.get(), cell_idx - simData->width)
				&& IsLiquidPermeable(simData->updatedCells.get(), cellLeft - simData->width)
				&& simEvents->SpawnFallingLiquid(simData, cellLeft, elemIdx, massFlow, simData->cells->temperature[cell_idx], simData->cells->diseaseIdx[cell_idx], diseaseMove, simData->debugProperties.isDebugEditing))
				updatFlag = true;
			else if (UpdateNeighbourLiquidMass(simData, simEvents, cell_idx, elemIdx, cellLeft, elemIdxL, massFlow, simData->cells->diseaseIdx[cell_idx], diseaseMove))
				updatFlag = true;

			if (updatFlag) {
				simData->cells->diseaseCount[cell_idx] -= diseaseMove;
				if (simData->updatedCells->diseaseIdx[cell_idx] == simData->cells->diseaseIdx[cell_idx])
					simData->updatedCells->ModifyDiseaseCount(cell_idx, -diseaseMove);

				massOrg_C                             -= massFlow;
				simData->updatedCells->mass[cell_idx] -= massFlow;
				simData->flow[cell_idx].x             += massFlow;
			}
		}
	}
	if (massOrg_C <= 0.0) return;

	// Right cell
	int cellRight = cell_idx + 1;
	int elemIdxR  = simData->cells->elementIdx[cellRight];
	if ((elemIdx == elemIdxR || (gElementLiquidData[elemIdxR].state & 3) <= 1)
		&& (gElementLiquidData[simData->updatedCells->elementIdx[cellRight]].state & 3) != 3
		&& (simData->cells->properties[cellRight] & 2) == 0) 
	{
		float massOrg_R = 0;
		if ((gElementLiquidData[elemIdxR].state & 3) == 2)
			massOrg_R = simData->cells->mass[cellRight];

		float viscosity = MIN_F(gElementLiquidData[elemIdx].viscosity, massOrg_C);
		float massFlow  = MIN_F(viscosity, (massOrg_C - massOrg_R) * 0.25f);
		if (massFlow >= gElementLiquidData[elemIdx].minHorizontalFlow && massFlow > 0) {
			int diseaseMove = (int)(massFlow / massOrg_C * simData->cells->diseaseCount[cell_idx]);

			bool updatFlag = false;
			if (IsSolid(simData->cells.get(), cell_idx - simData->width) 
				&& IsLiquidPermeable(simData->updatedCells.get(), cellRight - simData->width)
				&& simEvents->SpawnFallingLiquid(simData, cellRight, elemIdx, massFlow, simData->cells->temperature[cell_idx], simData->cells->diseaseIdx[cell_idx], diseaseMove, simData->debugProperties.isDebugEditing))
				updatFlag = true;
			else if (UpdateNeighbourLiquidMass(simData, simEvents, cell_idx, elemIdx, cellRight, elemIdxR, massFlow, simData->cells->diseaseIdx[cell_idx], diseaseMove))
				updatFlag = true;

			if (updatFlag) {
				simData->cells->diseaseCount[cell_idx] -= diseaseMove;
				if (simData->updatedCells->diseaseIdx[cell_idx] == simData->cells->diseaseIdx[cell_idx])
					simData->updatedCells->ModifyDiseaseCount(cell_idx, -diseaseMove);

				massOrg_C                             -= massFlow;
				simData->updatedCells->mass[cell_idx] -= massFlow;
				simData->flow[cell_idx].y             += massFlow;
			}
		}
	}
	if (massOrg_C <= 0.0) return;

	// Top cell
	int cellTop  = cell_idx + simData->width;
	int elemIdxT = simData->cells->elementIdx[cellTop];
	// Same liquid or gas or vacuum
	if ((elemIdx == elemIdxT || (gElementLiquidData[elemIdxT].state & 3) <= 1)
		&& (gElementLiquidData[simData->updatedCells->elementIdx[cellTop]].state & 3) != 3
		&& (simData->cells->properties[cellTop] & 2) == 0)
	{
		float massOrg_T = 0;
		if ((gElementLiquidData[elemIdxT].state & 3) == 2)
			massOrg_T = simData->cells->mass[cellTop];
		// $mass_Flow = (mass_Center - maxMass_Stack) / 2; mass_Flow > 0.01
#ifdef __SIMDLL_PLUS__ // Stacked liquid mass: Fixed 1.01 -> Variable
        float massMax  = MAX_F(massOrg_T * gElementLiquidData[elemIdx].compression, gElementLiquidData[elemIdx].maxMass);
#else
        float massMax  = MAX_F(massOrg_T * 1.01f, gElementLiquidData[elemIdx].maxMass);
#endif
		float massFlow = MAX_F(massOrg_C - massMax, 0.0) * 0.5f;
		if (massOrg_C < massMax * 2) {
			massFlow = MIN_F(massFlow, massOrg_C);
			massFlow = MIN_F(massFlow, gElementLiquidData[elemIdx].viscosity);
		}
		if (massFlow > 0.01) {
			int diseaseMove = (int)(massFlow / massOrg_C * simData->cells->diseaseCount[cell_idx]);
			if (UpdateNeighbourLiquidMass(simData, simEvents, cell_idx, elemIdx, cellTop, elemIdxT, massFlow, simData->cells->diseaseIdx[cell_idx], diseaseMove)) {
				if (simData->updatedCells->diseaseIdx[cell_idx] == simData->cells->diseaseIdx[cell_idx])
					simData->updatedCells->ModifyDiseaseCount(cell_idx, -diseaseMove);

				massOrg_C                             -= massFlow;
				simData->updatedCells->mass[cell_idx] -= massFlow;
				simData->flow[cell_idx].w             += massFlow;
			}
		}
	}
	if (massOrg_C <= 0.0) return;

	// Bottom Right cell
	int cellBR = cell_idx - simData->width + 1;
#ifndef __DEBUGED__
	// BottomRight gas, Right solid / impermeable, Bottom not solid: Swap Center <-> BottomRight
	if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBR]].state & 3) == 1) {
#else
	if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBR]].state & 3) <= 1) {
#endif
		if ((gElementLiquidData[simData->updatedCells->elementIdx[cellRight]].state & 3) == 3
			|| (simData->updatedCells->properties[cellRight] & 1) != 0)
		{
			if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBottom]].state & 3) != 3
				&& (simData->updatedCells->properties[cell_idx]   & 1) == 0
				&& (simData->updatedCells->properties[cellBR]     & 2) == 0
				&& (simData->updatedCells->properties[cellBottom] & 3) == 0)
			{
				simData->updatedCells->SwapCells(cell_idx, cellBR);
				simEvents->ChangeSubstance(simData, cell_idx);
				simEvents->ChangeSubstance(simData, cellBR);
				return;
			}
		}
	}

	// Bottom Left cell
	int cellBL = cell_idx - simData->width - 1;
#ifndef __DEBUGED__
	if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBL]].state & 3) == 1) {
#else
	if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBL]].state & 3) <= 1) {
#endif
		if ((gElementLiquidData[simData->updatedCells->elementIdx[cellLeft]].state & 3) == 3
			|| (simData->updatedCells->properties[cellLeft] & 1) != 0)
		{
			if ((gElementLiquidData[simData->updatedCells->elementIdx[cellBottom]].state & 3) != 3
				&& (simData->updatedCells->properties[cell_idx]   & 1) == 0
				&& (simData->updatedCells->properties[cellBL]     & 2) == 0
				&& (simData->updatedCells->properties[cellBottom] & 3) == 0)
			{
				simData->updatedCells->SwapCells(cell_idx, cellBL);
				simEvents->ChangeSubstance(simData, cell_idx);
				simEvents->ChangeSubstance(simData, cellBL);
				return;
			}
		}
	}
}

void UpdateTemperature(SimData* simData, SimEvents* simEvents, const int cell, const int ncell)
{
	ASSERT_TEMP(simData->cells->mass[cell ], simData->cells->temperature[cell ]);
	ASSERT_TEMP(simData->cells->mass[ncell], simData->cells->temperature[ncell]);

	ElementTemperatureData* tempData_c = &gElementTemperatureData[simData->cells->elementIdx[cell ]];
	ElementTemperatureData* tempData_n = &gElementTemperatureData[simData->cells->elementIdx[ncell]];

	float insulation_c = simData->cells->insulation[cell ] * simData->cells->insulation[cell ] * 0.0000153787f; //1/255/255
	float insulation_n = simData->cells->insulation[ncell] * simData->cells->insulation[ncell] * 0.0000153787f;
	float conduct_c    = insulation_c * tempData_c->thermalConductivity;
	float conduct_n    = insulation_n * tempData_n->thermalConductivity;
	float conduct_f    = MIN_F(conduct_c, conduct_n);
	if (insulation_c >= 1 && insulation_n >= 1) {
		conduct_f = expf((logf(conduct_c) + logf(conduct_n)) * 0.5f);
	}

	//Calculate surfaceAreaMultiplier
	float surfaceAreaMultiplier_c = 1.0f;
	switch (tempData_n->state & 3) {
		case 1: surfaceAreaMultiplier_c = tempData_c->gasSurfaceAreaMultiplier;		break;
		case 2: surfaceAreaMultiplier_c = tempData_c->liquidSurfaceAreaMultiplier;	break;
		case 3: surfaceAreaMultiplier_c = tempData_c->solidSurfaceAreaMultiplier;	break;
	}
	float surfaceAreaMultiplier_n = 1.0f;
	switch (tempData_c->state & 3) {
		case 1: surfaceAreaMultiplier_n = tempData_n->gasSurfaceAreaMultiplier;		break;
		case 2: surfaceAreaMultiplier_n = tempData_n->liquidSurfaceAreaMultiplier;	break;
		case 3: surfaceAreaMultiplier_n = tempData_n->solidSurfaceAreaMultiplier;	break;
	}

	//Calculate heatTransfer value
	float tempOrg_c      = simData->cells->temperature[cell];
	float tempOrg_n      = simData->cells->temperature[ncell];
	float tempDelta      = tempOrg_n - tempOrg_c;
	float heatCapacity_c = simData->cells->mass[cell ] * tempData_c->specificHeatCapacity;
	float heatCapacity_n = simData->cells->mass[ncell] * tempData_n->specificHeatCapacity;
    float heatClamp      = fabsf(tempOrg_n * heatCapacity_n - tempOrg_c * heatCapacity_c);
    float heatTransfer   = fabsf(tempDelta * conduct_f);
	heatTransfer        *= surfaceAreaMultiplier_c * surfaceAreaMultiplier_n * 0.2f;
	heatTransfer         = MIN_F(heatClamp, heatTransfer);

    if (heatTransfer < 0.0001f) return;
	if (tempDelta <= 0) heatTransfer *= -1;

	//Avoid temperature reverse
	float tempFin_c = tempOrg_c + (heatTransfer / heatCapacity_c);
	float tempFin_n = tempOrg_n - (heatTransfer / heatCapacity_n);
	if ((tempFin_n - tempFin_c) * tempDelta < 0) {
		tempFin_n = tempFin_c = (tempOrg_n * heatCapacity_n + tempOrg_c * heatCapacity_c) / (heatCapacity_n + heatCapacity_c);
	}
	float tempMax = MAX_F(tempOrg_n, tempOrg_c);
	float tempMin = MIN_F(tempOrg_n, tempOrg_c);
	float tempDelta_c = CLAMP_F(tempFin_c, tempMax, tempMin) - tempOrg_c;
	float tempDelta_n = CLAMP_F(tempFin_n, tempMax, tempMin) - tempOrg_n;
	tempDelta_c  = CLAMP_F(tempDelta_c, (tempMax - tempOrg_c) * 0.25f, (tempMin - tempOrg_c) * 0.25f);
	tempDelta_n  = CLAMP_F(tempDelta_n, (tempMax - tempOrg_n) * 0.25f, (tempMin - tempOrg_n) * 0.25f);
	heatTransfer = MIN_F(fabsf(tempDelta_n * heatCapacity_n), fabsf(tempDelta_c * heatCapacity_c));

	//Start to use updatedCells
	if (tempDelta_c > 0)	simData->updatedCells->temperature[cell ] += heatTransfer / heatCapacity_c;
	else					simData->updatedCells->temperature[cell ] -= heatTransfer / heatCapacity_c;
	if (tempDelta_n > 0)	simData->updatedCells->temperature[ncell] += heatTransfer / heatCapacity_n;
	else					simData->updatedCells->temperature[ncell] -= heatTransfer / heatCapacity_n;
	//Update new temperature to updatedCells
	simData->updatedCells->temperature[cell]  = CLAMP_F(simData->updatedCells->temperature[cell ], SIM_MAX_TEMPERATURE, 1);
	simData->updatedCells->temperature[ncell] = CLAMP_F(simData->updatedCells->temperature[ncell], SIM_MAX_TEMPERATURE, 1);

	//Final check
	ASSERT_TEMP(simData->updatedCells->mass[cell ], simData->updatedCells->temperature[cell ]);
	ASSERT_TEMP(simData->updatedCells->mass[ncell], simData->updatedCells->temperature[ncell]);
	DoStateTransition(simData, simEvents, cell,  &gElementTemperatureData[simData->updatedCells->elementIdx[cell ]]);
	DoStateTransition(simData, simEvents, ncell, &gElementTemperatureData[simData->updatedCells->elementIdx[ncell]]);
}

void FloodRemoved(SimData* simData, float* mass_to_remove, int16_t remove_elem_idx, ConsumedMassInfo* removed_info, 
    uint16_t x, uint16_t y, int max_depth, std::vector<int>* visited, std::queue<FloodFillInfo>* next)
{
    FloodRemovedInfo info = { .remove_elem_idx = remove_elem_idx , .simData = simData,.mass_to_remove = mass_to_remove, .removed_info = removed_info };
    uint8_t state = gElementPostProcessData[remove_elem_idx].state & 3;
    Flood(simData, x, y, max_depth, state != 3, state == 3, state == 3, visited, next, info);
}