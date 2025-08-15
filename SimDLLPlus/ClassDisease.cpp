#include "global.h"
#include "ClassCell.h"
#include "ClassDisease.h"
#include "ClassSim.h"

//----- (New add) ----------------------------------------------------
int getTempLevel(float temperature, DiseaseInfo::RangeInfo* info)
{
    if (temperature > info->maxViable) return 4;
    if (temperature > info->maxGrowth) return 3;
    if (temperature > info->minGrowth) return 2;
    if (temperature > info->minViable) return 1;
    return 0;
}

float getTempValue(int level, DiseaseInfo::RangeInfo* info)
{
    switch (level) {
        case 0: return info->minViable;
        case 1: return info->minGrowth;
        case 2: return info->maxGrowth;
        case 3: return info->maxViable;
        default:
            ASSERT_TEXT("Assert Unknow Temperature Level");
            return 0;
    }
}

//----- (Decompiled) -------------------------------------------------
DiseaseResult* Disease::CalculateFinalDiseaseCount(DiseaseResult* result, 
    uint8_t src1_idx, int src1_count, uint8_t src2_idx, int src2_count)
{
    result->idx   = -1;
    result->count = 0;
    if (src1_idx == src2_idx) {
        result->idx   = src1_idx;
        result->count = src1_count + src2_count;
        return result;
    }
    if (src1_idx == 0xFF) {
        result->idx   = src2_idx;
        result->count = src2_count;
        return result;
    }
    if (src2_idx == 0xFF) {
        result->idx   = src1_idx;
        result->count = src1_count;
        return result;
    }

    float src1_strength = src1_count * this->diseases[src1_idx].strength;
    float src2_strength = src2_count * this->diseases[src2_idx].strength;
    if (src1_strength <= src2_strength) {
        if (src1_count < 0)
        {
            result->idx   = src2_idx;
            result->count = -src1_count;
            return result;
        }
        result->idx   = src1_idx;
        result->count = src1_count;
        return result;
    }
    int src2sub1_count = src2_count - (int)(src1_strength / src2_strength * src1_count);
    if (src2sub1_count >= 0) {
        result->idx   = src2_idx;
        result->count = src2sub1_count;
        return result;
    }
    result->idx = src1_idx;
    result->count = -src2sub1_count;
    return result;
}

Disease::Disease(BinaryBufferReader* reader)
{
    int num_diseases, num_elements;
    *reader >> num_diseases >> num_elements;
    if (num_elements != gElements.size())
        ASSERT_TEXT("Elements table and diseases growth info entries are dismatch");

    this->diseases    .clear();
    this->diseaseNames.clear();
    this->diseases    .resize(num_diseases);
    this->diseaseNames.resize(num_diseases);

    for (int idx_dise = 0; idx_dise < num_diseases; idx_dise++) {
        DiseaseInfo& info = this->diseases[idx_dise];
        reader->ReadString(&this->diseaseNames[idx_dise]);
        *reader >> info.hashID >> info.strength >> info.temperatureRange 
                >> info.temperatureHalfLives    >> info.pressureRange 
                >> info.pressureHalfLives       >> info.radiationKillRate;

        info.elemGrowthInfo.underPopulationDeathRate        .resize(num_elements);
        info.elemGrowthInfo.populationHalfLife              .resize(num_elements);
        info.elemGrowthInfo.overPopulationHalfLife          .resize(num_elements);
        info.elemGrowthInfo.diffusionScale                  .resize(num_elements);
        info.elemGrowthInfo.minCountPerKG                   .resize(num_elements);
        info.elemGrowthInfo.maxCountPerKG                   .resize(num_elements);
        info.elemGrowthInfo.minDiffusionCount               .resize(num_elements);
        info.elemGrowthInfo.minDiffusionInfestationTickCount.resize(num_elements);
        for (int idx_elem = 0; idx_elem < num_elements; idx_elem++) {
            *reader >> info.elemGrowthInfo.underPopulationDeathRate        [idx_elem]
                    >> info.elemGrowthInfo.populationHalfLife              [idx_elem]
                    >> info.elemGrowthInfo.overPopulationHalfLife          [idx_elem]
                    >> info.elemGrowthInfo.diffusionScale                  [idx_elem]
                    >> info.elemGrowthInfo.minCountPerKG                   [idx_elem]
                    >> info.elemGrowthInfo.maxCountPerKG                   [idx_elem]
                    >> info.elemGrowthInfo.minDiffusionCount               [idx_elem]
                    >> info.elemGrowthInfo.minDiffusionInfestationTickCount[idx_elem];
        }
        LOGGER_PRINT("Disease %s\t, id:%11d, radiationKillRate: %.2f\n", this->diseaseNames[idx_dise].c_str(), info.hashID, info.radiationKillRate);
    }
    LOGGER_PRINT("Create Disease. Count:%d , Element Count%d\n", num_diseases, num_elements);
}

void Disease::AddDiseaseToCell(CellSOA* cells, int cell_idx, uint8_t new_disease_idx, int new_disease_count)
{
//    LOGGER_PRINT2("%s\n", __func__);
    if (new_disease_count == 0)  return; // Improve performance

    int     org_disease_count = cells->diseaseCount[cell_idx];
    int     fin_disease_count = cells->diseaseCount[cell_idx];
    uint8_t org_disease_idx   = cells->diseaseIdx  [cell_idx];
    uint8_t fin_disease_idx   = cells->diseaseIdx  [cell_idx];

    if (org_disease_idx == new_disease_idx) {
        fin_disease_count = org_disease_count + new_disease_count;
        fin_disease_idx   = org_disease_idx;
    }
    else if (org_disease_idx == 0xFF) {
        fin_disease_count = new_disease_count;
        fin_disease_idx   = new_disease_idx;
    }
    else if (new_disease_idx == 0xFF) {
        fin_disease_count = org_disease_count;
        fin_disease_idx   = org_disease_idx;
    }
    else {
        float org_disease_strength = org_disease_count * this->diseases[org_disease_idx].strength;
        float new_disease_strength = new_disease_count * this->diseases[new_disease_idx].strength;

        if (org_disease_strength <= new_disease_strength && org_disease_count < 0) {
            fin_disease_count = -org_disease_count;
            fin_disease_idx = new_disease_idx;
        }
        else if (org_disease_strength <= new_disease_strength) {
            fin_disease_count = org_disease_count;
            fin_disease_idx = org_disease_idx;
        }
        else {
            int newSubOrg_count = new_disease_count - (int)(org_disease_strength / new_disease_strength * org_disease_count);
            if (newSubOrg_count < 0) {
                fin_disease_count = -newSubOrg_count;
                fin_disease_idx = org_disease_idx;
            }
            else if (newSubOrg_count >= 0) {
                fin_disease_count = newSubOrg_count;
                fin_disease_idx = new_disease_idx;
            }
        }
    }

    cells->diseaseCount[cell_idx] = fin_disease_count;
    cells->diseaseIdx  [cell_idx] = fin_disease_idx;
    if (org_disease_count <= 0) {
        cells->diseaseIdx                   [cell_idx] = -1;
        cells->diseaseCount                 [cell_idx] = 0;
        cells->diseaseInfestationTickCount  [cell_idx] = 0;
        cells->diseaseGrowthAccumulatedError[cell_idx] = 0;
        return;
    }
    if (org_disease_idx != fin_disease_idx) {
        cells->diseaseInfestationTickCount[cell_idx] = 0;
    }
}

float Disease::GetDiffusionScale(const SimData* simData, int cell_idx) 
{
    uint16_t elem_idx = simData->cells->elementIdx[cell_idx];
    DiseaseInfo::ElemGrowthInfo* info = &this->diseases[simData->cells->diseaseIdx[cell_idx]].elemGrowthInfo;
    if (simData->cells->diseaseCount[cell_idx] >= info->minDiffusionCount[elem_idx])
        if (simData->cells->diseaseInfestationTickCount[cell_idx] >= info->minDiffusionInfestationTickCount[elem_idx])
            if (elem_idx < info->diffusionScale.size())
                return info->diffusionScale[elem_idx];

    return 0.0f;
}

uint8_t Disease::GetDiseaseIndex(uint32_t disease_hash)
{
    for (uint8_t idx = 0; idx < this->diseases.size() && idx < 0xFF; idx++)
        if (this->diseases[idx].hashID == disease_hash)
            return idx;
    return 0xFF;
}

void Disease::PostProcess(SimData* simData, SimEvents* simEvents, int x_start, const int x_end, int y_start, int y_end)
{
    for (int y_pos = y_start; y_pos < y_end; y_pos++) {
        uint64_t cell     = y_pos * simData->width + x_start;
        uint64_t cell_end = y_pos * simData->width + x_end;
        if (cell >= cell_end)
            continue;

        for (; cell < cell_end; cell++) {
            //No disease, skip
            if (simData->updatedCells->diseaseIdx[cell] == 0xFF)
                continue;

            //Influence by temperature
            DiseaseInfo* info = &this->diseases[simData->updatedCells->diseaseIdx[cell]];
            int level = getTempLevel(simData->updatedCells->temperature[cell], &info->temperatureRange);
            int lev_H = std::min(level    , 3);
            int lev_L = std::max(level - 1, 0);

            //Calculate change_Rate by temperature
            float changeRate_Temp = 0;
            float halfLife_L = getTempValue(lev_L, &info->temperatureHalfLives);
            float halfLife_H = getTempValue(lev_H, &info->temperatureHalfLives);
            if (level == 2 || halfLife_L == INFINITY || halfLife_H == INFINITY) {
                changeRate_Temp = 1;
            }
            else {
                float temp_Ratio = 0;
                float temp_L     = getTempValue(lev_L, &info->temperatureRange);
                float temp_H     = getTempValue(lev_H, &info->temperatureRange);
                if (temp_H > temp_L)
                    temp_Ratio   = (simData->updatedCells->temperature[cell] - temp_L) / (temp_H - temp_L);

                changeRate_Temp = halfLife_L + temp_Ratio * (halfLife_H - halfLife_L);
                if (changeRate_Temp == INFINITY) changeRate_Temp = 1;
                else                             changeRate_Temp = powf(2.0f, -0.2f / changeRate_Temp);
            }

            int diseaseCount = simData->updatedCells->diseaseCount[cell];
            //Calculate change count by disease population
            float changeCount_Popu = 0;
            DiseaseInfo::ElemGrowthInfo* elem_info = &info->elemGrowthInfo;
            uint16_t elem_idx = simData->updatedCells->elementIdx[cell];
            if (diseaseCount >= simData->updatedCells->mass[cell] * elem_info->minCountPerKG[elem_idx]) {
                //Judge crowded
                float halfLife_Popu = elem_info->populationHalfLife[elem_idx];
                if (diseaseCount > simData->updatedCells->mass[cell] * elem_info->maxCountPerKG[elem_idx])
                    halfLife_Popu   = elem_info->overPopulationHalfLife[elem_idx];

                float changeRate_Popu = 0;
                if      (halfLife_Popu == 0.0)      changeRate_Popu = 0.0;
                else if (halfLife_Popu == INFINITY) changeRate_Popu = 1;
                else                                changeRate_Popu = powf(2.0f, -0.2f / halfLife_Popu);

                changeCount_Popu = (changeRate_Popu - 1.0f) * diseaseCount;
            }
            else {
                changeCount_Popu = elem_info->underPopulationDeathRate[elem_idx] * -0.2f;
            }

            //Apply count change
            float changeCount = diseaseCount * changeRate_Temp + simData->updatedCells->diseaseGrowthAccumulatedError[cell] - diseaseCount + changeCount_Popu;
            //Apply count change by Radiation
            if (simData->radiationEnabled)           
                changeCount -= (simData->updatedCells->radiation[cell] * info->radiationKillRate);
            //Record Integer Round-off error
            simData->updatedCells->diseaseGrowthAccumulatedError[cell] = changeCount - (float)(int)changeCount;
            simData->updatedCells->diseaseCount[cell] += (int)changeCount;
        }
    }

    for (int y_pos = y_start; y_pos < y_end; y_pos++) {
        uint64_t cell     = y_pos * simData->width + x_start;
        uint64_t cell_end = y_pos * simData->width + x_end;
        if (cell >= cell_end)
            continue;

        for (; cell < cell_end; cell++) {
            if (simData->updatedCells->diseaseCount[cell] < 0) {
                simData->updatedCells->diseaseIdx[cell]                    = -1;
                simData->updatedCells->diseaseCount[cell]                  = 0;
                simData->updatedCells->diseaseInfestationTickCount[cell]   = 0;
                simData->updatedCells->diseaseGrowthAccumulatedError[cell] = 0;
            }
            if (simData->updatedCells->diseaseInfestationTickCount[cell] < 0xFD)
                simData->updatedCells->diseaseInfestationTickCount[cell]++;
        }
    }
}

void Disease::UpdateCells(float dt, SimData* simData, SimEvents* simEvents, int cell, int ncell)
{
    if (((gElementPostProcessData[simData->cells->elementIdx[cell]].state ^ gElementPostProcessData[simData->cells->elementIdx[ncell]].state) & 3) != 0)
        return;

    int     count_C = simData->cells->diseaseCount[cell];
    int     count_N = simData->cells->diseaseCount[ncell];
    uint8_t idx_C   = simData->cells->diseaseIdx  [cell];
    uint8_t idx_N   = simData->cells->diseaseIdx  [ncell];

    uint8_t idx_fin = 0xFF;
    int     delta_C = 0;
    if (idx_C == 0xFF && idx_N == 0xFF) {
        // Both have no disease
        return;
    }
    else if (idx_C == idx_N) { 
        // Same disease
        int max_cell = count_N > count_C ? ncell : cell;
        delta_C = (int)((count_N - count_C) * 0.125f * this->GetDiffusionScale(simData, max_cell));
        idx_fin = idx_C;
    }
    else if (idx_C == 0xFF) { 
        // Only ncell has disease
        delta_C = (int)(count_N * 0.125f * this->GetDiffusionScale(simData, ncell));
        idx_fin = idx_N;
    }
    else if (idx_N == 0xFF) { 
        // Only cell has disease
        delta_C = -(int)(count_C * 0.125f * this->GetDiffusionScale(simData, cell));
        idx_fin = idx_C;
    }
    else if (count_C * this->diseases[idx_C].strength <= count_N * this->diseases[idx_N].strength) { 
        // Different diseases. Ncell is stronger
        delta_C = (int)(count_N * 0.125f * this->GetDiffusionScale(simData, ncell));
        idx_fin = idx_N;
    }
    else { 
        // Different diseases. Cell is stronger
        delta_C = -(int)(count_C * 0.125f * this->GetDiffusionScale(simData, cell));
        idx_fin = idx_C;
    }

    if ( delta_C) this->AddDiseaseToCell(simData->updatedCells.get(), cell,  idx_fin,  delta_C);
    if (-delta_C) this->AddDiseaseToCell(simData->updatedCells.get(), ncell, idx_fin, -delta_C);
}