#include "ClassSim.h"
#include "ClassDisease.h"

#define _USE_MATH_DEFINES
#include "math.h"

//----- (Local Struct) -----------------------------------------------
//Separated from RadiationAbsorptionAlongLine
struct RadiationAbsorption {
    SimData* simData;
    float    absorption;

    void doAbsorb(int locX, int locY) {
        int      cell = this->simData->width * locY + locX;
        uint16_t elem = simData->updatedCells->elementIdx[cell];
        float    factor = gElementRadiationData[elem].factor;
        float    absorb = factor * simData->RADIATION_CONSTRUCTED_FACTOR;
        // Not ConstructedTile
        if (!(simData->updatedCells->properties[cell] & 0x80)) {
            absorb = factor * (simData->RADIATION_BASE_WEIGHT +
                simData->updatedCells->mass[cell] / simData->RADIATION_MAX_MASS * simData->RADIATION_DENSITY_WEIGHT);
        }
        this->absorption *= 1 - CLAMP_F(absorb, 1, 0);
    }
};
//----- (New add) ----------------------------------------------------
//Flood__lambda_9390f724927d7d28953bae36843a8903_
void FindCellbyDepth(std::vector<uint16_t> elemIdx, int width, uint16_t src_x, uint16_t src_y, short remain_depth, 
    bool stop_at_solid, bool stop_at_liquid, bool stop_at_gas, std::unordered_map<int, short>* cellDic) /* cell, depth */
{
    int cell = src_y * width + src_x;
    uint8_t state = gElementPostProcessData[elemIdx[cell]].state & 3;

    if (state == 3 && stop_at_solid ) return;
    if (state == 2 && stop_at_liquid) return;
    if (state == 1 && stop_at_gas   ) return;

    if (cellDic->find(cell) != cellDic->end()) {
        //cell already in cellDic and has lower depth
        if (remain_depth <= (*cellDic)[cell]) return;
        //cell already in cellDic but has bigger depth
        (*cellDic)[cell] = remain_depth;
    }
    else {
        cellDic->insert(std::pair<int, short>(cell, remain_depth));
    }
    if (remain_depth > 1) {
        FindCellbyDepth(elemIdx, width, src_x, src_y - 1, remain_depth - 1, stop_at_solid, stop_at_liquid, stop_at_gas, cellDic);
        FindCellbyDepth(elemIdx, width, src_x - 1, src_y, remain_depth - 1, stop_at_solid, stop_at_liquid, stop_at_gas, cellDic);
        FindCellbyDepth(elemIdx, width, src_x, src_y + 1, remain_depth - 1, stop_at_solid, stop_at_liquid, stop_at_gas, cellDic);
        FindCellbyDepth(elemIdx, width, src_x + 1, src_y, remain_depth - 1, stop_at_solid, stop_at_liquid, stop_at_gas, cellDic);
    }
}

bool cmpDepth(std::pair<int, short> a, std::pair<int, short> b)
{
    if (a.second != b.second)
        return a.second < b.second;
    return a.first < b.first;
}

//----- (Decompiled) -------------------------------------------------
void GetReachableCells(SimData* simData, uint16_t x, uint16_t y, int max_depth, std::vector<int>* visited, 
    std::queue<FloodFillInfo>* next, std::vector<int>* reachable_cells, bool stop_at_solids)
{
    visited->clear();
    reachable_cells->clear();

    //Flood__lambda_9390f724927d7d28953bae36843a8903_(simData, x, y, max_depth, 
    // stop_at_solids, 0, 0, visited, next, reachable_cells);
    std::unordered_map<int, short> cellDic;
    for (; !next->empty(); next->pop()) {
        int cell = next->front().x + next->front().y * simData->width;
        cellDic[cell] = max_depth - next->front().depth;
    }
    FindCellbyDepth(simData->updatedCells->elementIdx, simData->width, x, y, max_depth, stop_at_solids, false, false, &cellDic);
    // LOGGER_PRINT2("%s-Pos[%d,%d], radius %d, count %llu\n", __func__, x, y, max_depth, cellDic.size());
#if 0 // unsorted
    for (const auto& pair : cellDic) {
#else // sorted
    std::vector<std::pair<int, short>> cellSorted(cellDic.begin(), cellDic.end());
    std::sort(cellSorted.begin(), cellSorted.end(), cmpDepth);

    for (const auto& pair : cellSorted) {
#endif
        visited->push_back(pair.first);
        reachable_cells->push_back(pair.first);
    }
}

bool AnyInputCellHasState(const SimData* simData, const std::vector<int>* valid_cells,
    uint16_t* desired_cell_elem_idx, ElementState::State state, int offset_base)
{
    for (int i = 0; i < valid_cells->size(); i++) {
        uint64_t offset   = (offset_base + i) % valid_cells->size();
        uint16_t elem_idx = simData->updatedCells->elementIdx[valid_cells->at(offset)];
        if ((state & gElements[elem_idx].state & 3) == state) {
            *desired_cell_elem_idx = elem_idx;
            return true;
        }
    }
    return false;
}

void InitializeFromMessageData(BuildingHeatExchangeData* result, BuildingHeatExchangeMessageData* data)
{
    int celllCount = (data->maxY - data->minY) * (data->maxX - data->minX);
    if (celllCount <= 0) ASSERT_TEXT("Size zero building passed into BuildingHeatExchange");
    ASSERT_TEMP(1, data->temperature);

    result->simMin.x = data->minX + 1;
    result->simMin.y = data->minY + 1;
    result->simMax.x = data->maxX + 1;
    result->simMax.y = data->maxY + 1;

    float HC = gElementTemperatureData[data->elemIdx].specificHeatCapacity * data->mass;
    float TC = gElementTemperatureData[data->elemIdx].thermalConductivity  * data->thermalConductivity;

    result->temperature         = data->temperature;
    result->overheatTemperature = data->overheatTemperature;
    result->operatingKilowatts  = data->operatingKilowatts;
    result->totalHeatCapacity   = HC;
    result->perCellHeatCapacity = HC / celllCount;
    result->thermalConductivity = TC;
    result->highTemp            = gElementTemperatureData[data->elemIdx].highTemp;
}

float ConductTemperature(float dt, float scale, float temp1, float HC1, float TC1, float insul1, 
    float temp2, float HC2, float TC2, float insul2, float* tempFin1, float* tempFin2)
{
    float heatTrans = fabsf(temp2 - temp1) * MIN_F(TC1 * insul1, TC2 * insul2) * scale * dt * 0.001f;
    heatTrans = MIN_F(heatTrans, fabsf(temp2 * HC2 - temp1 * HC1));
    if (temp2 <= temp1) heatTrans = -heatTrans;

    float tempMax = MAX_F(temp1, temp2);
    float tempMin = MIN_F(temp1, temp2);
    *tempFin1 = CLAMP_F(temp1 + heatTrans / HC1, tempMax, tempMin);
    *tempFin2 = CLAMP_F(temp2 - heatTrans / HC2, tempMax, tempMin);

    // Temperature inverse
    if ((*tempFin2 - *tempFin1) * (temp2 - temp1) < 0)    {
        *tempFin1 = (HC1 * temp1 + HC2 * temp2) / (HC1 + HC2);
        *tempFin2 = *tempFin1;
    }
    heatTrans = MIN_F(fabsf(*tempFin1 - temp1) * HC1, fabsf(*tempFin2 - temp2) * HC2);
    if (temp2 < temp1) heatTrans = -heatTrans;

    *tempFin1 = CLAMP_F(temp1 + heatTrans / HC1, tempMax, tempMin);
    *tempFin2 = CLAMP_F(temp2 - heatTrans / HC2, tempMax, tempMin);
    return -heatTrans;
}

//float ExchangeHeatEnergyWithWorld(float dt, SimData* simData, Handle* sim_handle, CellAccessor* sim_cell, float scale_factor, ElementChunkData* item)
float ExchangeHeatEnergyWithWorld(float dt, CellAccessor* sim_cell, float scale_factor, ElementChunkData* item)
{
    float massCell = sim_cell->cells->mass[sim_cell->cellIdx];
    if (massCell <= 0) return 0;
    ASSERT_TEMP(massCell, item->temperature);
    ASSERT_TEMP(massCell, sim_cell->cells->temperature[sim_cell->cellIdx]);

    //TemperatureInsulated
    uint16_t elem = sim_cell->cells->elementIdx[sim_cell->cellIdx];
    if (gElementTemperatureData[elem].state & 0x10) return 0;

    float tempFinCell, tempFinItem;
    float scale     = scale_factor * item->maxEnergyTransferScaleFactor;
    float HcCell    = massCell * gElementTemperatureData[elem].specificHeatCapacity;
    float TcCell    = gElementTemperatureData[elem].thermalConductivity;
    float insulCell = sim_cell->cells->insulation[sim_cell->cellIdx];

    float heatTrans = ConductTemperature(dt, scale, sim_cell->temperature(), HcCell, TcCell, insulCell * insulCell / 255 / 255,
        item->temperature, item->heatCapacity, item->thermalConductivity, 1, &tempFinCell, &tempFinItem);
    //LOGGER_PRINT2("%s cell %f->%f, item %f->%f, heat %f\n", __func__, sim_cell->temperature(), tempFinCell, item->temperature, tempFinItem, heatTrans);
    sim_cell->cells->temperature[sim_cell->cellIdx] = tempFinCell;
    item->temperature                               = tempFinItem;
    return heatTrans;
}

float RadiationAbsorptionAlongLine(SimData* simData, int startX, int startY, int endX, int endY)
{
    RadiationAbsorption data = { simData, 1 };
    int deltaX = std::abs(startX - endX);
    int deltaY = std::abs(startY - endY);
    // Rewritten
    int step = std::max(deltaX, deltaY);
    for (int idx = 0; idx <= step; idx++) {
        int locX = step ? (startX + (endX - startX) * idx / step) : startX;
        int locY = step ? (startY + (endY - startY) * idx / step) : startY;
        data.doAbsorb(locX, locY);
        if (data.absorption < 0)
            break;
    }
    return CLAMP_F(data.absorption, 1, 0);
}

bool inRadialRange(int originX, int originY, int x, int y, float emitAngle, float emitDirection)
{
    if (emitAngle == 360)             return true;
    if (originX == x && originY == y) return true;

    double angN = fmod(emitDirection - emitAngle * 0.5 + 360, 360);
    double angP = fmod(emitDirection + emitAngle * 0.5 + 360, 360);
    double tan  = atan2((float)(y - originY), (float)(x - originX));
    double ang  = fmod(tan * 180.0 / M_PI + 360.0, 360.0);
    if (angP <= angN && angN <= ang) return true;
    if (angP >  angN && angN >  ang) return false;
    return ang <= angP;
}

//void setAttractor(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads, float* out_rads)
float setAttractor(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads)
{
    if (cellX < 1)                 return 0;
    if (simData->width  <= cellX)  return 0;
    if (cellY < 1)                 return 0;
    if (simData->height <= cellY)  return 0;
    if (cellX < region->minimum.x) return 0;
    if (cellY < region->minimum.y) return 0;
    if (cellX > region->maximum.x) return 0;
    if (cellY > region->maximum.y) return 0;
    if (!inRadialRange(centerX, centerY, cellX, cellY, item->emitAngle, item->emitDirection)) return 0;
    if (centerX == cellX && centerY == cellY) return 0;

    int cell = cellX + cellY * simData->width;
    if (simData->updatedCells->radiation[cell] < 0.01f) return 0;
    rads = MIN_F(rads, simData->updatedCells->radiation[cell]);
    simData->updatedCells->radiation[cell] -= rads;

    return rads;
}

float setPulsing(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads)
{
    if (rads <= 0)                 return 0;
    if (cellX < 1)                 return 0;
    if (simData->width  <= cellX)  return 0;
    if (cellY < 1)                 return 0;
    if (simData->height <= cellY)  return 0;
    if (cellX < region->minimum.x) return 0;
    if (cellY < region->minimum.y) return 0;
    if (cellX > region->maximum.x) return 0;
    if (cellY > region->maximum.y) return 0;
    if (!inRadialRange(centerX, centerY, cellX, cellY, item->emitAngle, item->emitDirection)) return 0;
    if (centerX == cellX && centerY == cellY) return 0;

    int   cell   = cellX + cellY * simData->width;
    float absorb = RadiationAbsorptionAlongLine(simData, centerX, centerY, cellX, cellY);
    simData->updatedCells->radiation[cell] += absorb * rads;
    return 0;
}

float setSimplePulse(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads)
{
    if (cellX < 1)                 return 0;
    if (simData->width  <= cellX)  return 0;
    if (cellY < 1)                 return 0;
    if (simData->height <= cellY)  return 0;
    if (cellX < region->minimum.x) return 0;
    if (cellY < region->minimum.y) return 0;
    if (cellX > region->maximum.x) return 0;
    if (cellY > region->maximum.y) return 0;
    if (!inRadialRange(centerX, centerY, cellX, cellY, item->emitAngle, item->emitDirection)) return 0;

    int cell = cellX + cellY * simData->width;
    simData->updatedCells->radiation[cell] += rads;
    return 0;
}

float SetPixel4(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int deltaX, int deltaY, float rads,
    float(*func)(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads))
{
    float out_rads = 0;
    out_rads += func(simData, region, item, centerX, centerY, centerX + deltaX, centerY + deltaY, rads);
    out_rads += func(simData, region, item, centerX, centerY, centerX - deltaX, centerY + deltaY, rads);
    out_rads += func(simData, region, item, centerX, centerY, centerX + deltaX, centerY - deltaY, rads);
    out_rads += func(simData, region, item, centerX, centerY, centerX - deltaX, centerY - deltaY, rads);
    return out_rads;
}

float SetCircleAA(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int radiusX, int radiusY, float rads, 
    float(*func)(SimData* simData, Region* region, RadiationEmitterData* item, int centerX, int centerY, int cellX, int cellY, float rads))
{
    if (radiusX <= 0 || radiusY <= 0)
        return -1;

    float radius  = std::sqrtf((float)(radiusX * radiusX + radiusY * radiusY));
    float out_rads = 0;

    float lengthX = std::roundf(radiusX * radiusX / radius);
    for (float deltaX = 0; deltaX < lengthX; deltaX++) {
        float deltaY   = std::sqrtf(1 - (deltaX / radiusX) * (deltaX / radiusX)) * radiusY;
        float radLocal = std::roundf(deltaY * rads);
        out_rads += SetPixel4(simData, region, item, centerX, centerY, (int)deltaX, (int)deltaY    ,        radLocal, func);
        out_rads += SetPixel4(simData, region, item, centerX, centerY, (int)deltaX, (int)deltaY - 1, rads - radLocal, func);
    }
    float lengthY = std::roundf(radiusY * radiusY / radius);
    for (float deltaY = 0; deltaY < lengthY; deltaY++) {
        float deltaX   = std::sqrtf(1 - (deltaY / radiusY) * (deltaY / radiusY)) * radiusX;
        float radLocal = std::roundf(deltaX * rads);
        out_rads += SetPixel4(simData, region, item, centerX, centerY, (int)deltaX    , (int)deltaY,        radLocal, func);
        out_rads += SetPixel4(simData, region, item, centerX, centerY, (int)deltaX - 1, (int)deltaY, rads - radLocal, func);
    }

    return out_rads;
}

void tickConstant(SimData* simData, Region* region, RadiationEmitterData* item, float dt)
{
    int cenX = item->cell % simData->width;
    int cenY = item->cell / simData->width;
#ifdef __DEBUGED__
    for (int locY = cenY - item->radiusY; locY <= cenY + item->radiusY; locY++) {
#else
    for (int locY = cenY - item->radiusY; locY <= cenY - item->radiusY + 2 * item->radiusX; locY++) {
#endif
        if (locY <  1)                 continue;
        if (locY >= simData->height)   continue;
        if (locY <  region->minimum.y) continue;
        if (locY >  region->maximum.y) continue;

        for (int locX = cenX - item->radiusX; locX <= cenX + item->radiusX; locX++) {
            if (locX <  1)                 continue;
            if (locX >= simData->width)    continue;
            if (locX <  region->minimum.x) continue;
            if (locX >  region->maximum.x) continue;

            if (!inRadialRange(cenX, cenY, locX, locY, item->emitAngle, item->emitDirection)) continue;
            float ratioX = (float)(locX - cenX) / item->radiusX;
            float ratioY = (float)(locY - cenY) / item->radiusY;
            if (ratioY * ratioY + ratioX * ratioX > 1) continue;

            float loss   = 1 - fabsf(ratioX) - fabsf(ratioY) + fabsf(ratioX * ratioY);
            float lossEx = 0;
            float absorb = RadiationAbsorptionAlongLine(simData, cenX, cenY, locX, locY) * loss * item->emitRads / simData->RADIATION_LINGER_RATE;
            if (loss < 0.25)
                lossEx = (float)(RAND_FLOAT_01(simData->randomSeed) * 0.25 - absorb * 0.125);
            
            simData->updatedCells->radiation[locX + simData->width * locY] += lossEx + absorb;
        }
    }
}

void tickPulsing(SimData* simData, Region* region, RadiationEmitterData* item, float dt, int stepCount)
{
    float stepRatio = (float)(item->emitStep + 1) / stepCount;
    int   cenX      = item->cell % simData->width;
    int   cenY      = item->cell / simData->width;
    int   xRatio    = (int)(stepRatio * item->radiusX);
    int   yRatio    = (int)(stepRatio * item->radiusY);

    float emitRads = item->emitRads;
    if (item->emitType == RadiationEmitterType::PulsingAveraged)
        emitRads /= MAX_F((float)((xRatio - 1) * (yRatio = 1)), 1);

    if (!item->emitStep)
        simData->updatedCells->radiation[item->cell] += emitRads;
    else if (xRatio > 1 || cenY > 1)
        SetCircleAA(simData, region, item, cenX, cenY, xRatio, yRatio, emitRads, setPulsing);

    return;
}

void tickSimplePulse(SimData* simData, Region* region, RadiationEmitterData* item, float dt, int stepCount)
{
    float stepRatio = (float)(item->emitStep + 1) / stepCount;
    int   cenX      = item->cell % simData->width;
    int   cenY      = item->cell / simData->width;
    int   xRatio    = (int)MAX_F(stepRatio * item->radiusX, 1);
    int   yRatio    = (int)MAX_F(stepRatio * item->radiusY, 1);

    SetCircleAA(simData, region, item, cenX, cenY, xRatio, yRatio, item->emitRads, setSimplePulse);
    return;
}

void tickAttractor(SimData* simData, Region* region, RadiationEmitterData* item, float dt, int stepCount)
{
    float stepRatio = (float)(item->emitStep + 1) / stepCount;
    int   cenX      = item->cell % simData->width;
    int   cenY      = item->cell / simData->width;
    int   xRatio    = (int)MAX_F(stepRatio * item->radiusX, 1);
    int   yRatio    = (int)MAX_F(stepRatio * item->radiusY, 1);

    float rads = SetCircleAA(simData, region, item, cenX, cenY, xRatio, yRatio, item->emitRads, setAttractor);
    if (rads > 0)
        simData->updatedCells->radiation[item->cell] += rads;
    return;
}

//----- (CompactedVector) --------------------------------------------
template <typename T>
Handle CompactedVector<T>::Free(int handle)
{
    int     handleVal = handle & 0xFFFFFF;
    uint8_t handleVer = (handle >> 24) + 1;
    if (handleVal >= 0x100000)
        return handle;

    this->handles.versions[handleVal] = handleVer;
    this->handles.freeHandles.push_back(handleVal | (handleVer << 24));
    int handleItm = this->handles.items[handleVal];
    this->handles.items[handleVal] = 0;
    // Replace the target handle with the last handle
    int posEnd = (int)(this->data.size() - 1);
    if (handleItm < posEnd) {
        this->data[handleItm] = this->data[posEnd];
        int handleVal = this->dataHandleIndices[posEnd] & 0xFFFFFF;
        this->handles.items[handleVal] = handleItm;
        this->dataHandleIndices[handleItm] = this->dataHandleIndices[posEnd];
    }
    this->data.pop_back();
    this->dataHandleIndices.pop_back();
    return 0xFFFFFF;
}

template <typename T>
Handle CompactedVector<T>::AddData(T* new_data)
{
    Handle result = -1;
    if (this->handles.freeHandles.size()) {
        result = this->handles.freeHandles.back();
        this->handles.freeHandles.pop_back();
        this->handles.items[result & 0xFFFFFF] = (Handle)this->data.size();
    }
    else {
        this->handles.versions.push_back(0);
        result = this->handles.items.size() & 0xFFFFFF;
        this->handles.items.push_back((Handle)this->data.size());
    }
    this->dataHandleIndices.push_back(result);
    this->data.push_back(*new_data);

    return result;
}

template <typename T>
void CompactedVector<T>::SetData(int handle, T* new_data)
{
    if (!HANDLE_AVAILABLE(handle, this->handles.versions)) {
        ASSERT_TEXT("Illegal Handle or Version");
    }
    int handleItm = this->handles.items[handle & 0xFFFFFF];
    memcpy(&this->data[handleItm], new_data, sizeof(T));
}

template <typename T>
T* CompactedVector<T>::GetData(int handle)
{
    if (!HANDLE_AVAILABLE(handle, this->handles.versions)) {
        ASSERT_TEXT("Illegal Handle or Version");
    }
    return &this->data[this->handles.items[handle & 0xFFFFFF]];
}

//----- (SimComponent class) -----------------------------------------
Handle BuildingHeatExchange::Register(SimData* simData, AddBuildingHeatExchangeMessage* msg)
{
    BuildingHeatExchangeData new_data;
    InitializeFromMessageData(&new_data, &msg->data);

    return this->AddData(&new_data);
}

void BuildingHeatExchange::Modify(SimData* simData, ModifyBuildingHeatExchangeMessage* msg)
{
    BuildingHeatExchangeData new_data;
    InitializeFromMessageData(&new_data, &msg->data);

    BuildingHeatExchangeData* Data = this->GetData(msg->handle);
    if (Data->overheatTemperature <= Data->temperature && new_data.temperature < Data->overheatTemperature) {
        MeltedInfo melted_info = { .handle = msg->handle };
        simData->simEvents->buildingNoLongerOverheatedInfo.push_back(melted_info);
    }

    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) {
        ASSERT_TEXT("Illegal Handle or Version");
    }
    int handleItm = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[handleItm] = new_data;
}

void BuildingHeatExchange::Unregister(int handle)
{
    BuildingHeatExchangeData* data = this->GetData(handle);
    data->simMin.x = 0;
    data->simMin.y = 0;
    data->simMax.x = 0;
    data->simMax.y = 0;
    data->temperature         = -1;
    data->perCellHeatCapacity = -1;
    data->thermalConductivity = -1;

    this->temperatureInfo[handle & 0xFFFFFF].temperature = -1;
    this->Free(handle);
}

void BuildingHeatExchange::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "BuildingHeatExchange %s-%llu\n", __func__, this->data.size());
    this->temperatureInfo.resize(this->handles.items.size());
    for (int i = 0; i < this->data.size(); i++) {
        if (this->data[i].simMin.x < region->minimum.x) continue;
        if (this->data[i].simMin.y < region->minimum.y) continue;
        if (this->data[i].simMax.x > region->maximum.x) continue;
        if (this->data[i].simMax.y > region->maximum.y) continue;

        if (this->data[i].perCellHeatCapacity > 0 && this->data[i].temperature > 0) {
            float tempBulid      = this->data[i].temperature;
            float tempBuildMax   = tempBulid;
            float tempBuildMin   = tempBulid;
            float HeatTransBuild = 0;

            for (int locX = this->data[i].simMin.x; locX < this->data[i].simMax.x; locX++) {
                for (int locY = this->data[i].simMin.y; locY < this->data[i].simMax.y; locY++) {
                    int cell = locX + locY * simData->width;
                    if (simData->updatedCells->mass[cell] <= 0) continue;

                    ElementTemperatureData* elem = &gElementTemperatureData[simData->updatedCells->elementIdx[cell]];
                    float HcCell    = simData->updatedCells->mass[cell] * elem->specificHeatCapacity;
                    float tempCell  = simData->updatedCells->temperature[cell];
                    float HcLocal   = tempCell >= tempBulid ? HcCell : this->data[i].perCellHeatCapacity;
                    float HeatTrans = HcLocal * (tempCell - tempBulid) * elem->thermalConductivity * this->data[i].thermalConductivity
                                              * simData->updatedCells->insulation[cell] * simData->updatedCells->insulation[cell] 
                                              * 0.0000000768935f * dt //  0.0000000768935 = 1 / 255 / 255 / 200
                                              * simData->debugProperties.buildingTemperatureScale; // 100

                    float tempLocalMax = MAX_F(tempBulid,    tempCell);
                    float tempLocalMin = MIN_F(tempBulid,    tempCell);
                    tempBuildMax       = MAX_F(tempLocalMax, tempBuildMax);
                    tempBuildMin       = MIN_F(tempLocalMin, tempBuildMin);
                    float tempCellFin  = CLAMP_F(tempCell - (HeatTrans / HcCell), tempLocalMax, tempLocalMin);

                    // Temperature inverse
                    float tempBulildLocal = CLAMP_F(tempBulid + HeatTrans / this->data[i].perCellHeatCapacity, tempLocalMax, tempLocalMin);
                    if ((tempBulildLocal - tempCellFin) * (tempBulid - tempCell) < 0) {
                        tempCellFin = (this->data[i].perCellHeatCapacity * tempBulid + HcCell * tempCell) / (this->data[i].perCellHeatCapacity + HcCell);
                        tempCellFin = CLAMP_F(tempCellFin, tempLocalMax, tempLocalMin);
                        HeatTrans   = (tempCellFin - tempBulid) * this->data[i].perCellHeatCapacity;
                    }
                    ASSERT_TEMP(1, tempCellFin);
                    if (std::isnan(tempCellFin)) {
                        LOGGER_PRINT("NAN temperature. Cell %d, Capacity %f/%f\n", cell, this->data[i].perCellHeatCapacity, HcCell);
                        if (HcCell == 0)
                            LOGGER_PRINT("Cell Elem %s. mass %f, SHC %f\n", gElementNames[simData->updatedCells->elementIdx[cell]].c_str(), simData->updatedCells->mass[cell], elem->specificHeatCapacity);
                    }
                    simData->updatedCells->temperature[cell] = tempCellFin;
                    DoStateTransition(simData, simData->simEvents.get(), cell, elem);
                    HeatTransBuild += HeatTrans;
                }
            }

            int   areaBuild    = (this->data[i].simMax.y - this->data[i].simMin.y) * (this->data[i].simMax.x - this->data[i].simMin.x);
            float HcBuild      = areaBuild * this->data[i].perCellHeatCapacity;
            float tempBuildFin = CLAMP_F(this->data[i].temperature + HeatTransBuild / HcBuild, tempBuildMax, tempBuildMin);
            tempBuildFin      += this->data[i].operatingKilowatts * dt / HcBuild;
            ASSERT_TEMP(1, tempBuildFin);

            if (tempBuildFin >= this->data[i].overheatTemperature) {
                MeltedInfo info = { .handle = this->dataHandleIndices[i] };
                simData->simEvents->buildingOverheatInfo.push_back(info);
            }
            else if (tempBulid >= this->data[i].overheatTemperature) {
                MeltedInfo info = { .handle = this->dataHandleIndices[i] };
                simData->simEvents->buildingNoLongerOverheatedInfo.push_back(info);
            }
            if (tempBuildFin >= this->data[i].highTemp) {
                MeltedInfo info = { .handle = this->dataHandleIndices[i] };
                simData->simEvents->buildingMeltedInfo.push_back(info);
            }
            this->data[i].temperature = tempBuildFin;
        }

        int handleVal = this->dataHandleIndices[i] & 0xFFFFFF;
        this->temperatureInfo[handleVal].handle      = this->dataHandleIndices[i];
        this->temperatureInfo[handleVal].temperature = this->data[i].temperature;
    }

    LOGGER_LEVEL(1, "BuildingHeatExchange %s-done\n", __func__);
}

void BuildingHeatExchange::UpdateDataListOnly(SimData* simData)
{
    this->temperatureInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
        this->temperatureInfo[handleVal].handle      = this->dataHandleIndices[idx];
        this->temperatureInfo[handleVal].temperature = this->data[idx].temperature;
    }
}

void BuildingHeatExchange::ChangeBuildingTemperature(int buildingHandle, float newTemperature)
{
    if (!HANDLE_AVAILABLE(buildingHandle, this->handles.versions)) {
        ASSERT_TEXT("Illegal Handle or Version");
    }
    int handleItm = this->handles.items[buildingHandle & 0xFFFFFF];
    int handleVal = this->dataHandleIndices[handleItm] & 0xFFFFFF;
    this->temperatureInfo[handleVal].temperature = newTemperature;
    this->temperatureInfo[handleVal].handle      = this->dataHandleIndices[handleItm];
}

Handle BuildingHeatExchange::GetSimHandleForBuildingHandle(int buildingHandle)
{
    if (!HANDLE_AVAILABLE(buildingHandle, this->handles.versions)) {
        ASSERT_TEXT("Illegal Handle or Version");
    }
    int handleItm = this->handles.items[buildingHandle & 0xFFFFFF];
    return this->dataHandleIndices[handleItm];
}

Handle BuildingToBuildingHeatExchange::Register(SimData* simData, RegisterBuildingToBuildingHeatExchangeMessage* msg)
{
    Handle handle = HANDLE_AVAILABLE(msg->heatExchange_handle, simData->buildingHeatExchange.handles.versions) ? msg->heatExchange_handle : -1;
    if (msg->callbackIdx != -1) {
        ComponentStateChangedMessage callback = { .callbackIdx = msg->callbackIdx, .simHandle = handle };
        simData->simEvents->componentStateChangedMessages.push_back(callback);
    }
    return handle;
}

void BuildingToBuildingHeatExchange::Unregister(int handle)
{
    BuildingToBuildingHeatExchangeData* data = this->GetData(handle);
    data->heatExchange_handle = 0;
    this->Free(handle);
}

void BuildingToBuildingHeatExchange::Add(SimData* simData, AddInContactBuildingToBuildingToBuildingHeatExchangeMessage* msg)
{
    BuildingToBuildingHeatExchangeData* data = this->GetData(msg->self_handle);
    BuildingToBuildingHeatExchangeData result = { 
        .heatExchange_handle    = data->heatExchange_handle , 
        .inContactBuildingsData = std::vector<InContactBuildingData>(data->inContactBuildingsData) 
    };
    InContactBuildingData info = { .inContactBuildingHandler = msg->buildingInContact, .cellsInContact = msg->cellsInContact };
    result.inContactBuildingsData.push_back(info);
    this->SetData(msg->self_handle, &result);
}

void BuildingToBuildingHeatExchange::Remove(SimData* simData, RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage* msg)
{
    BuildingToBuildingHeatExchangeData* data = this->GetData(msg->self_handle);
    BuildingToBuildingHeatExchangeData current_data = {
        .heatExchange_handle    = data->heatExchange_handle ,
        .inContactBuildingsData = std::vector<InContactBuildingData>(data->inContactBuildingsData)
    };

    for (auto it = current_data.inContactBuildingsData.begin(); it != current_data.inContactBuildingsData.end(); it++) {
        if (it->inContactBuildingHandler == msg->buildingInContact) {
            current_data.inContactBuildingsData.erase(it);
            this->SetData(msg->self_handle, &current_data);
            return;
        }
    }
}

void BuildingToBuildingHeatExchange::Update(float dt, SimData* simData, Region* region)
{
#define GET_BUILD_BY_HANDLE(_handle) \
    (simData->buildingHeatExchange.data[simData->buildingHeatExchange.handles.items[_handle & 0xFFFFFF]])

    LOGGER_LEVEL(1, "BuildingToBuildingHeatExchange %s-%llu\n", __func__, this->data.size());
    for (BuildingToBuildingHeatExchangeData data : this->data) {
        if (!HANDLE_AVAILABLE(data.heatExchange_handle, simData->buildingHeatExchange.handles.versions))
            continue;

        BuildingHeatExchangeData& build = GET_BUILD_BY_HANDLE(data.heatExchange_handle);
        if (build.simMin.x < region->minimum.x) continue;
        if (build.simMin.y < region->minimum.y) continue;
        if (build.simMax.x > region->maximum.x) continue;
        if (build.simMax.y > region->maximum.y) continue;

        float tempMin = build.temperature;
        float tempMax = build.temperature;
        for (InContactBuildingData dataCont : data.inContactBuildingsData) {
            if (!HANDLE_AVAILABLE(dataCont.inContactBuildingHandler, simData->buildingHeatExchange.handles.versions))
                continue;
            BuildingHeatExchangeData* dataExchange = simData->buildingHeatExchange.GetData(dataCont.inContactBuildingHandler);
            tempMin = MIN_F(dataExchange->temperature, tempMin);
            tempMax = MAX_F(dataExchange->temperature, tempMax);
        }

        for (InContactBuildingData dataCont : data.inContactBuildingsData) {
            if (!HANDLE_AVAILABLE(dataCont.inContactBuildingHandler, simData->buildingHeatExchange.handles.versions))
                continue;
            BuildingHeatExchangeData& buildCont = GET_BUILD_BY_HANDLE(dataCont.inContactBuildingHandler);
            if (buildCont.totalHeatCapacity <= 0) continue;

            float HcLocal   = buildCont.temperature >= build.temperature ? buildCont.totalHeatCapacity : build.totalHeatCapacity;
            float heatTrans = (buildCont.temperature - build.temperature) * buildCont.thermalConductivity * build.thermalConductivity
                * HcLocal * dt * simData->debugProperties.buildingToBuildingTemperatureScale * 0.005f;
            float tempCont  = CLAMP_F(buildCont.temperature - heatTrans / buildCont.totalHeatCapacity, tempMax, tempMin);
            float tempBuild = CLAMP_F(build    .temperature + heatTrans / build    .totalHeatCapacity, tempMax, tempMin);
            heatTrans = MIN_F(fabsf(tempCont - buildCont.temperature) * buildCont.totalHeatCapacity, fabsf(tempBuild - build.temperature) * build.totalHeatCapacity);
            if (buildCont.temperature < build.temperature)
                heatTrans = -heatTrans;
            tempCont  = buildCont.temperature - heatTrans / buildCont.totalHeatCapacity;
            tempBuild = build    .temperature + heatTrans / build    .totalHeatCapacity;

            // Temperature inverse
            if ((tempCont - tempBuild) * (buildCont.temperature - build.temperature) < 0) {
                tempCont  = (buildCont.temperature * buildCont.totalHeatCapacity + build.temperature * build.totalHeatCapacity) / (buildCont.totalHeatCapacity + build.totalHeatCapacity);
                tempCont  = CLAMP_F(tempCont, tempMax, tempMin);
                tempBuild = tempCont;
            }

            Handle handleCont = simData->buildingHeatExchange.GetSimHandleForBuildingHandle(dataCont.inContactBuildingHandler);
            if (tempCont >= buildCont.overheatTemperature) {
                MeltedInfo info = { .handle = handleCont };
                simData->simEvents->buildingOverheatInfo.push_back(info);
            }
            else if (buildCont.temperature >= buildCont.overheatTemperature) {
                MeltedInfo info = { .handle = handleCont };
                simData->simEvents->buildingNoLongerOverheatedInfo.push_back(info);
            }
            buildCont.temperature = tempCont;
            simData->buildingHeatExchange.ChangeBuildingTemperature(dataCont.inContactBuildingHandler, tempCont );
            simData->buildingHeatExchange.ChangeBuildingTemperature(data    .heatExchange_handle,      tempBuild);
        }
    }

    LOGGER_LEVEL(1, "BuildingToBuildingHeatExchange %s-done\n", __func__);
}

Handle ElementChunk::Register(SimData* simData, AddElementChunkMessage* msg)
{
    ElementTemperatureData* elemData = &gElementTemperatureData[msg->elementIdx];
    ElementChunkData new_data = {
        .cell                         = CELL_G2S(msg->gameCell, simData->width),
        .temperature                  = msg->temperature,
        .heatCapacity                 = msg->mass < 0.001 ? 0 : msg->mass * elemData->specificHeatCapacity,
        .thermalConductivity          = elemData->thermalConductivity,
        .maxEnergyTransferScaleFactor = msg->surfaceArea / msg->thickness,
        .groundTransferScale          = msg->groundTransferScale,
        .highTemp                     = elemData->highTemp,
        .lowTemp                      = elemData->lowTemp,
        .adjuster = {
            .temperature              = 0,
            .heatCapacity             = 0,
            .thermalConductivity      = 0
        }
    };

    return this->AddData(&new_data);
}

void ElementChunk::Unregister(int handle)
{
    this->Free(handle);
}

void ElementChunk::Modify(SimData* simData, ModifyElementChunkMessage* msg)
{
    ElementChunkData* Data = this->GetData(msg->handle);
    Data->temperature  = msg->temperature;
    Data->heatCapacity = msg->heatCapacity;
    ASSERT_TEMP(Data->heatCapacity, Data->temperature);
}

void ElementChunk::ModifyAdjuster(SimData* simData, ModifyElementChunkAdjusterMessage* msg)
{
    ElementChunkData::ExternalAdjuster* Data = &this->GetData(msg->handle)->adjuster;
    Data->temperature         = msg->temperature;
    Data->heatCapacity        = msg->heatCapacity;
    Data->thermalConductivity = msg->thermalConductivity;
}

void ElementChunk::ModifyEnergy(SimData* simData, ModifyElementChunkEnergyMessage* msg)
{
    ElementChunkData* Data = this->GetData(msg->handle);
    if (Data->heatCapacity > 0)
        Data->temperature = CLAMP_F(msg->deltaKJ / Data->heatCapacity + Data->temperature, SIM_MAX_TEMPERATURE, 0);
}

void ElementChunk::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "ElementChunk %s-%llu\n", __func__, this->data.size());
    this->chunkInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        float heatTrans = 0;
        int cell = this->data[idx].cell;
        if (cell % simData->width >= region->minimum.x &&
            cell / simData->width >= region->minimum.y &&
            cell % simData->width <= region->maximum.x &&
            cell / simData->width <= region->maximum.y &&
            this->data[idx].heatCapacity)
        {
            ElementChunkData::ExternalAdjuster& adjuster = this->data[idx].adjuster;
            if (adjuster.heatCapacity > 0) {
                float tempFinCell, tempFinItem;
                ConductTemperature(dt, 1, adjuster.temperature, adjuster.heatCapacity, adjuster.thermalConductivity, 1,
                    data[idx].temperature, this->data[idx].heatCapacity, adjuster.thermalConductivity, 1, &tempFinCell, &tempFinItem);
                data[idx].temperature = tempFinItem;
            }
            else {
                CellAccessor sim_cell (simData->updatedCells.get(), cell);
                CellAccessor sim_cellB(simData->updatedCells.get(), cell - simData->width);

                heatTrans      = ExchangeHeatEnergyWithWorld(dt, &sim_cell,  1,                             &this->data[idx]);
                if ((gElementStateData[simData->updatedCells->elementIdx[cell - simData->width]].state & 3) == 3)
                    heatTrans += ExchangeHeatEnergyWithWorld(dt, &sim_cellB, data[idx].groundTransferScale, &this->data[idx]);
            }
        }

        int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
        this->chunkInfo[handleVal].temperature = this->data[idx].temperature;
        this->chunkInfo[handleVal].deltaKJ    += heatTrans;

        if (this->data[idx].temperature > this->data[idx].highTemp + 3 || 
            this->data[idx].temperature < this->data[idx].lowTemp - 3)
        {
            MeltedInfo info = { .handle = this->dataHandleIndices[idx] };
            simData->simEvents->elementChunkMeltedInfo.push_back(info);
        }
    }

    LOGGER_LEVEL(1, "ElementChunk %s-done\n", __func__);
}

void ElementChunk::UpdateDataListOnly(SimData* simData)
{
    this->chunkInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
        this->chunkInfo[handleVal].temperature = this->data[idx].temperature;
        this->chunkInfo[handleVal].deltaKJ     = 0;

        ASSERT_TEMP(this->data[idx].heatCapacity, this->chunkInfo[handleVal].temperature);
    }
}

Handle ElementConsumer::Register(SimData* simData, AddElementConsumerMessage* msg)
{
    int cell = CELL_G2S(msg->gameCell, simData->width);
    ElementConsumerData new_data = {
        .consumptionRate = 0,
        .maxDepth        = msg->radius,
        .configuration   = msg->configuration,
        .elemIdx         = msg->elementIdx,
        .offsetIdx       = 0,
        .cell            = Vector2<uint16_t>(cell % simData->width, cell / simData->width)
    };

    return this->AddData(&new_data);
}

void ElementConsumer::Unregister(int handle)
{
    this->Free(handle);
}

void ElementConsumer::Modify(SimData* simData, ModifyElementConsumerMessage* msg)
{
    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) return;

    int item = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[item].cell.x = msg->gameCell % (simData->width - 2) + 1;
    this->data[item].cell.y = msg->gameCell / (simData->width - 2) + 1;
    this->data[item].consumptionRate = msg->consumptionRate;
}

void ElementConsumer::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "ElementConsumer %s. Total: %llu\n", __func__, this->data.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        ElementConsumerData& data = this->data[idx];
        //int cell = data.cell.y * simData->width + data.cell.x; //debug

        if (data.cell.x < region->minimum.x) continue;
        if (data.cell.y < region->minimum.y) continue;
        if (data.cell.x > region->maximum.x) continue;
        if (data.cell.y > region->maximum.y) continue;
        
        uint16_t elemRemove = data.elemIdx; // data.configuration == 0
        bool     findFlag   = true;         // data.configuration == 0
        if (data.configuration == 1) {
            GetReachableCells(simData, data.cell.x, data.cell.y, data.maxDepth, &this->visited, &this->next, &this->reachableCells, true);
            findFlag = AnyInputCellHasState(simData, &this->reachableCells, &elemRemove, ElementState::State::Liquid, data.offsetIdx);
        }
        else if (data.configuration == 2) {
            GetReachableCells(simData, data.cell.x, data.cell.y, data.maxDepth, &this->visited, &this->next, &this->reachableCells, true);
            findFlag = AnyInputCellHasState(simData, &this->reachableCells, &elemRemove, ElementState::State::Gas, data.offsetIdx);
        }

        this->reachableCells.clear();
        this->visited.clear();
        if (!findFlag)
            continue;

        ConsumedMassInfo info = {
            .simHandle      = this->dataHandleIndices[idx],
            .removedElemIdx = elemRemove,
            .diseaseIdx     = 0xFF,
            .mass           = 0,
            .temperature    = 0,
            .diseaseCount   = 0,
        };
        float massConsume = dt * data.consumptionRate;
        FloodRemoved(simData, &massConsume, elemRemove, &info, data.cell.x, data.cell.y, data.maxDepth, &this->visited, &this->next);

        if (info.mass > 0)
            this->consumedMassInfo.push_back(info);
        while (!this->next.empty())
            this->next.pop();
        data.offsetIdx++;
    }

    LOGGER_LEVEL(1, "ElementConsumer %s-done\n", __func__);
}

Handle ElementEmitter::Register(SimData* simData, AddElementEmitterMessage* msg)
{
    ElementEmitterData new_data = {
        .elapsedTime     = -1,
        .emitInterval    = 0,
        .emitMass        = FLT_MAX,
        .emitTemperature = 0,
        .maxPressure     = msg->maxPressure,
        .cell            = Vector2<uint16_t>(0, 0),
        .elemIdx         = 0,
        .maxDepth        = 0,
        .offsetIdx       = 0,
        .blockedState    = 0xFF,
        .blockedCBIdx    = (uint32_t)msg->onBlockedCBIdx,
        .unblockedCBIdx  = (uint32_t)msg->onUnblockedCBIdx
    };

    Handle result = this->AddData(&new_data);

    int item = this->handles.items[result & 0xFFFFFF];
    if (this->emittedMassInfo.size() > item) {
        EmittedMassInfo info = {
            .elemIdx      = 0xFFFF,
            .diseaseIdx   = 0xFF,
            .mass         = 0,
            .temperature  = 0,
            .diseaseCount = 0,
        };
        this->emittedMassInfo[item] = info;
    }
    return result;
}

void ElementEmitter::Unregister(int handle)
{
    this->Free(handle);
}

void ElementEmitter::Modify(SimData* simData, ModifyElementEmitterMessage* msg)
{
    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) return;

    int item = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[item].elapsedTime      = 0;
    this->data[item].emitInterval     = msg->emitInterval;
    this->data[item].emitMass         = msg->emitMass;
    this->data[item].emitTemperature  = msg->emitTemperature;
    this->data[item].maxPressure      = msg->maxPressure;
    this->data[item].emitDiseaseCount = msg->diseaseCount;
    this->data[item].cell.x           = msg->gameCell % (simData->width - 2) + 1;
    this->data[item].cell.y           = msg->gameCell / (simData->width - 2) + 1;
    this->data[item].elemIdx          = msg->elementIdx;
    this->data[item].maxDepth         = msg->maxDepth;
    this->data[item].diseaseIdx       = msg->diseaseIdx;
    if (msg->emitMass <= 0)
        this->data[item].blockedState = -1;
}

void ElementEmitter::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "ElementEmitter %s-%llu\n", __func__, this->data.size());
    this->emittedMassInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        ElementEmitterData& data = this->data[idx];

        if (data.cell.x < region->minimum.x) continue;
        if (data.cell.y < region->minimum.y) continue;
        if (data.cell.x > region->maximum.x) continue;
        if (data.cell.y > region->maximum.y) continue;

        if (data.elapsedTime >= data.emitInterval) {
            int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
            uint16_t elemIdx = this->emittedMassInfo[handleVal].elemIdx;
            if (elemIdx == 0xFFFF || elemIdx == simData->vacuumElementIdx || elemIdx == data.elemIdx) {
                GetReachableCells(simData, data.cell.x, data.cell.y, data.maxDepth, &this->visitedCells, &this->next, &this->reachableCells, true);
                bool emitFlag = false;
                for (int cell : this->reachableCells) {
                    if (data.maxPressure > simData->updatedCells->mass[cell]) {
                        emitFlag = true;
                        break;
                    }
                }
                int callbackIdx = -1;
                if (emitFlag && data.blockedState) {
                    data.blockedState = 0;
                    callbackIdx = data.unblockedCBIdx;
                }
                else if (data.blockedState != 1) {
                    data.blockedState = 1;
                    callbackIdx = data.blockedCBIdx;
                }
                if (callbackIdx != -1) {
                    CallbackInfo info{ .callbackIdx = callbackIdx };
                    simData->simEvents->callbackInfo.push_back(info);
                }

                EmittedMassInfo info = this->TryEmit(&data, this->reachableCells, data.offsetIdx, simData);
                if (info.mass > 0) {
                    this->emittedMassInfo[handleVal].temperature = CalculateFinalTemperature(
                        this->emittedMassInfo[handleVal].mass, this->emittedMassInfo[handleVal].temperature, info.mass, info.temperature);
                    this->emittedMassInfo[handleVal].elemIdx     = info.elemIdx;
                    this->emittedMassInfo[handleVal].mass       += info.mass;
                }
                this->visitedCells.clear();
                this->reachableCells.clear();
            }
            data.elapsedTime -= data.emitInterval;
        }
        data.elapsedTime += dt;
    }

    LOGGER_LEVEL(1, "ElementEmitter %s-done\n", __func__);
}

void ElementEmitter::UpdateDataListOnly(SimData* simData)
{
    this->emittedMassInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
        this->emittedMassInfo[handleVal].elemIdx      = -1;
        this->emittedMassInfo[handleVal].diseaseIdx   = -1;
        this->emittedMassInfo[handleVal].mass         = 0;
        this->emittedMassInfo[handleVal].temperature  = 0;
        this->emittedMassInfo[handleVal].diseaseCount = 0;
    }
}

EmittedMassInfo ElementEmitter::Emit(ElementEmitterData* data, int cell_idx, float tempEmit, SimData* simData)
{
    if (simData->updatedCells->elementIdx[cell_idx] == data->elemIdx ||
        simData->updatedCells->elementIdx[cell_idx] == simData->vacuumElementIdx)
    {
        float tempFin = CalculateFinalTemperature(simData->updatedCells->mass[cell_idx], simData->updatedCells->temperature[cell_idx], data->emitMass, tempEmit);
        simData->updatedCells->elementIdx [cell_idx]  = data->elemIdx;
        simData->updatedCells->temperature[cell_idx]  = tempFin;
        simData->updatedCells->mass       [cell_idx] += data->emitMass;
        gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell_idx, data->diseaseIdx, data->emitDiseaseCount);
        ASSERT_TEMP(simData->updatedCells->mass[cell_idx], simData->updatedCells->temperature[cell_idx]);

        if (simData->updatedCells->elementIdx[cell_idx] == simData->vacuumElementIdx) {
            int cellInner = CELL_S2G(cell_idx, simData->width);
            if (CELL_AVAILABLE(cellInner, simData)) {
                SubstanceChangeInfo info{ cellInner, (uint16_t)-1, (uint16_t)-1 };
                simData->simEvents->substanceChangeInfo.push_back(info);
            }
            simData->timers[cell_idx].stableCellTicks |= 0x1F;
        }
    }
    EmittedMassInfo result = {
        .elemIdx      = data->elemIdx,
        .diseaseIdx   = 0xFF,
        .mass         = data->emitMass,
        .temperature  = tempEmit,
        .diseaseCount = 0
    };
    return result;
}

EmittedMassInfo ElementEmitter::TryEmit(ElementEmitterData* data, std::vector<int>& reachable_cells, int offset, SimData* simData)
{
    float   emitTemperature = data->emitTemperature >= 0 ? data->emitTemperature : gElements[data->elemIdx].defaultValues.temperature;
    uint8_t state = gElements[data->elemIdx].state & 3;

    if (state != 0) {
        for (int idx = 0; idx < reachable_cells.size(); idx++) {
            int      cell = reachable_cells[(idx + offset) % reachable_cells.size()];
            uint16_t elem = simData->updatedCells->elementIdx[cell];
            if (data->maxPressure <= simData->updatedCells->mass[cell])
                continue;

            if (state == 1) {
                if (elem == data->elemIdx ||
                    elem == simData->vacuumElementIdx ||
                    DisplaceGas(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elem))
                {
                    return this->Emit(data, cell, emitTemperature, simData);
                }
            }
            else if (state == 2) {
                if (elem == data->elemIdx ||
                    elem == simData->vacuumElementIdx ||
                    DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elem) ||
                    DisplaceGas(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elem))
                {
                    return this->Emit(data, cell, emitTemperature, simData);
                }
            }
            else if (state == 3) {
                int cellInner = CELL_S2G(cell, simData->width);
                if (data->emitMass <= 0 && CELL_AVAILABLE(cellInner, simData) && CELL_VISIABLE(cellInner, simData)) {
                    SpawnOreInfo info = {
                        .cellIdx      = cellInner,
                        .elemIdx      = data->elemIdx,
                        .diseaseIdx   = 0xFF,
                        .mass         = data->emitMass,
                        .temperature  = emitTemperature,
                        .diseaseCount = 0
                    };
                    simData->simEvents->spawnOreInfo.push_back(info);
                }
#ifdef __DEBUGED__
                break;
#endif
            }
        }
    }

    EmittedMassInfo info = {
        .elemIdx      = 0xFFFF,
        .diseaseIdx   = 0xFF,
        .mass         = 0,
        .temperature  = 0,
        .diseaseCount = 0,
    };
    return info;
}

Handle RadiationEmitter::Register(SimData* simData, AddRadiationEmitterMessage* msg)
{
    RadiationEmitterData new_data = {
        .cell          = CELL_G2S(msg->gameCell, simData->width),
        .radiusX       = msg->radiusX,
        .radiusY       = msg->radiusY,
        .emitRads      = msg->emitRads,
        .emitRate      = msg->emitRate,
        .emitSpeed     = MIN_F(msg->emitSpeed, msg->emitRads),
        .emitDirection = msg->emitDirection,
        .emitAngle     = msg->emitAngle,
        .emitTimer     = 0,
        .emitStepTimer = 0,
        .emitType      = (RadiationEmitterType)msg->emitType,
        .emitStep      = 0
    };

    return this->AddData(&new_data);;
}

void RadiationEmitter::Unregister(int handle)
{
    this->Free(handle);
}

void RadiationEmitter::Modify(SimData* simData, ModifyRadiationEmitterMessage* msg)
{
    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) return;

    int item = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[item].cell          = CELL_G2S(msg->gameCell, simData->width);
    this->data[item].radiusX       = msg->radiusX;
    this->data[item].radiusY       = msg->radiusY;
    this->data[item].emitRads      = msg->emitRads;
    this->data[item].emitRate      = msg->emitRate;
    this->data[item].emitSpeed     = MIN_F(msg->emitRate, msg->emitSpeed);
    this->data[item].emitDirection = msg->emitDirection;
    this->data[item].emitAngle     = msg->emitAngle;
    this->data[item].emitType      = (RadiationEmitterType)msg->emitType;
}

void RadiationEmitter::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "RadiationEmitter %s-%llu\n", __func__, this->data.size());
    if (!simData->radiationEnabled) return;
    for (int idx = 0; idx < this->data.size(); idx++) {
        RadiationEmitterData& data = this->data[idx];

        int locX = data.cell % simData->width;
        int locY = data.cell / simData->width;
        if (locX < region->minimum.x) continue;
        if (locY < region->minimum.y) continue;
        if (locX > region->maximum.x) continue;
        if (locY > region->maximum.y) continue;
        if (data.emitRads <= 0)       continue;
        if (data.radiusX <= 0 && data.radiusY <= 0) continue;

        data.emitTimer += dt;
        bool isUpdate = data.emitRate == 0 || data.emitTimer > data.emitRate;
        if (isUpdate || data.emitTimer <= data.emitSpeed) {
            if (isUpdate) {
                data.emitTimer = 0;
                data.emitStep  = 0;
            }
            int   radius   = std::max(data.radiusX, data.radiusY);
            float stepLoop = dt * radius / data.emitSpeed;
            if (data.emitRate == 0 || stepLoop <= 1 && data.emitSpeed <= radius * data.emitStepTimer)
                stepLoop = 1;
            for (int loop = 0; loop < stepLoop; loop++) {
                switch (data.emitType) {
                case RadiationEmitterType::Constant:
                    tickConstant(simData, region, &data, dt);
                    break;
                case RadiationEmitterType::Pulsing:
                case RadiationEmitterType::PulsingAveraged:
                    tickPulsing(simData, region, &data, dt, radius);
                    break;
                case RadiationEmitterType::SimplePulse:
                    tickSimplePulse(simData, region, &data, dt, radius);
                    break;
                case RadiationEmitterType::Attractor:
                    tickAttractor(simData, region, &data, dt, radius);
                    break;
                default:
                    break;
                }
                data.emitStep       = (data.emitStep + 1) % radius;
                data.emitStepTimer -=  data.emitSpeed     / radius;
            }
            data.emitStepTimer += dt;
        }
    }

    LOGGER_LEVEL(1, "RadiationEmitter %s-done\n", __func__);
}

Handle DiseaseEmitter::Register(SimData* simData, AddDiseaseEmitterMessage* msg)
{
    DiseaseEmitterData new_data = {
        .cell         = -1,
        .diseaseIdx   = 0xFF,
        .range        = 0,
        .emitCount    = 0,
        .emitInterval = 0,
        .elapsedTime  = 0
    };

    return this->AddData(&new_data);;
}

void DiseaseEmitter::Unregister(int handle)
{
    this->Free(handle);
}

void DiseaseEmitter::Modify(SimData* simData, ModifyDiseaseEmitterMessage* msg)
{
    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) return;

    int item = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[item].cell         = CELL_G2S(msg->gameCell, simData->width);
    this->data[item].diseaseIdx   = msg->diseaseIdx;
    this->data[item].range        = msg->maxDepth;
    this->data[item].emitCount    = msg->emitCount;
    this->data[item].emitInterval = msg->emitInterval;
}

void DiseaseEmitter::Update(float dt, SimData* simData, Region* region)
{
    LOGGER_LEVEL(1, "DiseaseEmitter %s-%llu\n", __func__, this->data.size());
    this->emittedInfo.resize(this->handles.items.size());
    for (DiseaseEmitterData& data : this->data) {
        int locX = data.cell % simData->width;
        int locY = data.cell / simData->width;
        if (data.diseaseIdx == 0xFF)  continue;
        if (locX < region->minimum.x) continue;
        if (locY < region->minimum.y) continue;
        if (locX > region->maximum.x) continue;
        if (locY > region->maximum.y) continue;

        if (data.elapsedTime >= data.emitInterval) {
            GetReachableCells(simData, locX, locY, data.range, &this->visitedCells, &this->next, &this->reachableCells, false);
            for (int cell : this->reachableCells)
                if (simData->updatedCells->diseaseIdx[cell] == data.diseaseIdx)
                    gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell, data.diseaseIdx, data.emitCount);
#ifdef __DEBUGED__
                else if (simData->updatedCells->diseaseIdx[cell] == 0xFF && data.emitCount > 0)
                    gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell, data.diseaseIdx, data.emitCount);
#endif
            data.elapsedTime -= data.emitInterval;
        }
        data.elapsedTime += dt;
    }

    LOGGER_LEVEL(1, "DiseaseEmitter %s-done\n", __func__);
}

void DiseaseEmitter::UpdateDataListOnly(SimData* simData)
{
    this->emittedInfo.resize(this->handles.items.size());
    for (int idx = 0; idx < this->data.size(); idx++) {
        int handleVal = this->dataHandleIndices[idx] & 0xFFFFFF;
        this->emittedInfo[handleVal].diseaseIdx = -1;
        this->emittedInfo[handleVal].count      = 0;
    }
}

Handle DiseaseConsumer::Register(SimData* simData, AddDiseaseConsumerMessage* msg)
{
    DiseaseConsumerData new_data = DiseaseConsumerData();
    return this->AddData(&new_data);
}

void DiseaseConsumer::Unregister(int handle)
{
    this->Free(handle);
}

//__GUESS__
// Unused parameter: emitInterval
void DiseaseConsumer::Modify(SimData* simData, ModifyDiseaseConsumerMessage* msg)
{
    if (!HANDLE_AVAILABLE(msg->handle, this->handles.versions)) return;

    int item = this->handles.items[msg->handle & 0xFFFFFF];
    this->data[item].cell.x          = msg->gameCell % (simData->width - 2) + 1;
    this->data[item].cell.y          = msg->gameCell / (simData->width - 2) + 1;
    this->data[item].consumptionRate = (float)msg->emitCount;
    this->data[item].maxDepth        = msg->maxDepth;

    this->consumedInfo[msg->handle & 0xFFFFFF].diseaseIdx = msg->diseaseIdx;
}