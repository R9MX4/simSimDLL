#include "ClassCell.h"
#include "ClassSim.h"
#include <time.h>

SimData::SimData(int sim_width, int sim_height, bool radiation_enabled, bool headless)
{
    int cellCount = sim_width * sim_height;
    this->width   = sim_width;
    this->height  = sim_height;
    this->numGameCells = (sim_height - 2) * (sim_width - 2);
    this->savedOptions = 0;
    this->cells       .reset();
    this->updatedCells.reset();
    this->worldZones  .reset();
    this->cosmicRadiationOcclusion = std::vector<float>();
    this->accumulatedFlow .reset();
    this->flow            .reset();
    this->biomeTemperature.reset();
    this->timers          .reset();
    this->activeRegions    = std::vector<ActiveRegion>();
    this->visibleGrid     .reset();
    this->tickCount        = 0;
    this->simEvents       .reset();
    this->debugProperties  = {
        .buildingTemperatureScale           = 100.0f, // 0.001f,
        .buildingToBuildingTemperatureScale = 0.001f,
        .biomeTemperatureLerpRate           = 0.001f };
    this->components                     = std::vector<SimComponent*>();
    this->elementConsumer                = ElementConsumer();
    this->elementEmitter                 = ElementEmitter();
    this->radiationEmitter               = RadiationEmitter();
    this->elementChunk                   = ElementChunk();
    this->buildingHeatExchange           = BuildingHeatExchange();
    this->buildingToBuildingHeatExchange = BuildingToBuildingHeatExchange();
    this->diseaseEmitter                 = DiseaseEmitter();
    this->diseaseConsumer                = DiseaseConsumer();
    this->initSettleThermalBoundaries    = 0;
    this->headless                       = headless;
    this->radiationEnabled               = radiation_enabled;
    this->RADIATION_LINGER_RATE          = 1.1f;
    this->RADIATION_MAX_MASS             = 2000.0f;
    this->RADIATION_BASE_WEIGHT          = 0.3f;
    this->RADIATION_DENSITY_WEIGHT       = 0.7f;
    this->RADIATION_CONSTRUCTED_FACTOR   = 0.8f;
    this->worlds = std::vector<SimData::WorldOffsetData>();

    this->cells            = std::make_unique<CellSOA>(cellCount);
    this->updatedCells     = std::make_unique<CellSOA>(cellCount);
    this->accumulatedFlow  = std::make_unique<float[]>          (cellCount);
    this->flow             = std::make_unique<Vector4<float>[]> (cellCount);
    this->biomeTemperature = std::make_unique<float[]>          (cellCount);
    this->timers           = std::make_unique<SimData::Timers[]>(cellCount);
    this->visibleGrid      = std::make_unique<uint8_t[]>(this->numGameCells);
    this->simEvents        = std::make_unique<SimEvents>();

    this->voidElementIdx        = GetElementIndex(0xA9360B34);
    this->vacuumElementIdx      = GetElementIndex(0x2D39BF75);
    this->unobtaniumElementIdx  = GetElementIndex(0x6D95058C);
    this->randomSeed            = (uint32_t)_time64(0);

    for (int cell = 0; cell < cellCount; cell++) {
        this->timers[cell].stableCellTicks |= 0x1F;
    }

    this->components.push_back(&this->elementConsumer);
    this->components.push_back(&this->elementEmitter);
    this->components.push_back(&this->radiationEmitter);
    this->components.push_back(&this->elementChunk);
    this->components.push_back(&this->buildingHeatExchange);
    this->components.push_back(&this->buildingToBuildingHeatExchange);
    this->components.push_back(&this->diseaseEmitter);
    this->components.push_back(&this->diseaseConsumer);
}

void SimEvents::ChangeSubstance(const SimData* simData, uint64_t sim_cell)
{
    int innerCell = (int)CELL_S2G(sim_cell, simData->width);
    if (CELL_AVAILABLE(innerCell, simData)) {
        SubstanceChangeInfo info = { innerCell, (uint16_t)-1, (uint16_t)-1 };
        this->substanceChangeInfo.push_back(info);
    }
    simData->timers[sim_cell].stableCellTicks |= 0x1F; // simData->timers[sim_cell] |= 0x1Fu;
}

bool SimEvents::SpawnOre(const SimData* simData, uint64_t sim_cell, uint16_t elemIdx, float mass, 
    float temperature, uint8_t disease_idx, int disease_count, bool spawn_even_if_hidden)
{
    ASSERT_TEMP(mass, temperature);
    int innerCell = (int)CELL_S2G(sim_cell, simData->width);

    if (mass <= 0.0)                                               return false;
    if (!CELL_AVAILABLE(innerCell, simData))                       return false;
    if (!spawn_even_if_hidden && !simData->visibleGrid[innerCell]) return false;

    SpawnOreInfo spawnOreInfo = {
        .cellIdx      = innerCell,
        .elemIdx      = elemIdx,
        .diseaseIdx   = disease_idx,
        .mass         = mass,
        .temperature  = temperature,
        .diseaseCount = disease_count
    };
    this->spawnOreInfo.push_back(spawnOreInfo);
    return true;
}

bool SimEvents::SpawnFallingLiquid(const SimData* simData, uint64_t sim_cell, uint16_t elemIdx, float mass,
    float temperature, uint8_t disease_idx, int disease_count, bool spawn_even_if_hidden)
{
    //LOGGER_LEVEL(1, "%s-cell:%d, elem:%d, temp:%f, mass:%f\n", __func__, (int)sim_cell, elemIdx, temperature, mass);
    if (simData->headless || (simData->cells->properties[sim_cell - simData->width] & 2) != 0)
        return false;
    int innerCell = (int)CELL_S2G(sim_cell, simData->width);
    if (!CELL_AVAILABLE(innerCell, simData))                       return false;
    if (!spawn_even_if_hidden && !simData->visibleGrid[innerCell]) return false;
    if (temperature < gElements[elemIdx].lowTemp  - 3)             return false;
    if (temperature > gElements[elemIdx].highTemp + 3)             return false;

    SpawnFallingLiquidInfo spawnFallingLiquidInfo = {
        .cellIdx      = innerCell,
        .elementIdx   = elemIdx,
        .diseaseIdx   = disease_idx,
        .mass         = mass,
        .temperature  = temperature,
        .diseaseCount = disease_count
    };
    this->spawnLiquidInfo.push_back(spawnFallingLiquidInfo);
    return true;
}

void SimEvents::SpawnFX(const SimData* simData, uint64_t sim_cell, int fxid, float rotation)
{
    int innerCell = (int)CELL_S2G(sim_cell, simData->width);
    if (!CELL_AVAILABLE(innerCell, simData))
        return;

    SpawnFXInfo spawnFXInfo = {
        .cellIdx   = innerCell,
        .fxID      = fxid,
        .rotation = rotation
    };
    this->spawnFXInfo.push_back(spawnFXInfo);
}

void SimData::ApplySaveSettings(uint8_t newSaveOptions)
{
    this->savedOptions = newSaveOptions;
}

void SimData::ClearCell(CellAccessor* cell)
{
    cell->cells->elementIdx  [cell->cellIdx] = this->vacuumElementIdx;
    cell->cells->mass        [cell->cellIdx] = 0.0f;
    cell->cells->temperature [cell->cellIdx] = 0.0f;
    cell->cells->diseaseIdx  [cell->cellIdx] = -1;
    cell->cells->diseaseCount[cell->cellIdx] = 0;
    cell->cells->diseaseInfestationTickCount  [cell->cellIdx] = 0;
    cell->cells->diseaseGrowthAccumulatedError[cell->cellIdx] = 0;
}

uint8_t SimData::GetStableTicksRemaining(uint64_t sim_cell)
{
    uint8_t result = this->timers[sim_cell].stableCellTicks & 0x1F;
    if (result == 31) {
        this->timers[sim_cell].stableCellTicks = RAND_INT(this->randomSeed, 3, 6) & 0x1F;
    }
    else if (result) {
        this->timers[sim_cell].stableCellTicks = --result & 0x1F;
    }
    return this->timers[sim_cell].stableCellTicks;
}

void SimData::SetBoundary(int cell)
{
    this->cells->elementIdx     [cell] = gSimData->unobtaniumElementIdx;
    this->cells->mass           [cell] = 9999.0;
    this->cells->temperature    [cell] = 0.0;
    this->cells->diseaseCount   [cell] = 0;
    this->cells->diseaseIdx     [cell] = -1;

    this->updatedCells->elementIdx  [cell] = gSimData->unobtaniumElementIdx;
    this->updatedCells->mass        [cell] = 9999.0;
    this->updatedCells->temperature [cell] = 0.0;
    this->updatedCells->diseaseCount[cell] = 0;
    this->updatedCells->diseaseIdx  [cell] = -1;

    this->biomeTemperature[cell] = 0.0;
}

void SimData::InitializeBoundary()
{
    if (this->width > 0) {
        for (int idx_x = 0; idx_x < this->width; idx_x++) {
            this->SetBoundary(idx_x);
            this->SetBoundary(idx_x + this->width * (this->height - 1));
        }
    }
    if (this->height - 1 > 1) {
        for (int idx_y = 0; idx_y < this->height; idx_y++) {
            this->SetBoundary( idx_y      * this->width);
            this->SetBoundary((idx_y + 1) * this->width - 1);
        }
    }
}

void SimData::SettleThermalBoundaries(CellSOA* src_cells, CellSOA* dest_cells)
{
    int ncell_offsets[8] = { -gSimData->width - 1, -gSimData->width, -gSimData->width + 1 ,
                      -1, 1 , gSimData->width - 1,  gSimData->width,  gSimData->width + 1 };
    uint16_t ElementIndex = GetElementIndex(0x3FE0146Eu);
    for (int idx_y = 1; idx_y < gSimData->height - 2; idx_y++) {
        for (int idx_x = 1; idx_x < gSimData->width - 1; idx_x++) {
            int cell = idx_x + this->width * idx_y;
            if (src_cells->elementIdx[cell] != ElementIndex) continue;

            float cellCount = 0;
            float tempTotal = 0;
            for (int offset : ncell_offsets) {
                int cell_n = cell + offset;
                if (src_cells->elementIdx[cell_n] != ElementIndex
                    && src_cells->elementIdx[cell_n] != GetElementIndex(0x6D95058Cu)
                    && src_cells->mass[cell_n] > 0.0)
                {
                    cellCount = cellCount + 1.0f;
                    tempTotal = tempTotal + src_cells->temperature[cell_n];
                }
            }
            if (cellCount > 0) 
                dest_cells->temperature[cell] = tempTotal / cellCount;
        }
    }
}

void SimData::UpdateComponents(float dt, Region* region)
{
    LOGGER_PRINT2("%s-%llu\n", __func__, this->components.size());
    if (this->components.size() <= 0) return;

    for (SimComponent* component : this->components) {
        component->Update(dt, this, region);
    }
    LOGGER_PRINT2("%s done\n", __func__);
}

void SimData::UpdateComponentsDataListOnly()
{
    LOGGER_PRINT2("%s-%llu\n", __func__, this->components.size());
    if (this->components.size() <= 0) return;

    for (SimComponent* component : this->components) {
        component->UpdateDataListOnly(this);
    }
    LOGGER_PRINT2("%s done\n", __func__);
}

int SimData::FreeGridCells(BinaryBufferReader* reader)
{
    LOGGER_PRINT2("%s\n", __func__);
    int width, height, x_offset, y_offset;
    *reader >> width >> height >> x_offset >> y_offset;

    for (int locY = 0; locY < height + 2; locY++) {
        for (int locX = 0; locX < width + 2; locX++) {
            int cell = gSimData->width * (locY + y_offset) + locX + x_offset;

            gSimData->cells->elementIdx  [cell] = gSimData->vacuumElementIdx;
            gSimData->cells->mass        [cell] = 0;
            gSimData->cells->temperature [cell] = 0;
            gSimData->cells->diseaseCount[cell] = 0;
            gSimData->cells->diseaseIdx  [cell] = 0xFF;
            gSimData->cells->radiation   [cell] = 0;
            gSimData->updatedCells->elementIdx  [cell] = gSimData->vacuumElementIdx;
            gSimData->updatedCells->mass        [cell] = 0;
            gSimData->updatedCells->temperature [cell] = 0;
            gSimData->updatedCells->diseaseCount[cell] = 0;
            gSimData->updatedCells->diseaseIdx  [cell] = 0xFF;
            gSimData->updatedCells->radiation   [cell] = 0;
            gSimData->biomeTemperature[cell] = 0;
        }
    }
    for (auto it = gSimData->worlds.begin(); it != gSimData->worlds.end(); it++) {
        if (it->offsetX == x_offset && it->offsetY == y_offset) {
            gSimData->worlds.erase(it);
            return 0;
        }
    }
    LOGGER_PRINT2("%s done\n", __func__);
    return 0;
}

int SimData::InitializeFromCells(BinaryBufferReader* reader)
{
#define CALLOC_And_READ(_DEST, _TYPE, _LENGTH, _READER) \
		_DEST = (_TYPE*)calloc(_LENGTH, sizeof(_TYPE)); if (_DEST == NULL) ASSERT_TEXT("calloc memory failed"); _READER->ReadBytes(sizeof(_TYPE) * _LENGTH, _DEST);

    InitializeFromCellsMessage msg;
    *reader >> msg.width >> msg.height >> msg.randSeed >> msg.radiation_enabled >> msg.headless;
    LOGGER_PRINT("%s width: %d, height: %d, Radiaton:%s, Headless:%s\n", __func__, msg.width, msg.height, msg.radiation_enabled ? "Yes" : "No", msg.headless ? "True" : "False");

    int cellCnt = msg.width * msg.height;
    CALLOC_And_READ(msg.cells,            InitializeFromCellsMessage::Cell,        cellCnt, reader);
    CALLOC_And_READ(msg.biomeTemperature, float,                                   cellCnt, reader);
    CALLOC_And_READ(msg.disease,          InitializeFromCellsMessage::DiseaseCell, cellCnt, reader);

    gSimData = std::make_unique<SimData>(msg.width + 2, msg.height + 2, msg.radiation_enabled, msg.headless);

    size_t length = sizeof(float) * msg.width;
    for (int locY = 0; locY < msg.height; locY++) {
        int cell      = gSimData->width * (locY + 1) + 1;
        int cellInner = msg.width * locY;
        memmove(&gSimData->biomeTemperature[cell], &msg.biomeTemperature[cellInner], length);
    }

    for (int locY = 0; locY < msg.height; locY++) {
        for (int locX = 0; locX < msg.width; locX++) {
            int cell      = gSimData->width * (locY + 1) + locX + 1;
            int cellInner = msg.width * locY + locX;
            gSimData->updatedCells->elementIdx  [cell] = msg.cells[cellInner].elementIdx;
            gSimData->updatedCells->temperature [cell] = msg.cells[cellInner].temperature;
            gSimData->updatedCells->mass        [cell] = msg.cells[cellInner].mass;
            gSimData->updatedCells->properties  [cell] = msg.cells[cellInner].properties;
            gSimData->updatedCells->insulation  [cell] = msg.cells[cellInner].insulation;
            gSimData->updatedCells->strengthInfo[cell] = msg.cells[cellInner].strengthInfo;
            gSimData->updatedCells->radiation   [cell] = 0;

            gSimData->updatedCells->diseaseIdx                   [cell] = msg.disease[cellInner].diseaseIdx;
            gSimData->updatedCells->diseaseCount                 [cell] = msg.disease[cellInner].count;
            gSimData->updatedCells->diseaseInfestationTickCount  [cell] = msg.disease[cellInner].infestationTickCount;
            gSimData->updatedCells->diseaseGrowthAccumulatedError[cell] = 0;
        }
    }

    gSimData->SettleThermalBoundaries(gSimData->cells.get(), gSimData->updatedCells.get());
    gSimData->InitializeBoundary();
    LOGGER_PRINT2("%s done\n", __func__);
    return cellCnt;
}

int SimData::ResizeAndInitializeVacuumCells(BinaryBufferReader* reader)
{
    LOGGER_PRINT2("%s\n", __func__);
    ResizeAndInitializeVacuumCellsMessage msg;
    *reader >> msg.grid_width >> msg.grid_height >> msg.width >> msg.height >> msg.x_offset >> msg.y_offset;

    for (int locY = 0; locY < msg.height; locY++) {
        for (int locX = 0; locX < msg.width; locX++) {
            int cell = gSimData->width * (locY + msg.y_offset + 1) + locX + msg.x_offset + 1;
            
            gSimData->cells->elementIdx  [cell] = gSimData->vacuumElementIdx;
            gSimData->cells->mass        [cell] = 0;
            gSimData->cells->temperature [cell] = 0;
            gSimData->cells->diseaseCount[cell] = 0;
            gSimData->cells->diseaseIdx  [cell] = 0xFF;
            gSimData->cells->radiation   [cell] = 0;
            gSimData->updatedCells->elementIdx  [cell] = gSimData->vacuumElementIdx;
            gSimData->updatedCells->mass        [cell] = 0;
            gSimData->updatedCells->temperature [cell] = 0;
            gSimData->updatedCells->diseaseCount[cell] = 0;
            gSimData->updatedCells->diseaseIdx  [cell] = 0xFF;
            gSimData->updatedCells->radiation   [cell] = 0;
            gSimData->biomeTemperature[cell] = -1;
        }
    }
    for (int locX = 0; locX < msg.width + 2; locX++) {
        int cellB = gSimData->width * (msg.y_offset                 ) + locX + msg.x_offset;
        int cellT = gSimData->width * (msg.y_offset + msg.height + 1) + locX + msg.x_offset;

        gSimData->cells->elementIdx  [cellB] = gSimData->unobtaniumElementIdx;
        gSimData->cells->mass        [cellB] = 9999;
        gSimData->cells->temperature [cellB] = 0;
        gSimData->cells->diseaseCount[cellB] = 0;
        gSimData->cells->diseaseIdx  [cellB] = 0xFF;
        gSimData->cells->radiation   [cellB] = 0;
        gSimData->updatedCells->elementIdx  [cellB] = gSimData->unobtaniumElementIdx;
        gSimData->updatedCells->mass        [cellB] = 9999;
        gSimData->updatedCells->temperature [cellB] = 0;
        gSimData->updatedCells->diseaseCount[cellB] = 0;
        gSimData->updatedCells->diseaseIdx  [cellB] = 0xFF;
        gSimData->updatedCells->radiation   [cellB] = 0;
        gSimData->biomeTemperature[cellB] = 0;

        gSimData->cells->elementIdx  [cellT] = gSimData->unobtaniumElementIdx;
        gSimData->cells->mass        [cellT] = 9999;
        gSimData->cells->temperature [cellT] = 0;
        gSimData->cells->diseaseCount[cellT] = 0;
        gSimData->cells->diseaseIdx  [cellT] = 0xFF;
        gSimData->cells->radiation   [cellT] = 0;
        gSimData->updatedCells->elementIdx  [cellT] = gSimData->unobtaniumElementIdx;
        gSimData->updatedCells->mass        [cellT] = 9999;
        gSimData->updatedCells->temperature [cellT] = 0;
        gSimData->updatedCells->diseaseCount[cellT] = 0;
        gSimData->updatedCells->diseaseIdx  [cellT] = 0xFF;
        gSimData->updatedCells->radiation   [cellT] = 0;
        gSimData->biomeTemperature[cellT] = 0;
    }
    for (int locY = 1; locY < msg.height + 1; locY++) {
        int cellL = gSimData->width * (locY + msg.y_offset    ) + msg.x_offset;
        int cellR = cellL + msg.width;

        gSimData->cells->elementIdx  [cellL] = gSimData->unobtaniumElementIdx;
        gSimData->cells->mass        [cellL] = 9999;
        gSimData->cells->temperature [cellL] = 0;
        gSimData->cells->diseaseCount[cellL] = 0;
        gSimData->cells->diseaseIdx  [cellL] = 0xFF;
        gSimData->cells->radiation   [cellL] = 0;
        gSimData->updatedCells->elementIdx  [cellL] = gSimData->unobtaniumElementIdx;
        gSimData->updatedCells->mass        [cellL] = 9999;
        gSimData->updatedCells->temperature [cellL] = 0;
        gSimData->updatedCells->diseaseCount[cellL] = 0;
        gSimData->updatedCells->diseaseIdx  [cellL] = 0xFF;
        gSimData->updatedCells->radiation   [cellL] = 0;
        gSimData->biomeTemperature[cellL] = 0;

        gSimData->cells->elementIdx  [cellR] = gSimData->unobtaniumElementIdx;
        gSimData->cells->mass        [cellR] = 9999;
        gSimData->cells->temperature [cellR] = 0;
        gSimData->cells->diseaseCount[cellR] = 0;
        gSimData->cells->diseaseIdx  [cellR] = 0xFF;
        gSimData->cells->radiation   [cellR] = 0;
        gSimData->updatedCells->elementIdx  [cellR] = gSimData->unobtaniumElementIdx;
        gSimData->updatedCells->mass        [cellR] = 9999;
        gSimData->updatedCells->temperature [cellR] = 0;
        gSimData->updatedCells->diseaseCount[cellR] = 0;
        gSimData->updatedCells->diseaseIdx  [cellR] = 0xFF;
        gSimData->updatedCells->radiation   [cellR] = 0;
        gSimData->biomeTemperature[cellR] = 0;
    }
    SimData::WorldOffsetData new_world = {
        .offsetX = msg.x_offset,
        .offsetY = msg.y_offset,
        .width   = msg.width,
        .height  = msg.height
    };
    gSimData->worlds.push_back(new_world);

    LOGGER_PRINT2("%s Done\n", __func__);
    return msg.height * msg.width;
}

int SimData::SetWorldZones(BinaryBufferReader* reader)
{
    int cellCnt = gSimData->width * gSimData->height;
    gSimData->worldZones = std::make_unique<uint8_t[]>(cellCnt);
    memset(gSimData->worldZones.get(), 0, sizeof(uint8_t) * cellCnt);

    for (int locY = 1; locY < gSimData->height - 1; locY++) {
        reader->ReadBytes(gSimData->width - 2, &gSimData->worldZones[gSimData->width * locY + 1]);
    }
    return 0;
}