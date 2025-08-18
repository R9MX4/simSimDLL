#include "ClassBase.h"
#include <iterator>

//----- (Variable) ---------------------------------------------------
int noffset    = -1;
int tick_count = 0;

//----- (Decompiled Basic) -------------------------------------------
void UpdateFlowTexture(const SimData* simData, GameData* new_game_data)
{
	for (int yIdx = 0; yIdx < new_game_data->height; yIdx++) {
		for (int xIdx = 0; xIdx < new_game_data->width; xIdx++) {
			int cellS = (yIdx + 1) * simData->width       + xIdx + 1;
			int cellN =  yIdx      * new_game_data->width + xIdx;
			float ratio = simData->updatedCells->elementIdx[cellS] == simData->cells->elementIdx[cellS] ? 
				1.0f / MAX_F(simData->updatedCells->mass[cellS], 1) : 0;
			new_game_data->propertyTextureFlow[cellN].x = ratio * (simData->flow[cellS].x - simData->flow[cellS].y);
			new_game_data->propertyTextureFlow[cellN].y = ratio * (simData->flow[cellS].w - simData->flow[cellS].z);
		}
	}
}

void UpdateLiquidPropertyTexture(const SimData* simData, GameData* new_game_data)
{
	for (int yIdx = 0; yIdx < new_game_data->height-1; yIdx++) {
		for (int xIdx = 0; xIdx < new_game_data->width; xIdx++) {
			int cellS = (yIdx + 1) * simData->width       + xIdx + 1;
			int cellN =  yIdx      * new_game_data->width + xIdx;
			uint16_t elemS = simData->updatedCells->elementIdx[cellS];
			uint16_t elemT = simData->updatedCells->elementIdx[cellS + simData->width];

			if ((gElementPropertyTextureData[elemS].state & 3) != 2) {
				new_game_data->propertyTextureLiquid    [cellN] = 0;
				new_game_data->propertyTextureLiquidData[cellN] = 0;
				continue;
			}

			// Top Layer or not
			float massFactor = (gElementPropertyTextureData[elemT].state & 3) < 2 ?
				powf(MIN_F(simData->updatedCells->mass[cellS] * 0.001f, 1), 0.45f) : 1;
			uint32_t texture = (int)(massFactor * 255) & 0xff;
			new_game_data->propertyTextureLiquid[cellN] = texture << 24 | gElementPropertyTextureData[elemS].colour & 0xffffff;

			ElementTemperatureData* tempData = &gElementTemperatureData[elemS];
			float lowTrans    = tempData->lowTemp  - 3;
			float highTrans   = tempData->highTemp + 3;
			float tempFactor  = (simData->updatedCells->temperature[cellS] - lowTrans) / (highTrans - lowTrans);
			tempFactor        = CLAMP_F(tempFactor, 1, 0);
			uint16_t transIdx = tempFactor < 0.5 ? tempData->lowTempTransitionIdx : tempData->highTempTransitionIdx;
			uint32_t colour   = transIdx == 0xFFFF ? 0xFFFFFF : gElementPropertyTextureData[transIdx].colour;
			uint32_t texture2 = (int)(tempFactor * 255) & 0xff;
			new_game_data->propertyTextureLiquidData[cellN] = texture2 << 24 | colour & 0xffffff;
		}
	}
}

void UpdateExposedToSunPropertyTexture(const SimData* simData, GameData* new_game_data)
{
	for (SimData::WorldOffsetData world : simData->worlds) {
		if (world.width <= 0 || world.height < 2)
			continue;

		int cellWorld      = simData->width * (world.height + world.offsetY) + world.offsetX + 1;
		int cellWorldInner = CELL_S2G(cellWorld, simData->width);

		// Set top layer Grid.ExposedToSunlight to 255
		memset(&new_game_data->propertyTextureExposedToSunlight[cellWorldInner], 255, world.width);
		std::vector<float> column_light_factors(world.width, 1);

		ElementLightAbsorptionData* glassData = &gElementLightAbsorptionData[GetElementIndex(0x2531469C)];
		for (int yIdx = world.height - 2; yIdx >= 1; yIdx--) {
			for (int xIdx = 0; xIdx < world.width; xIdx++) {
				int cell      = simData->width * (yIdx + world.offsetY + 1) + world.offsetX + xIdx + 1;
				int cellTop   = cell + simData->width;
				int cellInner = CELL_S2G(cell, simData->width);
				if ((simData->updatedCells->properties[cell] & 0x30) == 0x30) {
					// Transparent & Opaque
					float massFactor = MIN_F(glassData->massScale * simData->updatedCells->mass[cellTop], 1);
					column_light_factors[xIdx] = MAX_F(column_light_factors[xIdx] - massFactor * glassData->factor, 0);
				}
				else if (simData->updatedCells->properties[cell] & 0x20) {
					// Opaque
					column_light_factors[xIdx] = 0;
				}
				else {
					uint16_t elemTop = simData->updatedCells->elementIdx[cellTop];
					ElementLightAbsorptionData* data = &gElementLightAbsorptionData[elemTop];
					float massFactor = MIN_F(data->massScale * simData->updatedCells->mass[cellTop], 1);
					column_light_factors[xIdx] = MAX_F(column_light_factors[xIdx] - massFactor * data->factor, 0);
				}
				new_game_data->propertyTextureExposedToSunlight[cellInner] = uint8_t(column_light_factors[xIdx] * 255);
			}
		}
	}
}

void ReplaceElement(int sim_cell_idx, uint16_t element_idx, float mass, float temperature, 
	uint8_t disease_idx, int disease_count, SimData* simData)
{
	simData->updatedCells->elementIdx[sim_cell_idx] = element_idx;
	if (gElements[element_idx].state & 3) {
		simData->updatedCells->temperature[sim_cell_idx] = temperature;
		simData->updatedCells->mass       [sim_cell_idx] = mass;
	}
	else { // Vacuum
		simData->updatedCells->temperature[sim_cell_idx] = 0;
		simData->updatedCells->mass       [sim_cell_idx] = 0;
	}
	ASSERT_TEMP(simData->updatedCells->mass[sim_cell_idx], simData->updatedCells->temperature[sim_cell_idx]);
	simData->simEvents->ChangeSubstance(simData, sim_cell_idx);
	simData->updatedCells->diseaseIdx                   [sim_cell_idx] = disease_idx;
	simData->updatedCells->diseaseCount                 [sim_cell_idx] = disease_count;
	simData->updatedCells->diseaseInfestationTickCount  [sim_cell_idx] = 0;
	simData->updatedCells->diseaseGrowthAccumulatedError[sim_cell_idx] = 0;
}

void ReplaceAndDisplaceElement(int sim_cell_idx, uint16_t element_idx, float mass, float temperature,
	uint8_t disease_idx, int disease_count, SimData* simData)
{
	uint16_t elemCell = simData->updatedCells->elementIdx[sim_cell_idx];
	uint8_t  state    = gElements[elemCell].state & 3;
	if      (state == 1) DisplaceGas(   simData, simData->simEvents.get(), simData->updatedCells.get(), sim_cell_idx, elemCell);
	else if (state == 2) DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), sim_cell_idx, elemCell);

	simData->updatedCells->elementIdx[sim_cell_idx] = element_idx;
	if (gElements[element_idx].state & 3) {
		simData->updatedCells->temperature[sim_cell_idx] = CLAMP_F(temperature, SIM_MAX_TEMPERATURE, 0);
		simData->updatedCells->mass       [sim_cell_idx] = mass;
	}
	else { // Vacuum
		simData->updatedCells->temperature[sim_cell_idx] = 0;
		simData->updatedCells->mass       [sim_cell_idx] = 0;
	}
	ASSERT_TEMP(simData->updatedCells->mass[sim_cell_idx], simData->updatedCells->temperature[sim_cell_idx]);
	simData->simEvents->ChangeSubstance(simData, sim_cell_idx);
	simData->updatedCells->diseaseIdx                   [sim_cell_idx] = disease_idx;
	simData->updatedCells->diseaseCount                 [sim_cell_idx] = disease_count;
	simData->updatedCells->diseaseInfestationTickCount  [sim_cell_idx] = 0;
	simData->updatedCells->diseaseGrowthAccumulatedError[sim_cell_idx] = 0;
}

void AddGas(int sim_cell_idx, uint16_t element_idx, float massDelta, float temperature, uint8_t disease_idx, 
	int disease_count, SimData* simData)
{
	// Same element
	if (simData->updatedCells->elementIdx[sim_cell_idx] == element_idx) {
		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			sim_cell_idx, massDelta, temperature, disease_idx, disease_count);
		return;
	}
	// Try to move liquid/gas
	uint16_t elemTarg = simData->updatedCells->elementIdx[sim_cell_idx];
	uint8_t  state    = gElements[elemTarg].state & 3;
	if ( state == 0 ||
		(state == 1 && DisplaceGas   (simData, simData->simEvents.get(), simData->updatedCells.get(), sim_cell_idx, elemTarg)) ||
		(state == 2 && DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), sim_cell_idx, elemTarg)))
	{
		simData->updatedCells->elementIdx[sim_cell_idx] = element_idx;
		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			sim_cell_idx, massDelta, temperature, disease_idx, disease_count);
		simData->simEvents->ChangeSubstance(simData, sim_cell_idx);
		ASSERT_TEMP(simData->updatedCells->mass[sim_cell_idx], simData->updatedCells->temperature[sim_cell_idx]);
		return;
	}

	int candidate_cells[3] = { sim_cell_idx - 1, sim_cell_idx + 1, sim_cell_idx + simData->width };
	for (int cell : candidate_cells) {
		// Try to find vacuum in bigger range
		if (simData->updatedCells->elementIdx[cell] == simData->vacuumElementIdx) {
			simData->updatedCells->elementIdx  [cell] = element_idx;
			simData->updatedCells->mass        [cell] = massDelta;
			simData->updatedCells->temperature [cell] = temperature;
			simData->updatedCells->diseaseIdx  [cell] = disease_idx;
			simData->updatedCells->diseaseCount[cell] = disease_count;
			simData->simEvents->ChangeSubstance(simData, cell);
			return;
		}
		// Try to find same element in bigger range
		if (simData->updatedCells->elementIdx[cell] == element_idx) {
			AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), 
				cell, massDelta, temperature, disease_idx, disease_count);
			return;
		}
	}
	// Mass elimination
	if (massDelta >= simData->updatedCells->mass[sim_cell_idx]) {
		float massAdd = massDelta - simData->updatedCells->mass[sim_cell_idx];
		simData->updatedCells->elementIdx [sim_cell_idx] = element_idx;
		simData->updatedCells->mass       [sim_cell_idx] = 0.0;
		simData->updatedCells->temperature[sim_cell_idx] = 0.0;

		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			sim_cell_idx, massAdd, temperature, disease_idx, (int)(massAdd / massDelta * disease_count));
		ASSERT_TEMP(simData->updatedCells->mass[sim_cell_idx], simData->updatedCells->temperature[sim_cell_idx]);

		if (simData->updatedCells->mass[sim_cell_idx] <= FLT_MIN) {
			CellAccessor cellAccessor(simData->updatedCells.get(), sim_cell_idx);
			simData->ClearCell(&cellAccessor);
		}
		simData->simEvents->ChangeSubstance(simData, sim_cell_idx);
	}
	else {
		simData->updatedCells->mass[sim_cell_idx] -= massDelta;
	}
}

void AddLiquid(int sim_cell_idx, uint16_t element_idx, float massDelta, float temperature, uint8_t disease_idx, 
	int disease_count, SimData* simData)
{
#ifndef __DEBUGED__
	// Same element
	if (simData->updatedCells->elementIdx[sim_cell_idx] == element_idx) {
		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			sim_cell_idx, mass, temperature, disease_idx, disease_count);
		return;
	}
#endif
	// Find a candidate cell which is not solid
	int     cellTarg = sim_cell_idx;
	uint8_t state    = GET_STATE(simData->updatedCells, cellTarg);
	if (state == 3) {
		int candidate_cells[4] = { cellTarg - 1, cellTarg + 1, cellTarg - simData->width , cellTarg + simData->width };
		for (int cell : candidate_cells) {
			if (GET_STATE(simData->updatedCells, cell) != 3) {
				cellTarg = cell;
				state    = GET_STATE(simData->updatedCells, cellTarg);
				break;
			}
		}
	}
#ifdef __DEBUGED__
	// Same element
	if (simData->updatedCells->elementIdx[cellTarg] == element_idx) {
		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			cellTarg, massDelta, temperature, disease_idx, disease_count);
		return;
	}
#endif

	uint16_t elemTarg = simData->updatedCells->elementIdx[cellTarg];
	if (state != 3 && gElements[elemTarg].id != simData->voidElementIdx) {
		// Vacuum
		if (state == 0) {
			simData->updatedCells->elementIdx[cellTarg] = element_idx;
			AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
				cellTarg, massDelta, temperature, disease_idx, disease_count);
			simData->simEvents->ChangeSubstance(simData, cellTarg);
			return;
		}
		// Try to move liquid/gas
		if ((state == 1 && DisplaceGas   (simData, simData->simEvents.get(), simData->updatedCells.get(), cellTarg, elemTarg)) ||
			(state == 2 && DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), cellTarg, elemTarg)))
		{
			simData->updatedCells->elementIdx[cellTarg] = element_idx;
			AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
				cellTarg, massDelta, temperature, disease_idx, disease_count);
			simData->simEvents->ChangeSubstance(simData, cellTarg);
			return;
		}
		// Target cell is gas. Check if liquid above gas
		int cellAbove = cellTarg + simData->width;
		if (state == 1 && GET_STATE(simData->updatedCells, cellAbove) == 2) {
			simData->updatedCells->SwapCells(cellAbove, cellTarg);
			simData->simEvents->ChangeSubstance(simData, cellAbove);
			simData->simEvents->ChangeSubstance(simData, cellTarg);
		}
		// Try to find same element in bigger range
		int candidate_cells[4] = { cellTarg + 1, cellTarg - 1, cellTarg + simData->width , cellTarg - simData->width };
		for (int cell : candidate_cells) {
			if (simData->updatedCells->elementIdx[cell] == element_idx) {
				AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(), 
					cell, massDelta, temperature, disease_idx, disease_count);
				return;
			}
		}
	}
	// Mass elimination
	if (massDelta >= simData->updatedCells->mass[cellTarg]) {
		float massAdd = massDelta - simData->updatedCells->mass[cellTarg];
		simData->updatedCells->elementIdx [cellTarg] = element_idx;
		simData->updatedCells->mass       [cellTarg] = 0.0;
		simData->updatedCells->temperature[cellTarg] = 0.0;

		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			cellTarg, massAdd, temperature, disease_idx, (int)(massAdd / massDelta * disease_count));
		ASSERT_TEMP(simData->updatedCells->mass[cellTarg], simData->updatedCells->temperature[cellTarg]);

		if (simData->updatedCells->mass[cellTarg] <= FLT_MIN) {
			//				simData->updatedCells->elementIdx[cellTarg] = simData->vacuumElementIdx;
			CellAccessor cellAccessor(simData->updatedCells.get(), cellTarg);
			simData->ClearCell(&cellAccessor);
		}
		simData->simEvents->ChangeSubstance(simData, cellTarg);
	}
	else {
		simData->updatedCells->mass[cellTarg] -= massDelta;
	}
}

void AddSolid(int base_sim_cell, uint16_t element_idx, float mass, float temperature, uint8_t disease_idx, 
	int disease_count, AddSolidMassSubType sub_type, SimData* simData)
{
	if (sub_type == DoVerticalDisplacement) {
		int   yOffset[2] = { -1, 0 };
		float massRemain = mass;
		for (int idx = 0; idx < 2 && massRemain > 0; idx++) {
			int cell = base_sim_cell + simData->width * yOffset[idx];
			if (simData->updatedCells->elementIdx[cell] != element_idx) 
				continue;
			float massMove = gElementPostProcessData[element_idx].maxMass - simData->updatedCells->mass[cell];
			if (massMove <= 0) 
				continue;

			massMove = MIN_F(massMove, massRemain);
			AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
				cell, massMove, temperature, disease_idx, (int)(massMove / mass * disease_count));
			massRemain -= massMove;
		}

		if (massRemain <= 0)
			return;
		
		int cellMax = simData->width * (simData->height - 1);
		for (int cell = base_sim_cell; cell < cellMax; cell += simData->width) {
			uint16_t elem  = simData->updatedCells->elementIdx[cell];
			uint8_t  state = gElements[elem].state & 3;
			if (state == 3) continue;
				
			if (gElements[elem].id != simData->voidElementIdx) {
				if (state == 1) DisplaceGas(   simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elem);
				if (state == 2) DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elem);

				simData->updatedCells->elementIdx [cell] = element_idx;
				simData->updatedCells->mass       [cell] = massRemain;
				simData->updatedCells->temperature[cell] = temperature;
				ASSERT_TEMP(simData->updatedCells->mass[cell], simData->updatedCells->temperature[cell]);
				simData->simEvents->ChangeSubstance(simData, cell);

				CellAccessor cellAccessor(simData->updatedCells.get(), cell);
#ifdef __DEBUGED__
				cellAccessor.SetDisease(disease_idx, (int)(massRemain / mass * disease_count));
#else
				cellAccessor.SetDisease(disease_idx, disease_count);
#endif
			}
			break;
		}
	}
	else if (sub_type == OnlyIfSameElement) {
		uint16_t elem_idx = simData->updatedCells->elementIdx[base_sim_cell];
		uint8_t  state    = gElements[elem_idx].state & 3;
		switch (state) {
		case 1: DisplaceGas(   simData, simData->simEvents.get(), simData->updatedCells.get(), base_sim_cell, elem_idx); break;
		case 2: DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), base_sim_cell, elem_idx); break;
		case 3: if (elem_idx != element_idx) return;
		}
		AddMassAndUpdateTemperature(simData->updatedCells.get(), simData->updatedCells.get(),
			base_sim_cell, mass, temperature, disease_idx, disease_count);
		ASSERT_TEMP(simData->updatedCells->mass[base_sim_cell], simData->updatedCells->temperature[base_sim_cell]);
		simData->simEvents->ChangeSubstance(simData, base_sim_cell);
	}
}

//void GameSyncFunction(BinaryBufferReader* reader, GameData* gameFrame, GameData* simFrame)
void GameSyncFunction(BinaryBufferReader* reader, GameData* simFrame)
{
	reader->ReadBytes(simFrame->width * simFrame->height, simFrame->visibleGrid.get());
}

//----- (ProcessFrame) -----------------------------------------------
void ProcessSetInsulation(SimFrameInfo* frame, SimData* simData)
{
	for (SetCellFloatValue msg : frame->setInsulationValues) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		simData->updatedCells->insulation[cell] = (uint8_t)(msg.value * 255);
	}
	frame->setInsulationValues.clear();
}

void ProcessSetStrength(SimFrameInfo* frame, SimData* simData)
{
	for (SetCellFloatValue msg : frame->setStrengthValues) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		simData->updatedCells->strengthInfo[cell] = (uint8_t)msg.value;
	}
	frame->setStrengthValues.clear();
}

void ProcessCellEnergyModifications(SimFrameInfo* frame, SimData* simData)
{
	for (CellEnergyModification msg : frame->cellEnergyModifications) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		ElementTemperatureData* data = &gElementTemperatureData[simData->updatedCells->elementIdx[cell]];
		if (!(data->state & 3))
			continue;
		if (simData->updatedCells->mass[cell] <= 0.001)
			continue;
		if (msg.maxTemperature <= 0.0)
			continue;

		float tempDelta = msg.kilojoules / (simData->updatedCells->mass[cell] * data->specificHeatCapacity);
		float tempFinal = MIN_F(simData->updatedCells->temperature[cell] + tempDelta, msg.maxTemperature);
		ASSERT_TEMP(simData->updatedCells->mass[cell], tempFinal);
		simData->updatedCells->temperature[cell] = CLAMP_F(tempFinal, SIM_MAX_TEMPERATURE, 1);
		DoStateTransition(simData, simData->simEvents.get(), cell, data);
	}
	frame->cellEnergyModifications.clear();
}

void ProcessCellProperties(SimFrameInfo* frame, SimData* simData)
{
	for (CellPropertiesChange msg : frame->setCellProperties) {
		int cell = CELL_G2S(msg.gameCell, simData->width);

		if (msg.set)
			simData->updatedCells->properties[cell] |= msg.properties;
		else
			simData->updatedCells->properties[cell] &= ~msg.properties;
		if (!simData->updatedCells->properties[cell])
			continue;
		if (simData->updatedCells->mass[cell] <= 0.0)
			continue;

		if (simData->updatedCells->properties[cell] & 1)
			DisplaceGas(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, simData->updatedCells->elementIdx[cell]);
		else if (simData->updatedCells->properties[cell] & 2)
			DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, simData->updatedCells->elementIdx[cell]);
	}
	frame->setCellProperties.clear();
}

void ProcessMassConsumption(SimFrameInfo* frame, SimData* simData)
{
	for (MassConsumption msg : frame->massConsumptionMessages) {
		std::vector<int>          visitedCells;
		std::queue<FloodFillInfo> nextCells;

		int cell = CELL_G2S(msg.gameCell, simData->width);
		int locX = cell % simData->width;
		int locY = cell / simData->width;

		ConsumedMassInfo info = { .simHandle = -1,.removedElemIdx = msg.elementIdx, .diseaseIdx = 0xFF, .mass = 0,.diseaseCount = 0 };
		FloodRemoved(simData, &msg.mass, info.removedElemIdx, &info, locX, locY, msg.radius, &visitedCells, &nextCells);

		if (msg.callbackIdx != -1){
			MassConsumedCallback callback = {
				.callbackIdx  = msg.callbackIdx,
				.elemIdx      = info.removedElemIdx,
				.diseaseIdx   = info.diseaseIdx,
				.mass         = info.mass,
				.temperature  = info.temperature,
				.diseaseCount = info.diseaseCount
			};
			simData->simEvents->massConsumedCallbacks.push_back(callback);
		}
	}
	frame->massConsumptionMessages.clear();
}

void ProcessMassEmission(SimFrameInfo* frame, SimData* simData)
{
	for (MassEmission msg : frame->massEmissionMessages) {
		int cell = CELL_G2S(msg.gameCell, simData->width);

		uint16_t elemCell  = simData->updatedCells->elementIdx[cell];
		uint16_t elemEmit  = msg.elementIdx;
		bool     displaced = false;

		if (elemCell != elemEmit && elemCell != simData->vacuumElementIdx) {
			// Different element, try displace
			uint8_t state = gElements[elemEmit].state & 3;
			if (state == 1) DisplaceGas   (simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elemCell);
			if (state == 2) DisplaceLiquid(simData, simData->simEvents.get(), simData->updatedCells.get(), cell, elemCell);
			if (!displaced) {
				if (msg.callbackIdx == -1)
					continue;
				MassEmittedCallback callBack = {
					.callbackIdx  = msg.callbackIdx,
					.elemIdx      = 0xFFFF,
					.emitted      = 0,
					.diseaseIdx   = 0xFF,
					.mass         = 0,
					.temperature  = 0,
					.diseaseCount = 0
				};
				simData->simEvents->massEmittedCallbacks.push_back(callBack);
				continue;
			}

			simData->updatedCells->temperature[cell] = CalculateFinalTemperature(
				simData->updatedCells->mass[cell], simData->updatedCells->temperature[cell], msg.mass, msg.temperature);
			simData->updatedCells->mass[cell] += msg.mass;
			if (msg.diseaseIdx != 0xFF)
				gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell, msg.diseaseIdx, msg.diseaseCount);
			
			if (elemCell != msg.elementIdx) {
				simData->updatedCells->elementIdx[cell] = msg.elementIdx;
				if (CELL_AVAILABLE(msg.gameCell, simData)) {
					SubstanceChangeInfo info = { .cellIdx = msg.gameCell, .oldElementIdx = (uint16_t)-1, .newElementIdx = (uint16_t)-1 };
					simData->simEvents->substanceChangeInfo.push_back(info);
				}
				simData->timers[cell].stableCellTicks |= 0x1F;
			}
			if (msg.callbackIdx != -1) {
				MassEmittedCallback callBack = {
					.callbackIdx  = msg.callbackIdx,
					.elemIdx      = msg.elementIdx,
					.emitted      = true,
					.diseaseIdx   = msg.diseaseIdx,
					.mass         = msg.mass,
					.temperature  = msg.temperature,
					.diseaseCount = msg.diseaseCount
				};
				simData->simEvents->massEmittedCallbacks.push_back(callBack);
			}
		}
	}
	frame->massEmissionMessages.clear();
}

void ProcessConsumeDisease(SimFrameInfo* frame, SimData* simData)
{
	for (ConsumeDisease msg : frame->consumeDisease) {
		int cell = CELL_G2S(msg.gameCell, simData->width);

		uint8_t diseIdx = -1;
		int     diseCnt = 0;
		if (simData->updatedCells->diseaseIdx[cell] != 0xFF) {
			diseCnt = std::max((int)(simData->updatedCells->diseaseCount[cell] * msg.percentToConsume + 0.5), msg.maxToConsume);
			simData->updatedCells->diseaseCount[cell] -= diseCnt;
			if (simData->updatedCells->diseaseCount[cell] <= 0) {
				simData->cells->diseaseIdx                   [cell] = -1;
				simData->cells->diseaseCount                 [cell] = 0;
				simData->cells->diseaseInfestationTickCount  [cell] = 0;
				simData->cells->diseaseGrowthAccumulatedError[cell] = 0;
			}
			diseIdx = simData->updatedCells->diseaseIdx[cell];
		}

		if (msg.callbackIdx != -1) {
			DiseaseConsumedCallback callBack = {
				.callbackIdx  = msg.callbackIdx,
				.diseaseIdx   = diseIdx,
				.diseaseCount = diseCnt
			};
			simData->simEvents->diseaseConsumedCallbacks.push_back(callBack);
		}
	}
	frame->consumeDisease.clear();
}

void ProcessCellDiseaseModifications(SimFrameInfo* frame, SimData* simData)
{
	for (CellDiseaseModification msg : frame->cellDiseaseModifications) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		if (msg.diseaseIdx == 0xFF) {
			simData->updatedCells->diseaseCount[cell] += msg.diseaseCount;
			if (simData->updatedCells->diseaseCount[cell] <= 0) {
				simData->cells->diseaseIdx                   [cell] = -1;
				simData->cells->diseaseCount                 [cell] = 0;
				simData->cells->diseaseInfestationTickCount  [cell] = 0;
				simData->cells->diseaseGrowthAccumulatedError[cell] = 0;
			}
		}
		else gDisease->AddDiseaseToCell(simData->updatedCells.get(), cell, msg.diseaseIdx, msg.diseaseCount);
	}
	frame->cellDiseaseModifications.clear();
}

void ProcessCellRadiationChanges(SimFrameInfo* frame, SimData* simData)
{
	for (CellRadiationModification msg : frame->cellRadiationModifications) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
#ifdef __DEBUGED__
		float radiaDelta = MAX_F(msg.radiationDelta, -simData->updatedCells->radiation[cell]);
		simData->updatedCells->radiation[cell] += radiaDelta;
#else
		float radiaDelta;
		if (simData->updatedCells->radiation[cell] + msg.radiationDelta > 0) {
			radiaDelta = msg.radiationDelta;
			simData->updatedCells->radiation[cell] += msg.radiationDelta;
		}
		else {
			radiaDelta = simData->updatedCells->radiation[cell];
			simData->updatedCells->radiation[cell] = 0;
		}
#endif
		if (msg.callbackIdx != -1) {
			RadiationConsumedCallback callback = {
				.callbackIdx = msg.callbackIdx,
				.gameCell    = msg.gameCell,
				.radiation   = radiaDelta
			};
			simData->simEvents->radiationConsumedCallbacks.push_back(callback);
		}
	}
	frame->cellRadiationModifications.clear();

	for (RadiationParamsModification msg : frame->radiationParamsModifications) {
		switch (msg.type) {
		case 0: simData->RADIATION_LINGER_RATE        = msg.value; break;
#ifdef __DEBUGED__
		case 1: simData->RADIATION_BASE_WEIGHT        = msg.value; break;
		case 2: simData->RADIATION_DENSITY_WEIGHT     = msg.value; break;
		case 3: simData->RADIATION_CONSTRUCTED_FACTOR = msg.value; break;
		case 4: simData->RADIATION_MAX_MASS           = msg.value; break;
#else
		case 2: simData->RADIATION_BASE_WEIGHT        = msg.value; break;
		case 3: simData->RADIATION_DENSITY_WEIGHT     = msg.value; break;
		case 4: simData->RADIATION_CONSTRUCTED_FACTOR = msg.value; break;
		case 5: simData->RADIATION_MAX_MASS           = msg.value; break;
#endif
		}
	}
	frame->radiationParamsModifications.clear();
}

void ProcessDigPoints(SimFrameInfo* frame, SimData* simData)
{
	for (DigPoint msg : frame->digPoints) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		// LOGGER_PRINT2("%s Cell %d->%d\n", __func__, msg.gameCell, cell);

		if (msg.callbackIdx != -1) {
			CallbackInfo info{ .callbackIdx = msg.callbackIdx };
			simData->simEvents->callbackInfo.push_back(info);
		}
		// Not Solid, skip
		if (GET_STATE(simData->updatedCells, cell) != 3)
			continue;

		if (CELL_AVAILABLE(msg.gameCell, simData)) {
			if (!msg.skipEvent && simData->updatedCells->mass[cell] > 0) {
				SpawnOreInfo info = {
					.cellIdx      = msg.gameCell,
					.elemIdx      = simData->updatedCells->elementIdx  [cell],
					.diseaseIdx   = simData->updatedCells->diseaseIdx  [cell],
					.mass         = simData->updatedCells->mass        [cell],
					.temperature  = simData->updatedCells->temperature [cell],
					.diseaseCount = simData->updatedCells->diseaseCount[cell]
				};
				simData->simEvents->digInfo.push_back(info);
			}
			SubstanceChangeInfo info = { .cellIdx = msg.gameCell, .oldElementIdx = (uint16_t)-1, .newElementIdx = (uint16_t)-1 };
			simData->simEvents->substanceChangeInfo.push_back(info);
		}
		CellAccessor cellAccessor = { simData->updatedCells.get(), cell };
		simData->ClearCell(&cellAccessor);
		simData->timers[cell].stableCellTicks |= 0x1F;
	}
	frame->digPoints.clear();
}

void ProcessCellModifications(SimFrameInfo* frame, SimData* simData)
{
	for (const CellModification& msg : frame->cellModifications) {
		int      cell          = CELL_G2S(msg.gameCell, simData->width);
		float    mass          = msg.mass;
		float    temperature   = msg.temperature;
		uint8_t  disease_idx   = msg.diseaseIdx;
		uint16_t elementIdx    = msg.elementIdx;
		int      disease_count = msg.diseaseCount;

		ASSERT_TEMP(mass, temperature);

		if (mass <= 0.0) {
			mass        = 0;
			temperature = 0;
		}
		else if (temperature <= 0.0) {
			temperature = 0.0;
		}
		switch (msg.replaceType) {
		case 0: // None
			switch (gElements[elementIdx].state & 3) {
			case 1:
				if (mass <= 0) {
					if (GET_STATE(simData->updatedCells, cell) != 1) break;

					simData->updatedCells->mass[cell] = MAX_F(mass + simData->updatedCells->mass[cell], 0);
					if (simData->updatedCells->mass[cell] <= FLT_MIN) {
						//CellAccessor cellAccessor = { simData->updatedCells.get(), cell };
						//simData->ClearCell(&cellAccessor);
						//simData->simEvents->ChangeSubstance(simData, cell);
						Evaporate(cell, simData, simData->simEvents.get());
					}
				}
				else if (mass > 0) {
					AddGas(cell, elementIdx, mass, temperature, disease_idx, disease_count, simData);
				}
				break;
			case 2:
				if (mass < 0) {
					if (GET_STATE(simData->updatedCells, cell) != 2) break;

					simData->updatedCells->mass[cell] = MAX_F(mass + simData->updatedCells->mass[cell], 0);
					if (simData->updatedCells->mass[cell] <= FLT_MIN) {
						Evaporate(cell, simData, simData->simEvents.get());
					}
				}
				else if (mass > 0) {
					AddLiquid(cell, elementIdx, mass, temperature, disease_idx, disease_count, simData);
				}
				break;
			case 3:
				if (mass < 0) {
					if (GET_STATE(simData->updatedCells, cell) != 3) break;

					simData->updatedCells->mass[cell] = MAX_F(mass + simData->updatedCells->mass[cell], 3);
					if (simData->updatedCells->mass[cell] <= FLT_MIN) {
						Evaporate(cell, simData, simData->simEvents.get());
					}
				}
				else if (mass > 0) {
					AddSolid(cell, elementIdx, mass, temperature, disease_idx, disease_count, 
						(AddSolidMassSubType)msg.addSubType, simData);
				}
				break;
			default:
				ASSERT_TEXT("Invalid replacement type");
				break;
			}
			break;
		case 1: // Replace
			ReplaceElement(cell, elementIdx, mass, temperature, disease_idx, disease_count, simData); 
			break;
		case 2: // ReplaceAndDisplace
			ReplaceAndDisplaceElement(cell, elementIdx, mass, temperature, disease_idx, disease_count, simData); 
			break;
		}
		
		if (msg.callbackIdx != -1) {
			CallbackInfo info{ .callbackIdx = msg.callbackIdx };
			simData->simEvents->callbackInfo.push_back(info);
		}
	}

	frame->cellModifications.clear();
}

void ProcessCellWorldZoneModification(SimFrameInfo* frame, SimData* simData)
{
	for (CellWorldZoneModification msg : frame->cellWorldZoneModifications) {
		int cell = CELL_G2S(msg.gameCell, simData->width);
		simData->worldZones[cell] = msg.zoneID;
	}
	frame->cellWorldZoneModifications.clear();
}

void AddBuildingHeatExchangeMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddBuildingHeatExchangeMessage msg : frame->buildingHeatExchangeMessages.adds) {
		Handle result = simData->buildingHeatExchange.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->buildingHeatExchangeMessages.adds.clear();
}

void ProcessBuildingHeatExchangeMessages(SimFrameInfo* frame, SimData* simData)
{
	for (ModifyBuildingHeatExchangeMessage msg : frame->buildingHeatExchangeMessages.modifies) {
		if (HANDLE_AVAILABLE(msg.handle, simData->buildingHeatExchange.handles.versions))
			simData->buildingHeatExchange.Modify(simData, &msg);
	}
	frame->buildingHeatExchangeMessages.modifies.clear();

	for (ModifyBuildingEnergyMessage msg : frame->modifyBuildingEnergyMessages) {
		if (HANDLE_AVAILABLE(msg.handle, simData->buildingHeatExchange.handles.versions))
			continue;
		int handleItm = simData->buildingHeatExchange.handles.items[msg.handle & 0xFFFFFF];
		BuildingHeatExchangeData* data = &simData->buildingHeatExchange.data[handleItm];
		float tempMin = MIN_F(msg.minTemperature, data->temperature);
		float tempMax = MAX_F(msg.maxTemperature, data->temperature);
		if (tempMin < 0 || tempMax > 10000)
			continue;

		int   cellCount = (data->simMax.y - data->simMin.y) * (data->simMax.x - data->simMin.x);
		float HcTotal   = cellCount * data->perCellHeatCapacity;
		float tempDelta = msg.deltaKJ / HcTotal;
		float tempFinal = CLAMP_F(data->temperature + tempDelta, tempMax, tempMin);
		if (tempFinal <= 0 || tempFinal >= 10000)
			ASSERT_TEXT("Invalid final temperature in ProcessBuildingHeatExchangeMessages");
		else
			data->temperature = tempFinal;
	}
	frame->modifyBuildingEnergyMessages.clear();
}

void RemoveBuildingHeatExchangeMessages(SimFrameInfo* frame, SimData* simData)
{
	for (RemoveBuildingHeatExchangeMessage msg : frame->buildingHeatExchangeMessages.removes) {
		simData->buildingHeatExchange.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->buildingHeatExchangeMessages.removes.clear();
}

void CombineB2BHeatExchangeMessages(SimFrameInfo* frame, SimData* simData)
{
	for (RegisterBuildingToBuildingHeatExchangeMessage msg : frame->buildingToBuildingHeatExchangeMessages.adds) {
		Handle result = simData->buildingToBuildingHeatExchange.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->buildingToBuildingHeatExchangeMessages.adds.clear();

	for (RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage msg : frame->buildingToBuildingHeatExchangeMessages.modifies) {
		if (HANDLE_AVAILABLE(msg.self_handle, simData->buildingToBuildingHeatExchange.handles.versions)) {
			simData->buildingToBuildingHeatExchange.Remove(simData, &msg);
		}
	}
	frame->buildingToBuildingHeatExchangeMessages.modifies.clear();

	for (AddInContactBuildingToBuildingToBuildingHeatExchangeMessage msg : frame->addBuildingInContactToBuildingToBuildingHeatExchangeMessages) {
		if (HANDLE_AVAILABLE(msg.self_handle, simData->buildingToBuildingHeatExchange.handles.versions)) {
			simData->buildingToBuildingHeatExchange.Add(simData, &msg);
		}
	}
	frame->addBuildingInContactToBuildingToBuildingHeatExchangeMessages.clear();

	for (RemoveBuildingToBuildingHeatExchangeMessage msg : frame->buildingToBuildingHeatExchangeMessages.removes) {
		simData->buildingToBuildingHeatExchange.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->buildingToBuildingHeatExchangeMessages.removes.clear();
}

void AddElementChunkMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddElementChunkMessage msg : frame->elementChunkMessages.adds) {
		Handle result = simData->elementChunk.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementChunkMessages.adds.clear();
}

void ProcessElementChunkMessages(SimFrameInfo* frame, SimData* simData)
{
	for (MoveElementChunkMessage msg : frame->moveElementChunkMessages) {
		if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
			simData->elementChunk.data[simData->elementChunk.handles.items[msg.handle & 0xFFFFFF]].cell = CELL_G2S(msg.gameCell, simData->width);
		}
	}
	frame->moveElementChunkMessages.clear();

	for (ModifyElementChunkMessage msg : frame->elementChunkMessages.modifies) {
		if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
			simData->elementChunk.Modify(simData, &msg);
		}
	}
	frame->elementChunkMessages.modifies.clear();

	for (ModifyElementChunkEnergyMessage msg : frame->modifyElementChunkEnergyMessages) {
		if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
			simData->elementChunk.ModifyEnergy(simData, &msg);
		}
	}
	frame->modifyElementChunkEnergyMessages.clear();

	for (ModifyElementChunkAdjusterMessage msg : frame->modifyElementChunkAdjusterMessages) {
		if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
			simData->elementChunk.ModifyAdjuster(simData, &msg);
		}
	}
	frame->modifyElementChunkAdjusterMessages.clear();
}

void RemoveElementChunkMessages(SimFrameInfo* frame, SimData* simData)
{
	for (RemoveElementChunkMessage msg : frame->elementChunkMessages.removes) {
		simData->elementChunk.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementChunkMessages.removes.clear();
}

void CombineElementConsumerMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddElementConsumerMessage msg : frame->elementConsumerMessages.adds) {
		Handle result = simData->elementConsumer.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementConsumerMessages.adds.clear();

	for (ModifyElementConsumerMessage msg : frame->elementConsumerMessages.modifies) {
		//if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
		simData->elementConsumer.Modify(simData, &msg);
		//}
	}
	frame->elementConsumerMessages.modifies.clear();

	for (RemoveElementConsumerMessage msg : frame->elementConsumerMessages.removes) {
		simData->elementConsumer.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementConsumerMessages.removes.clear();
}

void CombineElementEmitterMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddElementEmitterMessage msg : frame->elementEmitterMessages.adds) {
		Handle result = simData->elementEmitter.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementEmitterMessages.adds.clear();

	for (ModifyElementEmitterMessage msg : frame->elementEmitterMessages.modifies) {
		//if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
		simData->elementEmitter.Modify(simData, &msg);
		//}
	}
	frame->elementEmitterMessages.modifies.clear();

	for (RemoveElementEmitterMessage msg : frame->elementEmitterMessages.removes) {
		simData->elementEmitter.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->elementEmitterMessages.removes.clear();
}

void CombineRadiationEmitterMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddRadiationEmitterMessage msg : frame->radiationEmitterMessages.adds) {
		Handle result = simData->radiationEmitter.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->radiationEmitterMessages.adds.clear();

	for (ModifyRadiationEmitterMessage msg : frame->radiationEmitterMessages.modifies) {
		//if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
		simData->radiationEmitter.Modify(simData, &msg);
		//}
	}
	frame->radiationEmitterMessages.modifies.clear();

	for (RemoveRadiationEmitterMessage msg : frame->radiationEmitterMessages.removes) {
		simData->radiationEmitter.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->radiationEmitterMessages.removes.clear();
}

void CombineDiseaseEmitterMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddDiseaseEmitterMessage msg : frame->diseaseEmitterMessages.adds) {
		Handle result = simData->diseaseEmitter.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->diseaseEmitterMessages.adds.clear();

	for (ModifyDiseaseEmitterMessage msg : frame->diseaseEmitterMessages.modifies) {
		//if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
		simData->diseaseEmitter.Modify(simData, &msg);
		//}
	}
	frame->diseaseEmitterMessages.modifies.clear();

	for (RemoveDiseaseEmitterMessage msg : frame->diseaseEmitterMessages.removes) {
		simData->diseaseEmitter.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->diseaseEmitterMessages.removes.clear();
}

void CombineDiseaseConsumerMessages(SimFrameInfo* frame, SimData* simData)
{
	for (AddDiseaseConsumerMessage msg : frame->diseaseConsumerMessages.adds) {
		Handle result = simData->diseaseConsumer.Register(simData, &msg);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = result };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->diseaseConsumerMessages.adds.clear();

	for (ModifyDiseaseConsumerMessage msg : frame->diseaseConsumerMessages.modifies) {
		//if (HANDLE_AVAILABLE(msg.handle, simData->elementChunk.handles.versions)) {
		simData->diseaseConsumer.Modify(simData, &msg);
		//}
	}
	frame->diseaseConsumerMessages.modifies.clear();

	for (RemoveDiseaseConsumerMessage msg : frame->diseaseConsumerMessages.removes) {
		simData->diseaseConsumer.Unregister(msg.handle);
		if (msg.callbackIdx != -1) {
			ComponentStateChangedMessage info = { .callbackIdx = msg.callbackIdx, .simHandle = -1 };
			simData->simEvents->componentStateChangedMessages.push_back(info);
		}
	}
	frame->diseaseConsumerMessages.removes.clear();
}

//----- (Decompiled Main) --------------------------------------------
ParallelTaskQueue::ParallelTaskQueue(uint64_t threadCount)
{
	cnd_init(&this->mWorkCompleteCondition);
	cnd_init(&this->mTaskCondition);
	mtx_init(&this->mMutex, mtx_timed);

	this->mTasks        = std::queue<Task*>();
	this->mTaskCount    = 0;
	this->mWorkers      = std::vector<std::thread>();
	this->mShuttingDown = false;
	for (int i = 0; i < threadCount; i++) {
		this->mWorkers.push_back(std::thread(&ParallelTaskQueue::WorkerFunc,this));
	}
}

ParallelTaskQueue::~ParallelTaskQueue()
{
	if (!this->mShuttingDown) {
		mtx_lock(&this->mMutex);
		this->mShuttingDown = true;
		mtx_unlock(&this->mMutex);

		cnd_broadcast(&this->mTaskCondition);
		for (auto it = this->mWorkers.begin(); it != this->mWorkers.end(); it++)
			it->join();
	}
}

void ParallelTaskQueue::WorkerFunc()
{
	LOGGER_PRINT2("%s init\n", __func__);
	while (true) {
		mtx_lock(&this->mMutex);
		while (true) {
			if (this->mShuttingDown) {
				mtx_unlock(&this->mMutex);
				return;
			}

			if (this->mTasks.size()) break;

			cnd_wait(&this->mTaskCondition, &this->mMutex);
		}

		LOGGER_PRINT2("%s do task. Total %lld\n", __func__, this->mTasks.size());
		Task* task = this->mTasks.back();
		this->mTasks.pop();
		mtx_unlock(&this->mMutex);

		if (task == NULL) {
			ASSERT_TEXT("WorkerFunc Invalid Task");
			break;
		}
		task->InternalDoTask();

		mtx_lock(&this->mMutex);
		this->mTaskCount--;
		mtx_unlock(&this->mMutex);

		if (this->mTaskCount == 0)
			cnd_broadcast(&this->mWorkCompleteCondition);
		LOGGER_PRINT2("%s done task. mTaskCount %lld/%lld\n", __func__, this->mTaskCount, this->mTasks.size());
	}
}

SimFrameManager::SimFrameManager() 
{
	LOGGER_PRINT("%s\n", __func__);
	this->currentFrame    = NULL;
	this->framePool       = std::vector<SimFrameInfo*>();
	this->queuedFrames    = std::vector<SimFrameInfo*>();
	this->activeFrames    = std::vector<SimFrameInfo*>();
	this->processedFrames = std::vector<SimFrameInfo*>();
	this->activeRegions   = std::vector<ActiveRegion>();

	this->NewFrame();
}

int64_t SimFrameManager::BeginFrameProcessing()
{
	LOGGER_PRINT2("%s\n", __func__);
	mtx_lock(&this->frameMutex);

	this->framePool.insert(this->framePool.end(), this->processedFrames.begin(), this->processedFrames.end());
	this->processedFrames.clear();
	this->activeFrames.insert(this->activeFrames.end(), this->queuedFrames.begin(), this->queuedFrames.end());
	this->queuedFrames.clear();

	mtx_unlock(&this->frameMutex);
	this->numFramesProcessed = 0;
	LOGGER_PRINT2("%s done-%llu\n", __func__, this->activeFrames.size());
	return this->activeFrames.size();
}

void SimFrameManager::EndFrameProcessing()
{
	mtx_lock(&this->frameMutex);

	if (this->numFramesProcessed > 0) {
		auto it = this->activeFrames.begin();
		std::advance(it, this->numFramesProcessed);
		this->processedFrames.insert(this->processedFrames.end(), this->activeFrames.begin(), it);
		this->activeFrames.erase(this->activeFrames.begin(), it);
	}

	mtx_unlock(&this->frameMutex);
}

bool SimFrameManager::HandleMessage(Hashes::SimMessageHashes msg_id, int msg_length, char* msg)
{
	switch (msg_id) {
	case Hashes::SetDebugProperties:
		this->currentFrame->debugProperties = *(DebugProperties*)msg;
		return true;
	case Hashes::SimFrameManager_NewGameFrame:
		this->currentFrame->elapsedSeconds = *(float*)msg;
#ifdef __DEBUGED__
		mtx_lock(&this->regionMutex);
		LOGGER_PRINT2("%s lock\n", __func__);
#endif
		this->activeRegions.clear();
		for (int ptr = 0; ptr < msg_length; ptr += 0x1C) {
			ActiveRegion actRegion;
			actRegion.region.minimum.x = *(int*)(msg + ptr + 0x04);
			actRegion.region.minimum.y = *(int*)(msg + ptr + 0x08);
			actRegion.region.maximum.x = *(int*)(msg + ptr + 0x0C);
			actRegion.region.maximum.y = *(int*)(msg + ptr + 0x10);
			actRegion.currentSunlightIntensity        = *(float*)(msg + ptr + 0x14);
			actRegion.currentCosmicRadiationIntensity = *(float*)(msg + ptr + 0x18);
			this->activeRegions.push_back(actRegion);
		}
#ifdef __DEBUGED__
		LOGGER_PRINT2("%s unlock\n", __func__);
		mtx_unlock(&this->regionMutex);
#endif
		this->NewFrame();
		return true;
	case Hashes::Dig:
		this->currentFrame->digPoints.push_back(*(DigPoint*)msg);
		return true;
	case Hashes::ModifyCell:
		this->currentFrame->cellModifications.push_back(*(CellModification*)msg);
		return true;
	case Hashes::ModifyCellEnergy:
		this->currentFrame->cellEnergyModifications.push_back(*(CellEnergyModification*)msg);
		return true;
	case Hashes::SetInsulationValue:
		this->currentFrame->setInsulationValues.push_back(*(SetCellFloatValue*)msg);
		return true;
	case Hashes::SetStrengthValue:
		this->currentFrame->setStrengthValues.push_back(*(SetCellFloatValue*)msg);
		return true;
	case Hashes::ChangeCellProperties:
		this->currentFrame->setCellProperties.push_back(*(CellPropertiesChange*)msg);
		return true;
	case Hashes::ConsumeDisease:
		this->currentFrame->consumeDisease.push_back(*(ConsumeDisease*)msg);
		return true;
	case Hashes::CellDiseaseModification:
		this->currentFrame->cellDiseaseModifications.push_back(*(CellDiseaseModification*)msg);
		return true;
	case Hashes::CellRadiationModification:
		this->currentFrame->cellRadiationModifications.push_back(*(CellRadiationModification*)msg);
		return true;
	case Hashes::RadiationParamsModification:
		this->currentFrame->radiationParamsModifications.push_back(*(RadiationParamsModification*)msg);
		return true;
	case Hashes::ModifyCellWorldZone:
		this->currentFrame->cellWorldZoneModifications.push_back(*(CellWorldZoneModification*)msg);
		return true;
	case Hashes::MassConsumption:
		this->currentFrame->massConsumptionMessages.push_back(*(MassConsumption*)msg);
		return true;
	case Hashes::MassEmission:
		this->currentFrame->massEmissionMessages.push_back(*(MassEmission*)msg);
		return true;
	case Hashes::AddBuildingHeatExchange:
		this->currentFrame->buildingHeatExchangeMessages.adds.push_back(*(AddBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::ModifyBuildingHeatExchange:
		this->currentFrame->buildingHeatExchangeMessages.modifies.push_back(*(ModifyBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::RemoveBuildingHeatExchange:
		this->currentFrame->buildingHeatExchangeMessages.removes.push_back(*(RemoveBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::ModifyBuildingEnergy:
		this->currentFrame->modifyBuildingEnergyMessages.push_back(*(ModifyBuildingEnergyMessage*)msg);
		return true;
	case Hashes::AddBuildingToBuildingHeatExchange:
		this->currentFrame->buildingToBuildingHeatExchangeMessages.adds.push_back(*(RegisterBuildingToBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::RemoveBuildingInContactFromBuildingToBuildingHeatExchange:
		this->currentFrame->buildingToBuildingHeatExchangeMessages.modifies.push_back(*(RemoveBuildingInContactFromBuildingToBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::RemoveBuildingToBuildingHeatExchange:
		this->currentFrame->buildingToBuildingHeatExchangeMessages.removes.push_back(*(RemoveBuildingToBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::AddInContactBuildingToBuildingToBuildingHeatExchange:
		this->currentFrame->addBuildingInContactToBuildingToBuildingHeatExchangeMessages.push_back(*(AddInContactBuildingToBuildingToBuildingHeatExchangeMessage*)msg);
		return true;
	case Hashes::AddElementChunk:
		this->currentFrame->elementChunkMessages.adds.push_back(*(AddElementChunkMessage*)msg);
		return true;
	case Hashes::SetElementChunkData:
		this->currentFrame->elementChunkMessages.modifies.push_back(*(ModifyElementChunkMessage*)msg);
		return true;
	case Hashes::RemoveElementChunk:
		this->currentFrame->elementChunkMessages.removes.push_back(*(RemoveElementChunkMessage*)msg);
		return true;
	case Hashes::MoveElementChunk:
		this->currentFrame->moveElementChunkMessages.push_back(*(MoveElementChunkMessage*)msg);
		return true;
	case Hashes::ModifyElementChunkEnergy:
		this->currentFrame->modifyElementChunkEnergyMessages.push_back(*(ModifyElementChunkEnergyMessage*)msg);
		return true;
	case Hashes::ModifyChunkTemperatureAdjuster:
		this->currentFrame->modifyElementChunkAdjusterMessages.push_back(*(ModifyElementChunkAdjusterMessage*)msg);
		return true;
	case Hashes::AddElementConsumer:
		this->currentFrame->elementConsumerMessages.adds.push_back(*(AddElementConsumerMessage*)msg);
		return true;
	case Hashes::SetElementConsumerData:
		this->currentFrame->elementConsumerMessages.modifies.push_back(*(ModifyElementConsumerMessage*)msg);
		return true;
	case Hashes::RemoveElementConsumer:
		this->currentFrame->elementConsumerMessages.removes.push_back(*(RemoveElementConsumerMessage*)msg);
		return true;
	case Hashes::AddElementEmitter:
		this->currentFrame->elementEmitterMessages.adds.push_back(*(AddElementEmitterMessage*)msg);
		return true;
	case Hashes::ModifyElementEmitter:
		this->currentFrame->elementEmitterMessages.modifies.push_back(*(ModifyElementEmitterMessage*)msg);
		return true;
	case Hashes::RemoveElementEmitter:
		this->currentFrame->elementEmitterMessages.removes.push_back(*(RemoveElementEmitterMessage*)msg);
		return true;
	case Hashes::AddDiseaseEmitter:
		this->currentFrame->diseaseEmitterMessages.adds.push_back(*(AddDiseaseEmitterMessage*)msg);
		return true;
	case Hashes::ModifyDiseaseEmitter:
		this->currentFrame->diseaseEmitterMessages.modifies.push_back(*(ModifyDiseaseEmitterMessage*)msg);
		return true;
	case Hashes::RemoveDiseaseEmitter:
		this->currentFrame->diseaseEmitterMessages.removes.push_back(*(RemoveDiseaseEmitterMessage*)msg);
		return true;
	case Hashes::AddDiseaseConsumer:
		this->currentFrame->diseaseConsumerMessages.adds.push_back(*(AddDiseaseConsumerMessage*)msg);
		return true;
	case Hashes::ModifyDiseaseConsumer:
		this->currentFrame->diseaseConsumerMessages.modifies.push_back(*(ModifyDiseaseConsumerMessage*)msg);
		return true;
	case Hashes::RemoveDiseaseConsumer:
		this->currentFrame->diseaseConsumerMessages.removes.push_back(*(RemoveDiseaseConsumerMessage*)msg);
		return true;
	case Hashes::AddRadiationEmitter:
		this->currentFrame->radiationEmitterMessages.adds.push_back(*(AddRadiationEmitterMessage*)msg);
		return true;
	case Hashes::ModifyRadiationEmitter:
		this->currentFrame->radiationEmitterMessages.modifies.push_back(*(ModifyRadiationEmitterMessage*)msg);
		return true;
	case Hashes::RemoveRadiationEmitter:
		this->currentFrame->radiationEmitterMessages.removes.push_back(*(RemoveRadiationEmitterMessage*)msg);
		return true;
	case Hashes::SetLiquidCompression:
		for (int ptr = 0; ptr < msg_length; ptr += 0x6) {
			uint16_t element  = *(uint16_t*)(msg + ptr);
			float    compress = *(float   *)(msg + ptr + 0x2);
			LOGGER_PRINT("SetLiquidCompression %s:%f\n", gElementNames[element].c_str(), compress);
			gElementLiquidData     [element].compression = compress;
			gElementPostProcessData[element].compression = compress;
			gElementPostProcessData[element].maxCompression = MAX_F(powf(compress, 30), 1.01f);
		}
		return true;
	default:
		return false;
	}
}

void SimFrameManager::NewFrame()
{
	LOGGER_PRINT("%s\n", __func__);
	mtx_lock(&this->frameMutex);

	if (this->currentFrame) {
		std::vector<SimFrameInfo*>* p_framePool;
		if (this->queuedFrames.size() >= 0x1E)
			p_framePool = &this->framePool;
		else
			p_framePool = &this->queuedFrames;
		p_framePool->push_back(currentFrame);
	}

	if (this->framePool.size()) {
		this->currentFrame = this->framePool[this->framePool.size() - 1];
		this->framePool.pop_back();
	}
	else {
		this->currentFrame = new SimFrameInfo();// (SimFrameInfo*)calloc(1, sizeof(SimFrameInfo));
	}

	mtx_unlock(&this->frameMutex);
}

void SimFrameManager::ProcessFrame(SimFrameInfo* frame, SimData* simData)
{
	LOGGER_PRINT2("%s %lld/%lld\n", __func__, (uint64_t)frame, (uint64_t)simData);
	ProcessSetInsulation				(frame, simData); //Splited function
	ProcessSetStrength					(frame, simData); //Splited function
	ProcessCellEnergyModifications		(frame, simData);
	ProcessCellProperties				(frame, simData);
	ProcessMassConsumption				(frame, simData);
	ProcessMassEmission					(frame, simData);
	ProcessConsumeDisease				(frame, simData);
	ProcessCellDiseaseModifications		(frame, simData); //Splited function
	ProcessCellRadiationChanges			(frame, simData);
	ProcessDigPoints					(frame, simData);
	ProcessCellModifications			(frame, simData);
	ProcessCellWorldZoneModification	(frame, simData); //Splited function

	AddBuildingHeatExchangeMessages		(frame, simData); //Splited function
	ProcessBuildingHeatExchangeMessages	(frame, simData);
	RemoveBuildingHeatExchangeMessages	(frame, simData); //Splited function

	CombineB2BHeatExchangeMessages		(frame, simData); //Splited function

	AddElementChunkMessages				(frame, simData); //Splited function
	ProcessElementChunkMessages			(frame, simData);
	RemoveElementChunkMessages			(frame, simData); //Splited function

	CombineElementConsumerMessages		(frame, simData); //Splited function
	CombineElementEmitterMessages		(frame, simData); //Splited function
	CombineRadiationEmitterMessages		(frame, simData); //Splited function
	CombineDiseaseEmitterMessages		(frame, simData); //Splited function
	CombineDiseaseConsumerMessages		(frame, simData); //Splited function

	simData->debugProperties = frame->debugProperties;
	LOGGER_PRINT2("%s done\n", __func__);
}

float SimFrameManager::ProcessNextFrame(SimData* simData)
{
	SimFrameInfo* info = this->activeFrames[this->numFramesProcessed];
	this->ProcessFrame(info, simData);
#ifdef __DEBUGED__
	mtx_lock(&this->regionMutex);
	LOGGER_PRINT2("%s lock\n", __func__);
#endif
	if (&simData->activeRegions != &this->activeRegions)
		simData->activeRegions.assign(this->activeRegions.begin(), this->activeRegions.end());
#ifdef __DEBUGED__
	LOGGER_PRINT2("%s unlock\n", __func__);
	mtx_unlock(&this->regionMutex);
#endif
	this->numFramesProcessed++;
	return info->elapsedSeconds;
}

SimBase::SimBase() 
{
	//this->simFrameManager  = SimFrameManager();

	this->taskQueue        = new ParallelTaskQueue(1);
	this->temperatureTasks = std::vector<SimUpdateTask*>();
	this->simEvents        = std::vector<SimEvents*>();
	this->copyFlowTasks    = std::vector<UpdateCellTask<Vector4<float>>*>();
	this->copyToGameTasks  = std::vector<UpdateCellSOATask*>();
	this->elapsedSeconds   = 0;

	this->copyFlowTasks.push_back(new UpdateCellTask<Vector4<float>>());
	this->copyToGameTasks.push_back(new UpdateCellSOATask());
	this->InitializeUpdateTasks();
}

SimBase::~SimBase()
{
	// this->DestroyTasks();
	if (this->taskQueue) {
		this->taskQueue->~ParallelTaskQueue();
		delete(this->taskQueue);
	}

	for (auto& task : this->temperatureTasks) delete(task);
	for (auto& task : this->copyToGameTasks)  delete(task);
	for (auto& task : this->copyFlowTasks)    delete(task);
	for (auto& task : this->simEvents)        delete(task);
}

void SimBase::UpdateData(SimData* simData)
{
	LOGGER_PRINT2("%s, seed %d\n", __func__, simData->randomSeed);
	if (!simData->width)
		return;
	noffset      = -noffset;
	gGasDisplace = 1.01 + RAND_FLOAT_01(simData->randomSeed) * 0.99; // New_add. Make gas displacement stable.
	simData->cells->CopyFrom(simData->updatedCells.get());

	for (const ActiveRegion& activeRegion: simData->activeRegions) {
		LOGGER_PRINT2("\t%s: [%d-%d]x[%d-%d]\n", __func__, activeRegion.region.minimum.x, activeRegion.region.maximum.x, activeRegion.region.minimum.y, activeRegion.region.maximum.y);

		Region region = activeRegion.region;
		int yExtra = (region.maximum.y - region.minimum.y) % (int)this->taskQueue->mWorkers.size();
		int yInner = (region.maximum.y - region.minimum.y) / (int)this->taskQueue->mWorkers.size();

		for (int taskIdx = 0; taskIdx < this->temperatureTasks.size(); taskIdx++) {
			SimUpdateTask* task = this->temperatureTasks[taskIdx];

			region.maximum.y = yInner + region.minimum.y + (taskIdx < yExtra);
			task->simData	= simData;
			task->simEvents	= this->simEvents[taskIdx];
			task->start		= region.minimum;
			task->end		= region.maximum;
			
			mtx_lock(&this->taskQueue->mMutex);
			this->taskQueue->mTasks.push(task);
			this->taskQueue->mTaskCount++;
			mtx_unlock(&this->taskQueue->mMutex);
			if (cnd_signal(&this->taskQueue->mTaskCondition))
				throw (&this->taskQueue->mTaskCondition);
			region.minimum.y = region.maximum.y;
		}
		// Do temperatureTasks
		mtx_lock(&this->taskQueue->mMutex);
		while (this->taskQueue->mTaskCount)
			if (cnd_wait(&this->taskQueue->mWorkCompleteCondition, &this->taskQueue->mMutex))
				throw (&this->taskQueue->mMutex);
		mtx_unlock(&this->taskQueue->mMutex);
		LOGGER_PRINT2("\t%s: TemperatureTasks Done\n", __func__);

		region = activeRegion.region;
		simData->cells->CopyFrom(simData->updatedCells.get());
		// Process Gas Movement (not replace element)
		for (int yIdx = region.minimum.y; yIdx < region.maximum.y; yIdx++) {
			int cell     = simData->width * yIdx + region.minimum.x;
			int cell_end = simData->width * yIdx + region.maximum.x;

			// Replace cell<->cell_end
			if (noffset < 0) std::swap(cell, cell_end);

			// Traversal all cells
			for (; cell != cell_end; cell += noffset) {
				uint16_t elem_idx = simData->cells->elementIdx[cell];
				CellInfo cell_info = {
					.cell			= cell,
					.impermeable	= (simData->cells->properties[cell] & 1) == 1,
					.element_index	= elem_idx,
					.pressure		= &gElementPressureData[elem_idx],
					.element_state	=  gElementPressureData[elem_idx].state & 3u
				};
				if (cell_info.impermeable || simData->updatedCells->elementIdx[cell] != elem_idx  || cell_info.element_state > 1)
					continue;

				// Update Pressure with Next cell
				int      cell_next  = cell + noffset;
				uint16_t elem_next  = simData->cells->elementIdx[cell_next];
				uint8_t  state_next = gElementPressureData[elem_next].state & 3;

				if (   (simData->cells       ->properties[cell_next] & 1) == 0
					&&  simData->updatedCells->elementIdx[cell_next] == elem_next
					&&  state_next <= 1
					&& (elem_idx == elem_next || state_next != cell_info.element_state))
				{
					float massMove = UpdatePressure(simData,
						simData->simEvents.get(), simData->cells.get(), simData->updatedCells.get(), 
						cell, elem_idx, cell_info.element_state, gElementPressureData[elem_idx].flow, cell_next);

					simData->flow[cell     ].y += noffset * massMove;
					simData->flow[cell_next].x -= noffset * massMove;
				}

				// Update Pressure with Top cell
				int      cell_top  = cell + simData->width;
				uint16_t elem_top  = simData->cells->elementIdx[cell_top];
				uint8_t  state_top = gElementPressureData[elem_top].state & 3;
				if (    simData->updatedCells->elementIdx[cell    ] == elem_idx
					&& (simData->cells       ->properties[cell_top] & 1) == 0
					&&  simData->updatedCells->elementIdx[cell_top] == elem_top
					&&  state_top <= 1
					&& (elem_idx == elem_top || state_top != cell_info.element_state))
				{
					float massMove = UpdatePressure(simData, 
						simData->simEvents.get(), simData->cells.get(), simData->updatedCells.get(), 
						cell, elem_idx, cell_info.element_state, cell_info.pressure->flow, cell_top);
					
					simData->flow[cell    ].z += massMove;
					simData->flow[cell_top].w -= massMove;
				}

#ifndef __DEBUGED__
				if ((state_top | state_next) == 3) continue;
#else
				if (state_top == 3 || state_next == 3) continue;
#endif

				// Update Pressure with Diagonal cell
				int      cell_diag = cell + simData->width + noffset;
				uint16_t elem_diag = simData->cells->elementIdx[cell_diag];
				uint8_t  state_diag = gElementPressureData[elem_diag].state & 3;
				if (    simData->updatedCells->elementIdx[cell     ] == elem_idx
					&& (simData->cells       ->properties[cell_diag] & 1) == 0
					&&  simData->updatedCells->elementIdx[cell_diag] == elem_diag
					&&  state_diag <= 1
					&& (elem_idx == elem_diag || state_diag != cell_info.element_state))
				{
					UpdatePressure(simData, simData->simEvents.get(), simData->cells.get(), simData->updatedCells.get(),
						cell, elem_idx, cell_info.element_state, cell_info.pressure->flow, cell_diag);
				}
			}
		}
		LOGGER_PRINT2("\t%s: Gas Movement Done\n", __func__);

		int x_start = std::max(region.minimum.x,                   3);
		int x_end   = std::min(region.maximum.x, simData->width  - 3);
		int y_start = std::max(region.minimum.y,                   3);
		int y_end   = std::min(region.maximum.y, simData->height - 3);
		int x_step  = std::min(region.maximum.x - region.minimum.x, x_end - x_start);
		// Process Gas Pressure-based Element Replacement
		for (int yIdx = y_start; yIdx < y_end; yIdx++) {
			int cell     = simData->width * yIdx + x_start;
			int cell_end = cell + x_step;

			// Traversal all cells
			for (; cell != cell_end; cell++) {
				// Skip not Gas
				int16_t elem_idx = simData->cells->elementIdx[cell];
				if ((gElementPressureData[elem_idx].state & 3) != 1)
					continue;
				// Update Pressure Replace with 4 directions
				int offsets[4] = { -simData->width, -noffset, noffset, simData->width };
				for (int offset : offsets) {
					int     ncell  = cell + offset;
					int     nncell = cell + offset * 2;
					int16_t nelem  = simData->cells->elementIdx[ncell];

					if (nelem == elem_idx)								continue;
					if (simData->cells->properties[ncell] & 1)			continue;
					if ((gElementPressureData[nelem].state & 3) != 1)	continue;

					float massMove = DoGasPressureDisplacement(elem_idx, cell, ncell, nncell, simData);
					simData->flow[cell ].w -= massMove;
					simData->flow[ncell].z += massMove;
				}
			}
		}
		LOGGER_PRINT2("\t%s: Gas Replacement Done. [%d-%d]x[%d-%d]\n", __func__, x_start, x_end, y_start, y_end);

		simData->cells->CopyFrom(simData->updatedCells.get());
		this->ConsolidateEvents(simData);
		// Process Liquid Movement (not replace element)
		for (int yIdx = region.minimum.y; yIdx < region.maximum.y; yIdx++) {
			int cell     = simData->width * yIdx + region.minimum.x;
			int cell_end = simData->width * yIdx + region.maximum.x;

			// Traversal all cells
#ifdef __SIMDLL_PLUS__
			if (noffset < 0) std::swap(cell, cell_end);
			for (; cell != cell_end; cell += noffset) {
#else
			for (; cell != cell_end; cell++) {
#endif
				// Skip not liquid
				if ((gElementLiquidData[simData->cells->elementIdx[cell]].state & 3) != 2)
					continue;
				// void
				if (simData->cells->mass[cell] <= 0)
					continue;

				UpdateLiquid(simData, simData->simEvents.get(), cell);
			}
		}
		LOGGER_PRINT2("\t%s: Liquid Movement Done\n", __func__);

		// Process Liquid Pressure-based Element Replacement
		for (int yIdx = y_start; yIdx < y_end; yIdx++) {
			int cell     = simData->width * yIdx + region.minimum.x;
			int cell_end = simData->width * yIdx + region.maximum.x;
			// Traversal all cells
			for (; cell != cell_end; cell++) {
				// Skip not liquid
				int16_t elem_idx = simData->cells->elementIdx[cell];
				if ((gElementPressureData[elem_idx].state & 3) != 2)
					continue;
				// Update Pressure Replace with 3 directions
#ifdef __DEBUGED__
				if (simData->cells->elementIdx[cell - noffset] != elem_idx) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell - noffset, cell - noffset * 2, simData);
					simData->flow[cell        ].x += massMove;
					simData->flow[cell-noffset].y -= massMove;
				}
				if (simData->cells->elementIdx[cell + noffset] != elem_idx) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell + noffset, cell + noffset * 2, simData);
					simData->flow[cell        ].y += massMove;
					simData->flow[cell+noffset].x -= massMove;
				}
				if (simData->cells->elementIdx[cell + simData->width] != elem_idx && gElementPostProcessData[elem_idx].maxMass < simData->updatedCells->mass[cell]) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell + simData->width, cell + simData->width * 2, simData);
					simData->flow[cell               ].w += massMove;
					simData->flow[cell+simData->width].z -= massMove;
				}
#else
				if (simData->cells->elementIdx[cell - noffset] != elem_idx) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell - noffset, cell - noffset * 2, simData);
					simData->flow[cell  ].x += massMove;
					simData->flow[cell-1].y -= massMove;
				}
				if (simData->cells->elementIdx[cell + noffset] != elem_idx) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell + noffset, cell + noffset * 2, simData);
					simData->flow[cell  ].x += massMove;
					simData->flow[cell-1].y -= massMove;
				}
				if (simData->cells->elementIdx[cell + simData->width] != elem_idx && gElementPostProcessData[elem_idx].maxMass < simData->updatedCells->mass[cell]) {
					float massMove = DoLiquidPressureDisplacement(elem_idx, cell, cell + simData->width, cell + simData->width * 2, simData);
					simData->flow[cell               ].w += massMove;
					simData->flow[cell+simData->width].z -= massMove;
				}
#endif
			}
		}
		LOGGER_PRINT2("\t%s: Liquid Replacement Done\n", __func__);

		simData->cells->CopyFrom(simData->updatedCells.get());
		// Update Disease
		for (int yIdx = region.minimum.y; yIdx < region.maximum.y; yIdx++) {
			int cell     = simData->width * yIdx + region.minimum.x;
			int cell_end = simData->width * yIdx + region.maximum.x;
			// Traversal all cells
			for (; cell != cell_end; cell++) {
				if (simData->cells->diseaseIdx[cell] != 0xFF || simData->cells->diseaseIdx[cell + 1] != 0xFF)
					gDisease->UpdateCells(0.2f, simData, simData->simEvents.get(), cell, cell + 1);
				if (simData->cells->diseaseIdx[cell] != 0xFF || simData->cells->diseaseIdx[cell + simData->width] != 0xFF)
					gDisease->UpdateCells(0.2f, simData, simData->simEvents.get(), cell, cell + simData->width);
			}
		}
		LOGGER_PRINT2("\t%s: Update Disease Done\n", __func__);

		// Update Radiation
		if (simData->radiationEnabled) {
			uint8_t DiseaseIndex = gDisease->GetDiseaseIndex(0xD49F77D6);
			// Cosmic Radiation (Factor)
			simData->cosmicRadiationOcclusion.resize(simData->width * simData->height);
#ifdef __DEBUGED__
			float attenuation_factor = 1 - 1 / simData->RADIATION_LINGER_RATE;
			for (int yIdx = region.maximum.y - 1; yIdx >= region.minimum.y; yIdx--) {
#else
			for (int yIdx = region.maximum.y; yIdx >= region.minimum.y; yIdx--) {
#endif
				int cell     = simData->width * yIdx + region.minimum.x;
				int cell_end = simData->width * yIdx + region.maximum.x;
				for (; cell != cell_end; cell++) {
					// Factor of  up layer
					float inherit = 1.0f;
					if (yIdx + 1 < region.maximum.y)
						inherit = simData->cosmicRadiationOcclusion[cell + simData->width];
					// Decay by element of current layer
					float decay = gElementRadiationData[simData->updatedCells->elementIdx[cell]].factor;
					if (simData->updatedCells->properties[cell] & 0x80u) // ConstructedTile
						decay *= simData->RADIATION_CONSTRUCTED_FACTOR;
					else
						decay *= simData->updatedCells->mass[cell] / simData->RADIATION_MAX_MASS * simData->RADIATION_DENSITY_WEIGHT + simData->RADIATION_BASE_WEIGHT;
					inherit *= 1.0f - CLAMP_F(decay, 1, 0);
					if (inherit < 0.01) inherit = 0;
					simData->cosmicRadiationOcclusion[cell] = inherit;

					// Inherited from Pervious frame
#ifdef __DEBUGED__
					simData->updatedCells->radiation[cell] *= attenuation_factor;
#else
					if ((simData->updatedCells->radiation[cell] / simData->RADIATION_LINGER_RATE) == 0.0)
						simData->updatedCells->radiation[cell] = MAX_F(simData->updatedCells->radiation[cell] - 1.0f, 0.0f);
					else
						simData->updatedCells->radiation[cell] *= 1 - 1 / simData->RADIATION_LINGER_RATE;
#endif
				}
			}
			LOGGER_PRINT2("\t%s: Update Cosmic Radiation Done\n", __func__);

			Vector3<float> adjacentOffsets[25] = {
				{-2, -2, 0.1f }, {-1, -2, 0.15f}, {0, -2, 0.25f}, {1, -2, 0.15f}, {2, -2, 0.1f },
				{-2, -1, 0.15f}, {-1, -1, 0.5f }, {0, -1, 0.75f}, {1, -1, 0.5f }, {2, -1, 0.15f},
				{-2,  0, 0.25f}, {-1,  0, 0.75f}, {0,  0, 1    }, {1,  0, 0.75f}, {2,  0, 0.25f},
				{-2,  1, 0.15f}, {-1,  1, 0.5f }, {0,  1, 0.75f}, {1,  1, 0.5f }, {2,  1, 0.15f},
				{-2,  2, 0.1f }, {-1,  2, 0.15f}, {0,  2, 0.25f}, {1,  2, 0.15f}, {2,  2, 0.1f } };
			for (int yIdx = region.minimum.y; yIdx < region.maximum.y; yIdx++) {
				//LOGGER_PRINT2("\t%s: yIdx %d/%d\n", __func__, yIdx, region.maximum.y);
				for (int xIdx = region.minimum.x; xIdx < region.maximum.x; xIdx++) {
					int cell = simData->width * yIdx + xIdx;
					// Radiation diffusion to surrounding cells
					float massFactor = gElementRadiationData[simData->updatedCells->elementIdx[cell]].rads_per_1000;
					if (massFactor > 0) {
						for (Vector3<float> offset : adjacentOffsets) {
							int xIdx2 = xIdx + (int)offset.x;
							int yIdx2 = yIdx + (int)offset.y;
							if (xIdx2< region.minimum.x || xIdx2 >= region.maximum.x) continue;
							if (yIdx2< region.minimum.y || yIdx2 >= region.maximum.y) continue;

							int cella = simData->width * yIdx2 + xIdx2;
							simData->updatedCells->radiation[cella] += simData->updatedCells->mass[cell] * 0.001f * massFactor * offset.z;
						}
					}
					// Radioactive Contaminant (Germ)
					if (simData->updatedCells->diseaseIdx[cell] == DiseaseIndex) {
						simData->updatedCells->radiation[cell] += simData->updatedCells->diseaseCount[cell] * 0.001f;
					}
					// Cosmic Radiation (Value)
					if (simData->cosmicRadiationOcclusion[cell] > 0) {
						simData->updatedCells->radiation[cell] += activeRegion.currentCosmicRadiationIntensity / simData->RADIATION_LINGER_RATE * simData->cosmicRadiationOcclusion[cell];
					}
#ifdef __DEBUGED__
				}
			}
			for (int yIdx = region.minimum.y; yIdx < region.maximum.y; yIdx++) {
				for (int xIdx = region.minimum.x; xIdx < region.maximum.x; xIdx++) {
					int cell = simData->width * yIdx + xIdx;
#endif
					// Clamp
					simData->updatedCells->radiation[cell] = MIN_F(simData->updatedCells->radiation[cell], SIM_MAX_RADIATION);
					if (simData->updatedCells->radiation[cell] < 0.01)
						simData->updatedCells->radiation[cell] = 0;
				}
			}
		}
		LOGGER_PRINT2("\t%s: Update Radiation Done\n", __func__);

		simData->UpdateComponents(0.2f, &region);
		LOGGER_PRINT2("\t%s: UpdateComponents Done\n", __func__);

		x_start = std::max(x_start, 1);
		y_start = std::max(y_start, 1);
		x_end   = activeRegion.region.maximum.x + (x_end % 32 != 1);
		y_end   = activeRegion.region.maximum.y + (y_end % 32 != 1);
		x_end   = std::min(x_end, simData->width  - 1);
		y_end   = std::min(y_end, simData->height - 1);
		LOGGER_PRINT2("\t%s: Post Range[%d-%d]x[%d-%d]\n", __func__, x_start, x_end, y_start, y_end);

		// Space Biome Exposure
		if (simData->worldZones.get() && y_start < y_end && x_start < x_end) {
			for (int yIdx = y_start; yIdx < y_end; yIdx++) {
				int cell     = simData->width * yIdx + x_start;
				int cell_end = simData->width * yIdx + x_end;
				// Traversal all cells
				for (; cell != cell_end; cell++) {
					// Not space, skip
					if (simData->worldZones[cell] != 0xFF)
						continue;
					float lossRate = 0;
					switch (gElementPostProcessData[simData->updatedCells->elementIdx[cell]].state & 3) {
						case 1: lossRate = 1   ; break;
						case 2: lossRate = 1000; break;
						case 3: continue;
					}
					simData->updatedCells->mass[cell] = MAX_F(simData->updatedCells->mass[cell] - lossRate * 0.02f, 0);
					simData->flow[cell].x = -0.2f;
					simData->flow[cell].y =  0.2f;
					simData->flow[cell].z =  0.2f;
					simData->flow[cell].w = -0.2f;
				}
			}
		}
		LOGGER_PRINT2("\t%s: Space Biome Exposure Done\n", __func__);

		for (int yIdx = y_start; yIdx < y_end; yIdx++) {
			int cell     = simData->width * yIdx + x_start;
			int cell_end = simData->width * yIdx + x_end;
			// Traversal all cells
			for (; cell != cell_end; cell++) {
				PostProcessCell(simData, simData->simEvents.get(), cell);
			}
		}
		LOGGER_PRINT2("\t%s: PostProcessCell Done\n", __func__);

		gDisease->PostProcess(simData, simData->simEvents.get(), x_start, x_end, y_start, y_end);
		LOGGER_PRINT2("\t%s: PostProcess Done\n", __func__);

	}
	++simData->tickCount;
	LOGGER_PRINT2("%s done, tick %d\n", __func__, simData->tickCount);
}

void SimBase::InitializeUpdateTasks()
{
	for (int i = 0; i < this->taskQueue->mWorkers.size(); i++) {
		this->temperatureTasks.push_back(new SimUpdateTemperatureTask());
		this->simEvents.push_back(new SimEvents());
	}
}

void SimBase::ConsolidateEvents(SimData* simData)
{
#define VECTOR_MOVE(_src, _dest) for(const auto& _value : _src) {_dest.push_back(_value);} _src.clear()

	for (SimEvents* p_simEvent : this->simEvents) {
		VECTOR_MOVE(p_simEvent->substanceChangeInfo,			simData->simEvents->substanceChangeInfo);
		VECTOR_MOVE(p_simEvent->spawnOreInfo,					simData->simEvents->spawnOreInfo);
		VECTOR_MOVE(p_simEvent->spawnLiquidInfo,				simData->simEvents->spawnLiquidInfo);
		VECTOR_MOVE(p_simEvent->unstableCellInfo,				simData->simEvents->unstableCellInfo);
		VECTOR_MOVE(p_simEvent->elementChunkMeltedInfo,			simData->simEvents->elementChunkMeltedInfo);
		VECTOR_MOVE(p_simEvent->buildingMeltedInfo,				simData->simEvents->buildingMeltedInfo);
		VECTOR_MOVE(p_simEvent->buildingOverheatInfo,			simData->simEvents->buildingOverheatInfo);
		VECTOR_MOVE(p_simEvent->buildingNoLongerOverheatedInfo,	simData->simEvents->buildingNoLongerOverheatedInfo);
		VECTOR_MOVE(p_simEvent->cellMeltedInfo,					simData->simEvents->cellMeltedInfo);
		VECTOR_MOVE(p_simEvent->callbackInfo,					simData->simEvents->callbackInfo);
		VECTOR_MOVE(p_simEvent->worldDamageInfo,				simData->simEvents->worldDamageInfo);
		VECTOR_MOVE(p_simEvent->massConsumedCallbacks,			simData->simEvents->massConsumedCallbacks);
		VECTOR_MOVE(p_simEvent->massEmittedCallbacks,			simData->simEvents->massEmittedCallbacks);
		VECTOR_MOVE(p_simEvent->diseaseConsumedCallbacks,		simData->simEvents->diseaseConsumedCallbacks);
		VECTOR_MOVE(p_simEvent->spawnFXInfo,					simData->simEvents->spawnFXInfo);
		VECTOR_MOVE(p_simEvent->componentStateChangedMessages,	simData->simEvents->componentStateChangedMessages);
		VECTOR_MOVE(p_simEvent->radiationConsumedCallbacks,		simData->simEvents->radiationConsumedCallbacks);
		p_simEvent->digInfo.clear();
	}
}

void SimBase::CopyUpdatedCellsToCells(SimData* simData)
{
	simData->cells->CopyFrom(simData->updatedCells.get());
}

void SimBase::CopySimDataToGame(SimData* simData, GameData* new_game_data, const GameData* old_game_data, int num_frames_processed)
{
	LOGGER_PRINT2("%s. Processed frame:%d\n", __func__, num_frames_processed);
	int ySpace     = new_game_data->height / (int)this->copyFlowTasks.size();
	int yExtra     = new_game_data->height % this->copyFlowTasks.size();
	int srcYStart  = 1;
	int destYStart = 0;
	for (int i = 0; i < this->copyFlowTasks.size();i++) {
		LOGGER_PRINT2("%s. copyFlowTasks-%d/%llu\n", __func__, i+1, this->copyFlowTasks.size());
		UpdateCellTask<Vector4<float>>* task = this->copyFlowTasks[i];

		task->src			= simData->flow.get();
		task->dest			= new_game_data->flow.get();
		task->num_rows      = ySpace + (i < yExtra);
		task->src_y_start	= srcYStart;
		task->stride		= simData->width;
		task->dest_y_start	= destYStart;
		task->width			= new_game_data->width;

		mtx_lock(&this->taskQueue->mMutex);
		this->taskQueue->mTasks.push(task);
		this->taskQueue->mTaskCount++;
		mtx_unlock(&this->taskQueue->mMutex);
		if (cnd_signal(&this->taskQueue->mTaskCondition))
			throw (&this->taskQueue->mTaskCondition);

		srcYStart  += task->num_rows;
		destYStart += task->num_rows;
	}

	srcYStart  = 1;
	destYStart = 0;
	for (int i = 0; i < this->copyFlowTasks.size();i++) {
		LOGGER_PRINT2("%s. copyToGameTasks-%d/%llu\n", __func__, i+1, this->copyFlowTasks.size());
		UpdateCellSOATask* task = this->copyToGameTasks[i];

		task->src			= simData->updatedCells.get();
		task->dest			= new_game_data->cells.get();
		task->num_rows      = ySpace + (i < yExtra);
		task->src_y_start	= srcYStart;
		task->stride		= simData->width;
		task->dest_y_start	= destYStart;
		task->width			= new_game_data->width;

		mtx_lock(&this->taskQueue->mMutex);
		this->taskQueue->mTasks.push(task);
		this->taskQueue->mTaskCount++;
		mtx_unlock(&this->taskQueue->mMutex);
		if (cnd_signal(&this->taskQueue->mTaskCondition))
			throw (&this->taskQueue->mTaskCondition);

		srcYStart  += task->num_rows;
		destYStart += task->num_rows;
	}

	mtx_lock(&this->taskQueue->mMutex);
	while (this->taskQueue->mTaskCount)
		if (cnd_wait(&this->taskQueue->mWorkCompleteCondition, &this->taskQueue->mMutex))
			throw (&this->taskQueue->mMutex);
	mtx_unlock(&this->taskQueue->mMutex);

	new_game_data->elementChunkInfo.resize(simData->elementChunk.chunkInfo.size());
	new_game_data->elementChunkInfo.assign(simData->elementChunk.chunkInfo.begin(), simData->elementChunk.chunkInfo.end());
	for (int i = 0; i < simData->elementChunk.chunkInfo.size(); i++) {
		simData->elementChunk.chunkInfo[i].deltaKJ = 0.0;
	}

	std::swap(simData->elementConsumer.consumedMassInfo, new_game_data->consumedMassInfo);
	simData->elementConsumer.consumedMassInfo.clear();

	new_game_data->emittedMassInfo.resize(simData->elementEmitter.emittedMassInfo.size());
	new_game_data->emittedMassInfo.assign(simData->elementEmitter.emittedMassInfo.begin(), simData->elementEmitter.emittedMassInfo.end());
	for (int i = 0; i < simData->elementEmitter.emittedMassInfo.size(); i++) {
		simData->elementEmitter.emittedMassInfo[i].elemIdx = simData->vacuumElementIdx;
		simData->elementEmitter.emittedMassInfo[i].mass = 0;
	}

	new_game_data->buildingTemperatureInfo.resize(simData->buildingHeatExchange.temperatureInfo.size());
	new_game_data->buildingTemperatureInfo.assign(simData->buildingHeatExchange.temperatureInfo.begin(), simData->buildingHeatExchange.temperatureInfo.end());

	std::swap(simData->diseaseEmitter.emittedInfo,					new_game_data->diseaseEmittedInfo);
	std::swap(simData->diseaseConsumer.consumedInfo,				new_game_data->diseaseConsumedInfo);
	//std::swap(simData->simEvents->substanceChangeInfo,				new_game_data->substanceChangeInfo);
	std::swap(simData->simEvents->callbackInfo,						new_game_data->callbackInfo);
	std::swap(simData->simEvents->spawnLiquidInfo,					new_game_data->spawnFallingLiquidInfo);
	std::swap(simData->simEvents->unstableCellInfo,					new_game_data->unstableCellInfo);
	std::swap(simData->simEvents->elementChunkMeltedInfo,			new_game_data->elementChunkMeltedInfo);
	std::swap(simData->simEvents->buildingMeltedInfo,				new_game_data->buildingMeltedInfo);
	std::swap(simData->simEvents->buildingOverheatInfo,				new_game_data->buildingOverheatInfo);
	std::swap(simData->simEvents->buildingNoLongerOverheatedInfo,	new_game_data->buildingNoLongerOverheatedInfo);
	std::swap(simData->simEvents->cellMeltedInfo,					new_game_data->cellMeltedInfo);
	std::swap(simData->simEvents->worldDamageInfo,					new_game_data->worldDamageInfo);
	std::swap(simData->simEvents->massConsumedCallbacks,			new_game_data->massConsumedCallbacks);
	std::swap(simData->simEvents->radiationConsumedCallbacks,		new_game_data->radiationConsumedCallbacks);
	std::swap(simData->simEvents->massEmittedCallbacks,				new_game_data->massEmittedCallbacks);
	std::swap(simData->simEvents->diseaseConsumedCallbacks,			new_game_data->diseaseConsumedCallbacks);
	std::swap(simData->simEvents->spawnFXInfo,						new_game_data->spawnFXInfo);
	std::swap(simData->simEvents->componentStateChangedMessages,	new_game_data->componentStateChangedMessages);
	std::swap(simData->simEvents->digInfo,							new_game_data->digInfo);

	//spawnOreInfo: Combine infos in the same cellIdx and the same element
	new_game_data->spawnOreInfo.clear();
	std::sort(simData->simEvents->spawnOreInfo.begin(), simData->simEvents->spawnOreInfo.end(), 
		[](SpawnOreInfo a, SpawnOreInfo b){
			if (a.cellIdx != b.cellIdx) return a.cellIdx < b.cellIdx;
			if (a.elemIdx != b.elemIdx) return a.elemIdx < b.elemIdx;
			return a.diseaseCount < b.diseaseCount;
		});
	if (simData->simEvents->spawnOreInfo.size() > 0) {
		for (int i = 1; i < simData->simEvents->spawnOreInfo.size(); i++) {
			SpawnOreInfo* infoPre = &simData->simEvents->spawnOreInfo[i - 1];
			SpawnOreInfo* infoNxt = &simData->simEvents->spawnOreInfo[i];
			if (infoPre->cellIdx == infoNxt->cellIdx && infoPre->elemIdx == infoNxt->elemIdx) {
				infoNxt->temperature = CalculateFinalTemperature(infoPre->mass, infoPre->temperature, infoNxt->mass, infoNxt->temperature);
				infoNxt->mass += infoPre->mass;
			}
			else {
				new_game_data->spawnOreInfo.push_back(*infoPre);
			}
		}
		new_game_data->spawnOreInfo.push_back(simData->simEvents->spawnOreInfo.back());
	}

	//New Add, replace [std::swap(simData->simEvents->substanceChangeInfo, new_game_data->substanceChangeInfo);]
	//substanceChangeInfo: Combine infos in same cell
	new_game_data->substanceChangeInfo.clear();
	std::sort(simData->simEvents->substanceChangeInfo.begin(), simData->simEvents->substanceChangeInfo.end(),
		[](SubstanceChangeInfo a, SubstanceChangeInfo b) { return a.cellIdx < b.cellIdx; });
	if (simData->simEvents->substanceChangeInfo.size() > 0) {
		for (int i = 1; i < simData->simEvents->substanceChangeInfo.size(); i++) {
			int cellPre = simData->simEvents->substanceChangeInfo[i - 1].cellIdx;
			int cellNxt = simData->simEvents->substanceChangeInfo[i    ].cellIdx;
			//Add judgement: old_element != new_element
			if (cellPre != cellNxt && old_game_data->cells->elementIdx[cellPre] != new_game_data->cells->elementIdx[cellPre])
				new_game_data->substanceChangeInfo.push_back(simData->simEvents->substanceChangeInfo[i - 1]);
		}
		//Add judgement: old_element != new_element
		int cellEnd = simData->simEvents->substanceChangeInfo.back().cellIdx;
		if (old_game_data->cells->elementIdx[cellEnd] != new_game_data->cells->elementIdx[cellEnd])
			new_game_data->substanceChangeInfo.push_back(simData->simEvents->substanceChangeInfo.back());
	}

	simData->simEvents->substanceChangeInfo				.clear();
	simData->simEvents->spawnLiquidInfo					.clear();
	simData->simEvents->elementChunkMeltedInfo			.clear();
	simData->simEvents->buildingMeltedInfo				.clear();
	simData->simEvents->buildingOverheatInfo			.clear();
	simData->simEvents->buildingNoLongerOverheatedInfo	.clear();
	simData->simEvents->cellMeltedInfo					.clear();
	simData->simEvents->spawnOreInfo					.clear();
	simData->simEvents->spawnFXInfo						.clear();
	simData->simEvents->unstableCellInfo				.clear();
	simData->simEvents->callbackInfo					.clear();
	simData->simEvents->worldDamageInfo					.clear();
	simData->simEvents->massConsumedCallbacks			.clear();
	simData->simEvents->radiationConsumedCallbacks		.clear();
	simData->simEvents->massEmittedCallbacks			.clear();
	simData->simEvents->diseaseConsumedCallbacks		.clear();
	simData->simEvents->componentStateChangedMessages	.clear();
	simData->simEvents->digInfo							.clear();

	new_game_data->solidInfo							.clear();
	new_game_data->liquidChangeInfo						.clear();
	new_game_data->solidSubstanceChangeInfo				.clear();

	//solidInfo, liquidChangeInfo, solidSubstanceChangeInfo: Update base on substanceChangeInfo
	for (SubstanceChangeInfo& info : new_game_data->substanceChangeInfo) {
		int cellIdx = info.cellIdx;
		info.oldElementIdx = old_game_data->cells->elementIdx[cellIdx];
		info.newElementIdx = new_game_data->cells->elementIdx[cellIdx];
		uint8_t state_Old = gElementTemperatureData[info.oldElementIdx].state & 3;
		uint8_t state_New = gElementTemperatureData[info.newElementIdx].state & 3;

		//Remove judgement [old_element != new_element] to substanceChangeInfo
		if ((state_Old == 3) != (state_New == 3)) {
			SolidInfo infoAdd{ .cellIdx = cellIdx ,.solid = state_New == 3 };
			//LOGGER_PRINT2("%s. solidInfo:%d, %s->%s\n", __func__, cellIdx, gElementNames[info.oldElementIdx].c_str(), gElementNames[info.newElementIdx].c_str());
			new_game_data->solidInfo.push_back(infoAdd);
		}
		if (state_Old == 3 || state_New == 3) {
			SolidSubstanceChangeInfo infoAdd{ .cellIdx = cellIdx };
			//LOGGER_PRINT2("%s. solidSubstanceChangeInfo:%d, %s->%s\n", __func__, cellIdx, gElementNames[info.oldElementIdx].c_str(), gElementNames[info.newElementIdx].c_str());
			new_game_data->solidSubstanceChangeInfo.push_back(infoAdd);
		}
		if (state_Old == 2 || state_New == 2) {
			LiquidChangeInfo infoAdd{ .cellIdx = cellIdx };
			//LOGGER_PRINT2("%s. liquidChangeInfo::%d, %s->%s\n", __func__, cellIdx, gElementNames[info.oldElementIdx].c_str(), gElementNames[info.newElementIdx].c_str());
			new_game_data->liquidChangeInfo.push_back(infoAdd);
		}
	}

	for (int yIdx = 0; yIdx < new_game_data->height; yIdx++) {
		memmove(&new_game_data->accumulatedFlow[yIdx * new_game_data->width], 
			&simData->accumulatedFlow[(yIdx + 1) * simData->width + 1], 4 * new_game_data->width);
	}

	if (++tick_count >= 15)	{
		memset(simData->accumulatedFlow.get(), 0, 4 * simData->width * simData->height);
		tick_count = 0;
	}
	UpdateFlowTexture(simData, new_game_data);
	UpdateLiquidPropertyTexture(simData, new_game_data);
	UpdateExposedToSunPropertyTexture(simData, new_game_data);
	new_game_data->visibleGrid.swap(simData->visibleGrid); //new_game_data->swapVisibleGrid(&simData->visibleGrid);
	new_game_data->numFramesProcessed = num_frames_processed;

	LOGGER_PRINT2("%s done.\n", __func__);
}

void FrameSync::clear()
{
	this->mSimData.reset();
	this->mGameData.reset();
	this->mSimReadyToSwap = false;
}

void FrameSync::initGameData(int game_width, int game_height, int sim_width, int sim_height, CellSOA* updatedCells)
{
#define GameData_Move(_Key) memmove(&ptr[idx]->get()->cells->_Key[gameCell], &updatedCells->_Key[simCell], sizeof(ptr[idx]->get()->cells->_Key[0]) * game_width);

	LOGGER_PRINT2("%s %dx%d <- %dx%d\n", __func__, game_width, game_height, sim_width, sim_height);
	mtx_lock(&this->mMutex);
	if (this->mInitialized) ASSERT_TEXT("FrameSync Inited");

	std::unique_ptr<GameData>* ptr[2] = { &this->mSimData, &this->mGameData };
	for (int idx = 0; idx < 2; idx++) {
		*ptr[idx] = std::make_unique<GameData>(game_width, game_height);
		for (int locY = 0; locY < game_height; locY++) {
			int gameCell = game_width *  locY;
			int simCell  = sim_width  * (locY + 1) + 1;

			GameData_Move(elementIdx);
			GameData_Move(temperature);
			GameData_Move(mass);
			GameData_Move(properties);
			GameData_Move(insulation);
			GameData_Move(strengthInfo);
			GameData_Move(diseaseIdx);
			GameData_Move(diseaseCount);
			GameData_Move(diseaseInfestationTickCount);
			GameData_Move(diseaseGrowthAccumulatedError);
			GameData_Move(radiation);
		}
	}
	LOGGER_PRINT2("%s mSimData: %llu, mGameData: %llu\n", __func__, (uint64_t)this->mSimData.get(), (uint64_t)this->mGameData.get());
	this->mInitialized = true;
	mtx_unlock(&this->mMutex);
}

//void FrameSync::GameSync(void(*syncFunction)(GameData*, GameData*))
//{
//	mtx_lock(&this->mMutex);
//	while (!this->mSimReadyToSwap)
//		cnd_wait(&this->mGameCond, &this->mMutex);
//	syncFunction(this->mSimData.get(), this->mGameData.get());
//
//	mtx_lock(&this->mSimMutex);
//	std::swap(this->mSimData, this->mGameData);
//	this->mSimReadyToSwap = false;
//
//	mtx_unlock(&this->mSimMutex);
//	mtx_unlock(&this->mMutex);
//}
void FrameSync::GameSync_CleanUp()
{
	LOGGER_PRINT("%s\n", __func__);
	mtx_lock(&this->mMutex);
	while (!this->mSimReadyToSwap)
		cnd_wait(&this->mGameCond, &this->mMutex);
	gSim->mExitRequested = true; // Thread::Stop

	mtx_lock(&this->mSimMutex);
	std::swap(this->mSimData, this->mGameData);
	this->mSimReadyToSwap = false;

	mtx_unlock(&this->mSimMutex);
	mtx_unlock(&this->mMutex);

	cnd_signal(&this->mSimCond);
	LOGGER_PRINT("%s done\n", __func__);
}

void FrameSync::GameSync_PrepareGameData(BinaryBufferReader* reader)
{
	LOGGER_PRINT("%s init\n", __func__);
	mtx_lock(&this->mMutex);

	LOGGER_PRINT("%s 1 %s\n", __func__, this->mSimReadyToSwap ? "T" : "F");
	while (!this->mSimReadyToSwap)
		cnd_wait(&this->mGameCond, &this->mMutex);
	LOGGER_PRINT("%s 2\n", __func__);

	GameSyncFunction(reader, gSim->frameSync->mGameData.get());

	mtx_lock(&this->mSimMutex);
	std::swap(this->mSimData, this->mGameData);
	this->mSimReadyToSwap = false;

	mtx_unlock(&this->mSimMutex);
	mtx_unlock(&this->mMutex);

	cnd_signal(&this->mSimCond);
	LOGGER_PRINT("%s done\n", __func__);
}

void FrameSync::SimSync()
{
	LOGGER_PRINT2("%s init\n", __func__);
	mtx_lock(&this->mMutex);

	if (this->mSimReadyToSwap)
		ASSERT_TEXT("mSimReadyToSwap is true");
	this->mSimReadyToSwap = true;

	cnd_signal(&this->mGameCond);
	while (this->mSimReadyToSwap)
		cnd_wait(&this->mSimCond, &this->mMutex);

	mtx_unlock(&this->mMutex);
	LOGGER_PRINT2("%s done\n", __func__);
}

//----- (InternalDoTask) ---------------------------------------------
void SimUpdateTemperatureTask::InternalDoTask()
{
	LOGGER_PRINT2("SimUpdateTemperatureTask %s: [%d-%d]x[%d-%d]\n", __func__, this->start.x, this->end.x, this->start.y, this->end.y);
	for (int locY = this->start.y; locY < this->end.y; locY++) {
		for (int locX = this->start.x; locX < this->end.x; locX++) {
			int      cell = locX + simData->width * locY;
			uint16_t elem = simData->cells->elementIdx[cell];
			if (simData->cells->mass[cell] >= 0.001f) {
				int      cellR = cell + 1;
				int      cellT = cell + simData->width;
				uint16_t elemR = simData->cells->elementIdx[cellR];
				uint16_t elemT = simData->cells->elementIdx[cellT];
				if (locX < this->end.x - 1
					&& fabsf(simData->cells->temperature[cell] - simData->cells->temperature[cellR]) >= 1
					&& simData->cells->mass[cellR] >= 0.001f
					&& gElementTemperatureData[elem].thermalConductivity > 0
					&& gElementTemperatureData[elemR].thermalConductivity > 0
					&& !(gElementTemperatureData[elem].state & 0x10)
					&& !(gElementTemperatureData[elemR].state & 0x10))
				{
					UpdateTemperature(simData, simEvents, cell, cellR);
				}
				if (locY < this->end.y - 1
					&& fabsf(simData->cells->temperature[cell] - simData->cells->temperature[cellT]) >= 1
					&& simData->cells->mass[cellT] >= 0.001f
					&& gElementTemperatureData[elem].thermalConductivity > 0
					&& gElementTemperatureData[elemT].thermalConductivity > 0
					&& !(gElementTemperatureData[elem].state & 0x10)
					&& !(gElementTemperatureData[elemT].state & 0x10))
				{
					UpdateTemperature(simData, simEvents, cell, cellT);
				}
			}
			// Not TemperatureInsulated
			if (!(gElementTemperatureData[elem].state & 0x10)) {
				float biomeTemp = simData->biomeTemperature[cell];
				if (biomeTemp >= 0) {
					// biomeTemperatureLerpRate = 0.001
					simData->updatedCells->temperature[cell] += 
						(biomeTemp - simData->updatedCells->temperature[cell]) * simData->debugProperties.biomeTemperatureLerpRate;
					DoStateTransition(simData, simEvents, cell, &gElementTemperatureData[elem]);
				}
			}
		}
	}
	LOGGER_PRINT2("SimUpdateTemperatureTask %s done\n", __func__);
}

void UpdateCellSOATask::InternalDoTask()
{
	LOGGER_PRINT2("UpdateCellSOATask %s %lld->%lld\n", __func__, (uint64_t)this->src, (uint64_t)this->dest);
	const CellSOA* src = this->src;
	CellSOA* dest      = this->dest;

	int width          = this->width;
	int src_y_start    = this->src_y_start;
	int dest_y_start   = this->dest_y_start;
	for (int src_y_end = src_y_start + this->num_rows; src_y_start < src_y_end; src_y_start++, dest_y_start++) {
		uint64_t cell_src  = src_y_start  * this->stride + 1;
		uint64_t cell_dest = dest_y_start * width;
		memmove(&dest->elementIdx	[cell_dest], &src->elementIdx	[cell_src], width * 2);
		memmove(&dest->temperature	[cell_dest], &src->temperature	[cell_src], width * 4);
		memmove(&dest->mass			[cell_dest], &src->mass			[cell_src], width * 4);
		memmove(&dest->properties	[cell_dest], &src->properties	[cell_src], width);
		memmove(&dest->insulation	[cell_dest], &src->insulation	[cell_src], width);
		memmove(&dest->strengthInfo	[cell_dest], &src->strengthInfo	[cell_src], width);
		memmove(&dest->diseaseIdx	[cell_dest], &src->diseaseIdx	[cell_src], width);
		memmove(&dest->diseaseCount	[cell_dest], &src->diseaseCount	[cell_src], width * 4);
		memmove(&dest->radiation	[cell_dest], &src->radiation	[cell_src], width * 4);
		memmove(&dest->diseaseInfestationTickCount	[cell_dest], &src->diseaseInfestationTickCount	[cell_src], width);
		memmove(&dest->diseaseGrowthAccumulatedError[cell_dest], &src->diseaseGrowthAccumulatedError[cell_src], width * 4);
	}

	LOGGER_PRINT2("UpdateCellSOATask %s done\n", __func__);
}

template <typename T>
void UpdateCellTask<T>::InternalDoTask()
{
	LOGGER_PRINT2("UpdateCellTask %s %lld->%lld\n", __func__, (uint64_t)this->src, (uint64_t)this->dest);
	int src_y_start  = this->src_y_start;
	int dest_y_start = this->dest_y_start;
	for (int src_y_end = src_y_start + this->num_rows; src_y_start < src_y_end; src_y_start++, dest_y_start++) {
		uint64_t cell_src  = src_y_start * this->stride + 1;
		uint64_t cell_dest = dest_y_start * this->width;
		memmove(&this->dest[cell_dest], &this->src[cell_src], sizeof(T) * this->width);
	}
	LOGGER_PRINT2("UpdateCellTask %s done\n", __func__);
}
