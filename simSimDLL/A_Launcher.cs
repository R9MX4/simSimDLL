using HarmonyLib;
using System;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using System.IO;

namespace simSimDLL
{
    public class A_PI
    {
        [DllImport("SimDLL.sim")]
        public static extern void ConduitTemperatureManager_Initialize();

        [DllImport("SimDLL.sim")]
        public static extern void ConduitTemperatureManager_Shutdown();

        [DllImport("SimDLL.sim")]
        public static extern int ConduitTemperatureManager_Add(float contents_temperature, float contents_mass, int contents_element_hash, int conduit_structure_temperature_handle, float conduit_heat_capacity, float conduit_thermal_conductivity, bool conduit_insulated);

        [DllImport("SimDLL.sim")]
        public static extern int ConduitTemperatureManager_Set(int handle, float contents_temperature, float contents_mass, int contents_element_hash);

        [DllImport("SimDLL.sim")]
        public static extern void ConduitTemperatureManager_Remove(int handle);

        [DllImport("SimDLL.sim")]
        public static extern IntPtr ConduitTemperatureManager_Update(float dt, IntPtr building_conductivity_data);

        [DllImport("SimDLL.sim")]
        public static extern void ConduitTemperatureManager_Clear();

        [DllImport("SimDLL.sim")]
        public static extern void SIM_Initialize(Sim.GAME_MessageHandler callback);

        [DllImport("SimDLL.sim")]
        public static extern void SIM_Shutdown();

        [DllImport("SimDLL.sim")]
        public unsafe static extern IntPtr SIM_HandleMessage(int sim_msg_id, int msg_length, byte* msg);

        [DllImport("SimDLL.sim")]
        public unsafe static extern IntPtr SIM_HandleMessages(int sim_msg_id, int msg_length, int msg_count, byte* msg);

        [DllImport("SimDLL.sim")]
        public unsafe static extern byte* SIM_BeginSave(int* size, int x, int y);

        [DllImport("SimDLL.sim")]
        public static extern void SIM_EndSave();

        [DllImport("SimDLL.sim")]
        public static extern void SIM_DebugCrash();

        [DllImport("SimDLL.sim")]
        public unsafe static extern char* SYSINFO_Acquire();

        [DllImport("SimDLL.sim")]
        public static extern void SYSINFO_Release();
    }

    public unsafe class N_ReDirect
    {
        static readonly bool DEBUG_FLAGS = true;
        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Initialize")]
        public class EA_ConduitTemperatureManager_Initialize
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Initialize");
                A_PI.ConduitTemperatureManager_Initialize();
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Initialize done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Shutdown")]
        public class EA_ConduitTemperatureManager_Shutdown
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Shutdown");
                A_PI.ConduitTemperatureManager_Shutdown();
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Shutdown done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Add")]
        public class EA_ConduitTemperatureManager_Add
        {
            private static bool Prefix(float contents_temperature, float contents_mass, int contents_element_hash, int conduit_structure_temperature_handle, float conduit_heat_capacity, float conduit_thermal_conductivity, bool conduit_insulated, ref int __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Add");
                __result = A_PI.ConduitTemperatureManager_Add(contents_temperature, contents_mass, contents_element_hash, conduit_structure_temperature_handle, conduit_heat_capacity, conduit_thermal_conductivity, conduit_insulated);
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Add done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Set")]
        public class EA_ConduitTemperatureManager_Set
        {
            private static bool Prefix(int handle, float contents_temperature, float contents_mass, int contents_element_hash, ref int __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Set");
                __result = A_PI.ConduitTemperatureManager_Set(handle, contents_temperature, contents_mass, contents_element_hash);
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Set done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Remove")]
        public class EA_ConduitTemperatureManager_Remove
        {
            private static bool Prefix(int handle)
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Remove");
                A_PI.ConduitTemperatureManager_Remove(handle);
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Remove done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Update")]
        public class EA_ConduitTemperatureManager_Update
        {
            private static bool Prefix(float dt, IntPtr building_conductivity_data, ref IntPtr __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Update");
                __result = A_PI.ConduitTemperatureManager_Update(dt, building_conductivity_data);
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Update done");
                return false;
            }
        }

        [HarmonyPatch(typeof(ConduitTemperatureManager), "ConduitTemperatureManager_Clear")]
        public class EA_ConduitTemperatureManager_Clear
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Clear");
                A_PI.ConduitTemperatureManager_Clear();
                if (DEBUG_FLAGS) Console.WriteLine("ConduitTemperatureManager_Clear done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_Initialize")]
        public class EA_SIM_Initialize
        {
            private static bool Prefix(Sim.GAME_MessageHandler callback)
            {
                if (DEBUG_FLAGS) Console.WriteLine("SIM_Initialize");
                A_PI.SIM_Initialize(callback);
                if (DEBUG_FLAGS) Console.WriteLine("SIM_Initialize done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_Shutdown")]
        public class EA_SIM_Shutdown
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("SIM_Shutdown");
                A_PI.SIM_Shutdown();
                if (DEBUG_FLAGS) Console.WriteLine("SIM_Shutdown done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_HandleMessage")]
        public class EA_SIM_HandleMessage
        {
            private static bool Prefix(int sim_msg_id, int msg_length, byte* msg, ref IntPtr __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine(string.Format("SIM_HandleMessage {0}: {1} {2}", (SimMessageHashes)(sim_msg_id), sim_msg_id, msg_length));
                __result = A_PI.SIM_HandleMessage(sim_msg_id, msg_length, msg);
                if (sim_msg_id == 1078620451)
                {
                    Sim.GameDataUpdate* ptr = (Sim.GameDataUpdate*)((void*)__result);
                    if (DEBUG_FLAGS) Console.WriteLine("SIM_HandleMessage Done, result " + __result + "-" + ptr->propertyTextureExposedToSunlight);
                }
                else
                {
                    if (DEBUG_FLAGS) Console.WriteLine("SIM_HandleMessage done, result " + __result);
                }
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_HandleMessages")]
        public class EA_SIM_HandleMessages
        {
            private static bool Prefix(int sim_msg_id, int msg_length, int msg_count, byte* msg, ref IntPtr __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine(string.Format("SIM_HandleMessages: {0} {1}", msg_count, msg_length));
                __result = A_PI.SIM_HandleMessages(sim_msg_id, msg_length, msg_count, msg);
                if (DEBUG_FLAGS) Console.WriteLine("SIM_HandleMessages done, result " + __result);
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_BeginSave")]
        public class EA_SIM_BeginSave
        {
            private static bool Prefix(int* size, int x, int y, ref byte* __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine(string.Format("SIM_BeginSave: {0}, {1}", x, y));
                __result = A_PI.SIM_BeginSave(size, x, y);
                if (DEBUG_FLAGS) Console.WriteLine(string.Format("SIM_BeginSave done: {0}", *size));
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_EndSave")]
        public class EA_SIM_EndSave
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("SIM_EndSave");
                A_PI.SIM_EndSave();
                if (DEBUG_FLAGS) Console.WriteLine("SIM_EndSave done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SIM_DebugCrash")]
        public class EA_SIM_DebugCrash
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("SIM_DebugCrash");
                A_PI.SIM_DebugCrash();
                if (DEBUG_FLAGS) Console.WriteLine("SIM_DebugCrash done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SYSINFO_Acquire")]
        public class EA_SYSINFO_Acquire
        {
            private static bool Prefix(ref char* __result)
            {
                if (DEBUG_FLAGS) Console.WriteLine("SYSINFO_Acquire");
                __result = A_PI.SYSINFO_Acquire();
                if (DEBUG_FLAGS) Console.WriteLine("SYSINFO_Acquire done");
                return false;
            }
        }

        [HarmonyPatch(typeof(Sim), "SYSINFO_Release")]
        public class EA_SYSINFO_Release
        {
            private static bool Prefix()
            {
                if (DEBUG_FLAGS) Console.WriteLine("SYSINFO_Release");
                A_PI.SYSINFO_Release();
                if (DEBUG_FLAGS) Console.WriteLine("SYSINFO_Release done");
                return false;
            }
        }
    }
    public unsafe class N_SimDLLPlus
    {
        [HarmonyPatch(typeof(SimMessages), "CreateSimElementsTable")]
        public class SetLiquidCompression
        {
            private static void Postfix(List<Element> elements)
            {
                MemoryStream memoryStream = new MemoryStream(6 * elements.Count);
                BinaryWriter binaryWriter = new BinaryWriter(memoryStream);
                int count = 0;
                for (ushort i = 0; i < elements.Count; i++)
                {
                    if ((elements[i].state & Element.State.Solid) != Element.State.Liquid)
                        continue;
                    if (elements[i].maxCompression == 1.01f)
                        continue;

                    binaryWriter.Write(i);
                    binaryWriter.Write(elements[i].maxCompression);
                    count++;
                }

                Console.WriteLine("SetLiquidCompression Total: " + count);
                if (count == 0)
                    return;

                byte[] buffer = memoryStream.GetBuffer();
                fixed (byte* msg = &buffer[0])
                {
                    Sim.SIM_HandleMessage(1396974888, count * 6, msg);
                }
            }
        }
    }
}