using HarmonyLib;
using UnityEngine;
using static simSimDLL.S_Text;

namespace simSimDLL
{
    public class C_ShowMore
    {
        public class Toolbar_AddPos
        {
            [HarmonyPatch(typeof(SelectToolHoverTextCard), "UpdateHoverElements")]
            public class Get_Flag
            {
                private static bool Prefix(ref SelectToolHoverTextCard __instance)
                {
                    //init
                    if (DrawStyle == null)
                    {
                        DrawStyle = __instance.Styles_BodyText.Standard;
                        DrawIcon = __instance.iconDash;
                    }

                    DrawCell = Grid.PosToCell(Camera.main.ScreenToWorldPoint(KInputManager.GetMousePos()));
                    DrawFlag = Grid.IsVisible(DrawCell) && (int)Grid.WorldIdx[DrawCell] == ClusterManager.Instance.activeWorldId;

                    return true;
                }
            }

            [HarmonyPatch(typeof(HoverTextDrawer), "EndDrawing")]
            public class Add_Position
            {
                private static void Prefix(ref HoverTextDrawer __instance)
                {
                    if (DrawFlag)
                    {
                        DrawFlag = false;
                        __instance.BeginShadowBar(false);
                        __instance.DrawIcon(DrawIcon, 18);
                        int posX = DrawCell % Grid.WidthInCells;
                        int posY = DrawCell / Grid.WidthInCells;
                        int SimCell = posY * (Grid.WidthInCells + 2) + posX + Grid.WidthInCells + 3;

                        __instance.DrawText(string.Format(SHOWMORE.LOCATION, posX, posY), DrawStyle);
                        __instance.NewLine(18);
                        __instance.DrawIcon(DrawIcon, 18);
                        __instance.DrawText(string.Format(SHOWMORE.DEBUGCELL, DrawCell, SimCell), DrawStyle);

                        __instance.EndShadowBar();
                    }
                }
            }

            [HarmonyPatch(typeof(HoverTextHelper), "MassStringsReadOnly")]
            public class Get_Flag2
            {
                private static void Postfix(int cell)
                {
                    if (DrawFlag && Grid.Element[cell].IsGas && Grid.Element[cell].molarMass > 0)
                    {
                        DrawFlag2 = true;
                        molarVol = Grid.Mass[cell] / Grid.Element[cell].molarMass * 1000;
                    }
                }
            }


            [HarmonyPatch(typeof(HoverTextDrawer), "EndShadowBar")]
            public class Add_MolarVolume
            {
                private static void Prefix(ref HoverTextDrawer __instance)
                {
                    if (DrawFlag2 && molarVol > 0)
                    {
                        DrawFlag2 = false;
                        __instance.NewLine(26);
                        __instance.DrawIcon(DrawIcon, 18);
                        if (molarVol > 10)
                            __instance.DrawText(string.Format(SHOWMORE.MOLAR, molarVol.ToString("0.0")), DrawStyle);
                        else
                            __instance.DrawText(string.Format(SHOWMORE.MOLAR, molarVol.ToString("0.000")), DrawStyle);
                    }
                }
            }

            private static bool DrawFlag = false;
            private static bool DrawFlag2 = false;
            private static float molarVol = 0;
            private static int DrawCell = 0;
            private static TextStyleSetting DrawStyle = null;
            private static Sprite DrawIcon = null;
        }
    }
}
