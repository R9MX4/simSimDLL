using HarmonyLib;
using System.IO;
using System.Collections.Generic;
using STRINGS;
using System;
using System.Reflection;

namespace simSimDLL
{
    public class S_Text
    {
        public static class SHOWMORE
        {
            public static LocString LOCATION  = "Pos: X={0}, Y={1}";
            public static LocString DEBUGCELL = "GameCell: {0}, SimCell: {1}";
            public static LocString MOLAR     = "{0} Molar";
        }
    }
    public class Patch
    {
        [HarmonyPatch(typeof(Localization), "Initialize")]
        public class Localization_Initialize_Patch
        {
            private static void Postfix()
            {
                Localization.RegisterForTranslation(typeof(S_Text));
                string textpath = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), "string");
                Localization.Locale locale = Localization.GetLocale();

                string loc = locale == null ? "en" : locale.Code;
                string textpath2 = Path.Combine(textpath, loc + ".po");

                Dictionary<string, string> locDic = new Dictionary<string, string>();
                if (File.Exists(textpath2))
                {
                    locDic = Localization.LoadStringsFile(textpath2, false);
                    Localization.OverloadStrings(locDic);
                    Localization.GenerateStringsTemplate(typeof(S_Text), textpath);
                    Console.WriteLine("MOD-SimSimDLL_Text: Update language: " + loc);
                }
                else
                {
                    Localization.OverloadStrings(locDic);
                    Console.WriteLine("MOD-SimSimDLL_Text: !!! " + textpath2 + " not exist !!!");
                }
            }
        }
    }
}