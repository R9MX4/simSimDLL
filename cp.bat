set	gamedir=H:\SteamLibrary\steamapps\common\OxygenNotIncluded\OxygenNotIncluded_Data\Managed
set	moddir=%USERPROFILE%\Documents\Klei\OxygenNotIncluded\mods\dev
set	cppdir=%cd%\SimDLLPlus\x64\Debug
set	cspdir=%cd%\simSimDLL\bin\Debug

copy	%gamedir%\Assembly-CSharp.dll				%cd%
copy	%gamedir%\Assembly-CSharp-firstpass.dll		%cd%
copy	%gamedir%\0Harmony.dll						%cd%
copy	%gamedir%\UnityEngine.dll					%cd%
copy	%gamedir%\UnityEngine.CoreModule.dll		%cd%
copy	%gamedir%\System.dll						%cd%

copy	%cppdir%\SimDLLPlus.dll			%cd%\simSimDLL_lib\SimDLL.sim
copy	%cppdir%\SimDLLPlus.pdb			%cd%\simSimDLL_lib\SimDLLPlus.pdb
del	/s	/q	%moddir%\simSimDLL\*.*
xcopy	/s	%cd%\simSimDLL_lib\*.*		%moddir%\simSimDLL
copy	%cspdir%\simSimDLL.dll			%moddir%\simSimDLL
REM pause