This project contains following parts:

  SimDLLPlus      : Used for replacing origin SimDLL
  
  simSimDLL       : Link SimDLLPlus and main program
  
  simSimDLL_lib   : Text part
  
  cp.bat          : Combine the above parts together. Please configure %gamedir% to the correct location.
  
Compile environment:

  SimDLLPlus      : VS2022, c++, v143
  
  simSimDLL       : VS2022, NET Framework v4.7.2
  
Predefine option for compile SimDLLPlus (you cant find them in pch.h):

  __DEBUGED__     : Fix bugs like liquid duplication. Actually, for some reasons, you can't disable this option otherwise this mod will crash.
  
  __DEBUG_PRINT__ : Enable/Disable generate log.
  
  __SIMDLL_PLUS__ : Enable/Disable features metioned above.
