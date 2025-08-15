import time
def SDBMLower(s):
	import ctypes
	if s is None:
		return 0
			
	hashnum = 0
	for i in s:
		hashnum = ord(i.lower()) + (hashnum<<6) + \
			(hashnum<<16) - hashnum
		hashnum = hashnum & 0xffffffff
	return ctypes.c_int32(hashnum).value

print(SDBMLower("SetLiquidCompression"))
time.sleep(10)