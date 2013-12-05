import ctypes
#!gcc -shared cell.c -o cell.dll
ctt=ctypes.cdll.LoadLibrary('cell.dll')
