from paraview.simple import *
import re
import time
import os
from os import path

renderView = GetActiveViewOrCreate('RenderView')
source = GetActiveSource()

dir = path.join(path.dirname(source.FileName[0]), time.strftime('vid_%d%m%y_%H%M%S'))
if path.isdir(dir):
    os.rmdir(dir)
os.mkdir(dir)

dp = GetDisplayProperties(source, view=renderView)
(arr_type, name) = dp.ColorArrayName
prefix = re.match('([/a-zA-Z0-9]+)', name).group(1)

if arr_type == 'POINTS':
    arrays = list(source.PointArrayStatus)
else:
    arrays = list(source.CellArrayStatus)

i = arrays.index(name)
while i < len(arrays) and arrays[i].startswith(prefix):
    (arr_type, name) = dp.ColorArrayName
    print(arrays[i])
    dp.SetScalarColoring(arrays[i], servermanager.GetAssociationFromString(arr_type))
    UpdateScalarBars()

    SaveScreenshot(path.join(dir, arrays[i] + ".png"), renderView, ImageResolution=[1280, 1024])
    i += 1
