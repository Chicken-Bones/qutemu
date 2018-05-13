from paraview.simple import *
from collections import namedtuple
import vtk

VarType = namedtuple('VarType', ['name', 'colormap', 'colorfunc'])


def colorAPD(name, ctf):
    ctf.ApplyPreset('Blue to Red Rainbow', False)
    ctf.RescaleTransferFunction(0.0, 350.0)
    ctf.AutomaticRescaleRangeMode = "Never"


def colorData(name, ctf):
    ctf.ApplyPreset('Cool to Warm', False)
    ctf.RescaleTransferFunction(-100.0, 50.0)
    ctf.AutomaticRescaleRangeMode = "Never"
    ctf.NanColor = [0.5, 0.5, 0.5]


def colorActivation(name, ctf):
    for src in GetSources().values():
        info = src.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation(name, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        if arrayInfo is not None:
            range = arrayInfo.GetComponentRange(0)

    ctf.ApplyPreset(v.colormap, False)
    #arr = [0.0, 0.0]
    #ctf.GetRange(arr)
    ctf.RescaleTransferFunction(max(range[0], range[1] - 500), range[1])
    ctf.AutomaticRescaleRangeMode = "Never"


var_types = [
    VarType('/Activation', 'Blue to Red Rainbow', colorActivation),
    VarType('/APD', 'Blue to Red Rainbow', colorAPD),
    VarType('/Peak', 'Blue to Red Rainbow', None),
    VarType('/Data', 'Cool to Warm', colorData),
    VarType('V @ ', 'Cool to Warm', colorData)
]


def getVar(dataset_name):
    return next((v for v in var_types if dataset_name.startswith(v.name)), None)


datasets = []
for src in GetSources().values():
    for dset in list(src.PointArrayStatus) + list(src.CellArrayStatus):
        v = getVar(dset)
        if v is not None:
            datasets += [dset]

print datasets

for dset in datasets:
    v = getVar(dset)
    ctf = GetColorTransferFunction(dset)

    if v.colorfunc is not None:
        v.colorfunc(dset, ctf)
    else:
        ctf.ApplyPreset(v.colormap, True)