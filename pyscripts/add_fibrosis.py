import vtk
import sys

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(sys.argv[1])
reader.Update()
ugrid = reader.GetOutput()

fibrosis = vtk.vtkFloatArray()
fibrosis.SetName('fibrosis')
with open(sys.argv[2]) as f:
    for line in f.readlines()[1:]:
        fibrosis.InsertNextValue(float(line))

ugrid.GetCellData().AddArray(fibrosis)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(sys.argv[3])
writer.SetInputData(ugrid)
writer.Write()