from sys import argv
from os import path
import vtk
import numpy as np
import vtk.util.numpy_support as vtknp

#add_fibrosis.py mesh.vtu fibrosis.txt output.vtu

print('reading: ' + argv[1])
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(argv[1])
reader.Update()
ugrid = reader.GetOutput()

print('adding: ' + argv[2])
arr = np.genfromtxt(argv[2], skip_header=1)
fibrosis = vtknp.numpy_to_vtk(arr, array_type=vtk.VTK_FLOAT)
fibrosis.SetName('fibrosis')
ugrid.GetCellData().AddArray(fibrosis)

out_path = argv[3] if len(argv) > 3 else path.splitext(argv[2])[0] + ".vtu"
print('writing: ' + out_path)
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(out_path)
writer.SetInputData(ugrid)
writer.Write()