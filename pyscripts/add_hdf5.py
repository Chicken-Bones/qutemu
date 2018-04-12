from sys import argv
from os import path
import math
import vtk
import h5py
import numpy as np
import vtk.util.numpy_support as vtknp

#add_hdf5.py mesh.vtu fibrosis.txt output.vtu

print('reading: ' + argv[1])
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(argv[1])
reader.Update()
ugrid = reader.GetOutput()
pdata = ugrid.GetPointData()
num_points = pdata.GetNumberOfTuples()

perm_path = argv[4] if len(argv) > 4 else path.join(path.dirname(argv[2]), "permutation.txt")
print ('permuting: ' + perm_path)
if path.isfile(perm_path):
    perm = np.genfromtxt(perm_path, skip_header=1, dtype=np.int32)[:, 1]
else:
    perm = None

print('reading: ' + argv[2])
h5file = h5py.File(argv[2], 'r')
for key in h5file.keys():
    dataset = h5file[key]

    if len(dataset.shape) < 2 or np.prod(dataset.shape[1:]) != num_points:
        print(dataset.name + ' has wrong shape: ' + str(dataset.shape))
        continue

    print ('converting: ' + dataset.name)
    w = math.floor(math.log10(dataset.shape[0]))+1
    for i in range(0, dataset.shape[0]):
        subset = np.squeeze(dataset[i, ...])
        if perm is not None:
            subset = subset[perm]

        arr = vtknp.numpy_to_vtk(subset, deep=1, array_type=vtk.VTK_FLOAT)
        arr.SetName(dataset.name + '_' + str(i).zfill(w))
        pdata.AddArray(arr)
        print('added: ' + arr.GetName())

out_path = argv[3] if len(argv) > 3 else path.splitext(argv[2])[0] + ".vtu"
print('writing: ' + out_path)
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(out_path)
writer.SetInputData(ugrid)
writer.Write()