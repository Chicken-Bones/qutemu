from sys import argv
from os import path
import re
import math
import vtk
import h5py
import numpy as np
import vtk.util.numpy_support as vtknp

# add_hdf5.py mesh.vtu dataset.h5 output.vtu


def read_vtu():
    print('reading: ' + argv[1])
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(argv[1])
    reader.Update()
    return reader.GetOutput()


def write_vtu(ugrid):
    out_path = argv[3] if len(argv) > 3 else path.splitext(argv[2])[0] + ".vtu"
    print('writing: ' + out_path)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(out_path)
    writer.SetInputData(ugrid)
    writer.Write()


def read_permutation():
    perm_path = path.join(path.dirname(argv[2]), "permutation.txt")
    print('permuting: ' + perm_path)
    if not path.isfile(perm_path):
        return None

    return np.genfromtxt(perm_path, skip_header=1, dtype=np.int32)[:, 1]


def read_h5():
    print('reading: ' + argv[2])
    return h5py.File(argv[2], 'r')


def read_interval():
    log_path = path.join(path.dirname(argv[2]), "log.txt")
    if path.isfile(log_path):
        with open(log_path, 'r') as log_file:
            m = re.search("interval: (\d+)ms", log_file.read())
            if m is not None:
                return float(m.group(1))

    return None


def get_namer(dataset_name, dataset_count):
    if dataset_name == '/Data':
        dataset_name = 'V'
        if dataset_count > 1:
            interval = read_interval()
            if interval is not None:
                for i in range(2, -1, -1):
                    fmt = "%."+str(i)+"f"
                    s = fmt % interval
                    if not s.endswith('0'):
                        break

                print('time interval: ' + s)
                time_max = fmt % (interval * (dataset_count - 1))
                print('time max: ' + time_max)
                w = len(time_max)

                return lambda i: 'V @ ' + (fmt % (i*interval)).zfill(w) + 'ms'

    if dataset_count == 1:
        return lambda i: dataset_name

    w = int(math.floor(math.log10(dataset_count)) + 1)
    return lambda i: dataset_name + '_' + str(i).zfill(w)


def add_datasets(ugrid, perm, h5file):
    pdata = ugrid.GetPointData()
    cdata = ugrid.GetCellData()

    for key in h5file.keys():
        dataset = h5file[key]

        if len(dataset.shape) >= 2 and np.prod(dataset.shape[1:]) == pdata.GetNumberOfTuples():
            target = pdata
        elif len(dataset.shape) >= 2 and np.prod(dataset.shape[1:]) == cdata.GetNumberOfTuples():
            target = cdata
        else:
            print(dataset.name + ' has wrong shape: ' + str(dataset.shape))
            continue

        print('converting: ' + dataset.name)
        dataset_count = dataset.shape[0]
        namer = get_namer(dataset.name, dataset_count)
        for i in range(0, dataset_count):
            subset = np.squeeze(dataset[i, ...])
            if perm is not None:
                subset = subset[perm]

            arr = vtknp.numpy_to_vtk(subset, deep=1, array_type=vtk.VTK_FLOAT)
            arr.SetName(namer(i))
            target.AddArray(arr)
            print('added: ' + arr.GetName())


ugrid = read_vtu()
add_datasets(ugrid, read_permutation(), read_h5())
write_vtu(ugrid)
