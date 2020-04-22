import xacc

# Get access to the desired QPU and
# allocate some qubits to run on
qpu = xacc.getAccelerator('tnqvm', {'vqe-mode': True})
optimizer = xacc.getOptimizer('nlopt')
buffer = xacc.qalloc(4)

geom = '''
H  0.000000   0.0      0.0
H   0.0        0.0  .7474
'''
H = xacc.getObservable('pyscf', {'basis': 'sto-3g', 'geometry': geom})

adapt = xacc.getAlgorithm('adapt-vqe', {'accelerator': qpu,
                                  'optimizer': optimizer,
                                  'observable': H,
                                  'nElectrons': 2,
                                  'pool': "uccsd"})
# execute
adapt.execute(buffer)
