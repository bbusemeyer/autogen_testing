import cif2crystal
import runcrystal
import runqwalk
import copy
import job_control as jc
import os
import json

import veritas

element_list=[]
element_list.append(cif2crystal.Cif2Crystal())
element_list.append(runcrystal.RunCrystal(
  submitter=veritas.LocalVeritasCrystalSubmitter(
    nn=1,np=8,time="0:20:00",queue="batch")))
element_list.append(runcrystal.RunProperties(
  submitter=veritas.LocalVeritasPropertiesSubmitter(
    nn=1,np=1,time="1:00:00",queue="batch")))
element_list.append(runqwalk.Crystal2QWalk())
element_list.append(runqwalk.QWalkVarianceOptimize(
  submitter=veritas.LocalVeritasQwalkSubmitter(
    nn=1,np=8,time="0:10:00",queue="batch")))
#element_list.append(runqwalk.QWalkEnergyOptimize(
#  submitter=veritas.LocalVeritasQwalkSubmitter(
#    nn=1,np=8,time="0:10:00",queue="batch")))
#element_list.append(runqwalk.QWalkRunVMC(
#  submitter=veritas.LocalVeritasQwalkSubmitter(
#    nn=1,np=8,time="0:10:00",queue="batch")))
element_list.append(runqwalk.QWalkRunDMC(
  submitter=veritas.LocalVeritasQwalkSubmitter(
    nn=1,np=8,time="0:20:00",queue="batch")))
element_list.append(runqwalk.QWalkRunPostProcess(
  submitter=veritas.LocalVeritasQwalkSubmitter(
    nn=1,np=2,time="0:20:00",queue="batch")))

default_job=jc.default_job_record("si.cif")
default_job['dft']['kmesh'] = [2,2,2]
default_job['dft']['functional']['hybrid'] = 0
default_job['dft']['tolinteg'] = [6,6,6,6,12]
default_job['dft']['basis']=[0.4,1,1]
default_job['dft']['maxcycle'] = 100
default_job['dft']['fmixing'] = 80
default_job['dft']['edifftol'] = 6
default_job['dft']['broyden'] = [0.1,60,20]
default_job['qmc']['variance_optimize']['reltol']=0.1
default_job['qmc']['variance_optimize']['abstol']=10
default_job['qmc']['dmc']['save_trace'] = True
default_job['qmc']['dmc']['nblock']=5
default_job['qmc']['dmc']['target_error']=0.1
default_job['qmc']['dmc']['timestep']=[0.1,0.2]
default_job['total_spin'] = 0
idbase = "si_ag_"

# A demonstration of varying basis parameters.
count=1
results = []

# Simple run.
name = idbase+"simple"
job_record = copy.deepcopy(default_job)
job_record['control']['id']=name
job_record['qmc']['postprocess']['obdm'] = True
job_record['qmc']['postprocess']['basis'] = "../atomic.basis"
job_record['qmc']['postprocess']['orb'] = "../atomic.orb"
results.append(jc.execute(job_record,element_list))
count+=1

element_list[-2] = runqwalk.QWalkRunDMC(
  submitter=veritas.LocalVeritasBundleQwalkSubmitter(
    nn=1,np=8,time="0:20:00",queue="batch"))

# Restart and change something. Also check bundling.
name = idbase+"edit"
job_record = copy.deepcopy(default_job)
job_record['dft']['basis']=[0.6,1,1]
job_record['dft']['restart_from'] = "../%s/fort.9"%(idbase+"simple")
job_record['control']['id']=name
results.append(jc.execute(job_record,element_list))
count+=1

# Too-many cycles case.
name = idbase+"toomany"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['maxcycle'] = 10
job_record['dft']['edifftol'] = 10
results.append(jc.execute(job_record,element_list))
count+=1

name = idbase+"toomany_resumed"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['maxcycle'] = 10
job_record['dft']['edifftol'] = 10
job_record['dft']['resume_mode'] = 'stubborn'
results.append(jc.execute(job_record,element_list))
count+=1

# Reduce allowed run time.
element_list[1] = runcrystal.RunCrystal(
  submitter=veritas.LocalVeritasCrystalSubmitter(
    nn=1,np=8,time="0:00:20",queue="batch"
  ))

# Demonstrate dft killed-job error correction.
name = idbase+"killed"
job_record = copy.deepcopy(default_job)
job_record['control']['id']   = name
job_record['dft']['edifftol'] = 14
results.append(jc.execute(job_record,element_list))
count+=1

name = idbase+"killed_revived"
job_record = copy.deepcopy(default_job)
job_record['control']['id']      = name
job_record['dft']['resume_mode'] = 'optimistic'
job_record['dft']['edifftol'] = 14
results.append(jc.execute(job_record,element_list))
count+=1

# Reduce allowed run time.
element_list[1] = runcrystal.RunCrystal(
  submitter=veritas.LocalVeritasCrystalSubmitter(
    nn=1,np=8,time="0:20:00",queue="batch"
  ))

# Variance optimize takes two tries.
name = idbase+"varx2"
job_record = copy.deepcopy(default_job)
job_record['control']['id']=name
job_record['dft']['restart_from'] = "../si_ag_simple/fort.79"
job_record['qmc']['variance_optimize']['reltol']=0.001
results.append(jc.execute(job_record,element_list))
count+=1

print("Outputing results to test_results.json")
json.dump(results,open("test_results.json",'w'))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
