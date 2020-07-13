#script to remove jobs files

import os
import shutil

noofjobs = 10
#for jobno in range(201,noofjobs):
#    job = "job"+str(jobno)+"d.root"
  #  os.remove(job)

os.chdir("PandoraCondor/work/")

for jobno in range(0,noofjobs):
    job = "job"+str(jobno)
    if (os.path.exists(job)):
        shutil.rmtree(job)
        #os.remove(job)
