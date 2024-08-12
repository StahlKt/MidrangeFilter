#functions to submit jobs one by one with apply


#send off jobs one by one with sapply in script for parallelisation
#INPUT
#job.ids: integer, job.id as indicated by getJobTable(reg=registry)
#registry: filepath string to existing registry (R package batchtools)
#mem: integer indicating memory needed in MB
#wall.time: integer indicating wall time needed for one job in minutes
#partition.name: string indicating with partition on cluster to use (cluster specific)
#OUTPUT
#returns null, because it's only meant to execute the submitJob for each job individualls
#andivise to wrap the function in invisible(), which will prevent from showing the Null-return, 
#but not the  notification of jobs submitted
submit_jobs_wrap<- function(job.ids, registry, mem, wall.time, partition.name){
  
  job_ids<-data.table(job.id=job.ids, chunk=1)
  
  submitJobs(
    ids = job_ids,
    reg = registry,
    resources = list(
      ntasks = 1, ncpus = 1, memory = mem, 
      walltime = wall.time,
      partition = partition.name,
      chunks.as.arrayjobs = FALSE
    )
  )
  return()
}


#limits on medium: 
#gwdd[169-176] 	64 GB
#dmp[011-076] 128 GB 
#amp[001-092] 384 GB 
#walltime max on medium should be 24h



check_job_results<-function(check.id, registry){
  
  if(!file.exists(loadResult(check.id, reg=registry)[1])){
    return(check.id)
  }
  else if(!file.exists(loadResult(check.id, reg=registry)[2])){
    return(check.id)
  }
  return(NULL)
}
