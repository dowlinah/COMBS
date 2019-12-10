# @ shell = /bin/bash
#
# @ job_name = upwindFortran
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 00:30:00
#
# @ account_no = hpcf
#
# @ output               = $(job_name).$(schedd_host).$(jobid).o
# @ error                = $(job_name).$(schedd_host).$(jobid).e
# @ notification         = never
# @ class                = General
# @ network.MPI          = sn_all,not_shared,US
# @ task_affinity        = core(1)
# @ node                 = 1
# @ total_tasks          = 1
# @ queue

exe=/home/tolnaia/COMBS/benchmarks/fidibench/build/upwind/fortran/upwindFortran

# On the compute nodes /home is /gpfs_external/filesets/nesi/home
exe=`echo $exe | perl -ne "s#home#gpfs_external/filesets/nesi/home#;print;"`

cmd="time poe $exe 32 10"
echo "running..."
echo "$cmd"
$cmd
