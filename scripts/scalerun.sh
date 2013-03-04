mpirun -np 9 ./mpisolve -INITCONDTYPE 1 -NUM 64 -A 0.2 -EPS 0.001 -SAVEWAVEFNCS 1
mpirun --hostfile ~/my-hostfile.txt -np 9 ./mpisolve -INITCONDTYPE 0 -NUM 128 -A 0.1 -EPS 0.0001 -SAVEWAVEFNCS 0
#mpirun --hostfile ~/my-hostfile.txt -np 9 ./mpisolve -INITCONDTYPE 0 -NUM 128 -A 0.25 -EPS 0.001 -SAVEWAVEFNCS 1
#mpirun --hostfile ~/my-hostfile.txt -np 33 ./mpisolve -INITCONDTYPE 0 -NUM 256 -A 0.05 -EPS 0.0008 -SAVEWAVEFNCS 1
#mpirun --hostfile ~/my-hostfile.txt -np 33 ./mpisolve -INITCONDTYPE 0 -NUM 512 -A 0.025 -EPS 0.0002 -SAVEWAVEFNCS 0
