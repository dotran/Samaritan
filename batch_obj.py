# TODO: 
#   1. Change config.txt
#   2. Run the flippin' Samaritan until it ends. 
#   3. Run again and again...

#algorithm_name: MOEAD
#test_problem: DTLZ3
#number_variable 10
#number_objective 3
#popSize: 105
#max_evaluation: 52500
#runtime_output: 1
#output_interval: 10
#run_index_begin: 1
#run_index_end: 11
#analyse: FUN GD IGD HV PLOT

import subprocess

methods=['MOEAD','NSGA3']  # ,'IBEA','NSGA2','HypE','SMSEMOA'
num_methods=len(methods)
problem_suite=['WFG']
num_suite=len(problem_suite)
problem_number=[[43,44,45,47,48]]
num_obj=[3,4]
#num_var=[2,4,8,16]
popsize=300
output_interval=10
index_begin=1
index_end=1

 #num_var*6000 

for m in range(num_methods):
    for ps in range(num_suite):
        for pn in range(len(problem_number[ps])):
            for on in range(len(num_obj)):
                if methods[m] == 'MOEAD' and num_obj[on] == 4:
                    popsize=286
                else:
                    popsize=300
                evaluation= popsize*300
                
#                if problem_number[ps][pn]==1:
#                    nv=num_obj[on]+4
#                else:
#                    nv=num_obj[on]+9
                nv=100
                f=open('config.txt','w')
                f.write('algo_name: '+methods[m]+'\n')
                f.write('tst_prob: '+problem_suite[ps]+str(problem_number[ps][pn])+'\n')
                f.write('parameters: '+str(num_obj[on]-1)+'\n')
                f.write('num_var: '+str(nv)+'\n')
                f.write('num_obj: '+str(num_obj[on])+'\n')
                f.write('popSize: '+str(popsize)+'\n')
                f.write('max_eval: '+str(evaluation)+'\n')
                f.write('runtime_output: 1\n')
                f.write('output_interval: '+str(output_interval)+'\n')
                f.write('ind_begin: '+str(index_begin)+'\n')
                f.write('ind_end: '+str(index_end)+'\n')
                f.write('analyse: FUN HV PLOT') 
                f.close()
                subprocess.run("Samaritan.exe")
            
f.close()
#algorithm_name: MOEAD
#test_problem: DTLZ3
#number_variable 10
#number_objective 3
#popSize: 105
#max_evaluation: 52500
#runtime_output: 1
#output_interval: 10
#run_index_begin: 1
#run_index_end: 11
#analyse: FUN GD IGD HV PLOT    
        
        
        
        




















