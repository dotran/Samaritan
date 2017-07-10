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

methods=['NSGA3','MOEAD']  # ,'IBEA','NSGA2','HypE','SMSEMOA'
num_methods=len(methods)
problem_suite=['WFG']
num_suite=len(problem_suite)
problem_number=[[1,2,3,4,5,6,7,8,9]]
num_obj=[3]
#num_var=[2,4,8,16]
popsize=105
generations = [200,300,400,500]
output_interval=10
index_begin=1
index_end=1

#num_var*6000 

for ig, g in enumerate(generations):
    print 'Generation: ' + str(g)
    for im,m in enumerate(methods):
        for ips,ps in enumerate(problem_suite):
            for ipn,pn in enumerate(problem_number[ips]):
                for iobj,obj in enumerate(num_obj):
                    if m == 'MOEAD':
                        popsize=105
                    else:
                        popsize=100
                    evaluation= popsize*g

#                    if problem_number[ps][pn]==1:
#                        d=num_obj[on]+4
#                    else:
#                        d=num_obj[on]+9
                    d=24
                    k=4
                    f=open('config.txt','w')
                    f.write('algo_name: '+m+'\n')
                    f.write('tst_prob: '+ps+str(pn)+'\n')
                    f.write('parameters: '+str(k)+'\n')
                    f.write('num_var: '+str(d)+'\n')
                    f.write('num_obj: '+str(obj)+'\n')
                    f.write('popSize: '+str(popsize)+'\n')
                    f.write('max_eval: '+str(evaluation)+'\n')
                    f.write('runtime_output: 1\n')
                    f.write('output_interval: '+str(output_interval)+'\n')
                    f.write('ind_begin: '+str(index_begin)+'\n')
                    f.write('ind_end: '+str(index_end)+'\n')
                    f.write('analyse: FUN HV PLOT') 
                    f.close()
                    subprocess.call(['Samaritan.exe'])
                    plt_string='out\\'+ps+str(pn)+'_M'+str(obj)+'_D'+str(d)+'\\'+m+'\\FUN1.out'
                    subprocess.call(['python','out\\s3d.py',plt_string, 'g'+str(g), 'out\\'])
    
                
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
        
        
        
        




















