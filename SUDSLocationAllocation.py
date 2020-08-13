from gurobipy import *
import time
import os
path, filename = os.path.split(os.path.realpath(__file__))
sep=os.sep
#First include SETS
#bm1=10 la big m del tubo
#max_conections
#bm2=10000

def solve(N,O,D,T,S,E,demand,b,alpha,max_ar,min_ar,pot_wat,tub,budget,cu,evap,inf,beta1,beta2,dry,prob,dist,d_bound,bm1,max_conections,bm2,gapallow,sm1,maxtime):
    path, filename = os.path.split(os.path.realpath(__file__))
    sep=os.sep
    file2=open(path+sep+'results'+sep+'SUDS_db_'+str(d_bound)+'_mc_'+str(max_conections)+'_beta_'+str(int(beta1*100))+'_gapallow_'+str(gapallow)+'_budget_'+str(budget)+'.txt',"a")
    file2.write("Solution\n")
    #Create model
    creating_time=time.time()
    m=Model("SUDSLocationAllocation")
    #m.setParam("OutputFlag",0)
    m.setParam("LogFile","")
    #m.setParam("MIPGapAbs", gapallow)
    m.setParam("MIPGap", gapallow)
    m.setParam(GRB.Param.TimeLimit, maxtime)

    #Create variables
    x={}
    y={}
    for o in O:
        for t in T:
            if max_ar[o,t]>0:
                x[o,t] = m.addVar(vtype=GRB.BINARY, name='x(%d,%d)' % (o,t))
                y[o,t] = m.addVar(vtype=GRB.CONTINUOUS, name='y(%d,%d)' % (o,t))

    i={}
    teta={}
    z={}
    for o in O:
        for s in S:
            for e in E:
                i[o,s,e] = m.addVar(vtype=GRB.CONTINUOUS, name='i(%d,%d,%d)' % (o,s,e))
                teta[o,s,e] = m.addVar(vtype=GRB.BINARY, name='teta(%d,%d,%d)' % (o,s,e))
                z[o,s,e] = m.addVar(vtype=GRB.CONTINUOUS, name='z(%d,%d,%d)' % (o,s,e))

    q={}
    v={}
    for o in O:
        for n in N:
            if o!=n and dist[o,n]<=d_bound and tub[o,n]>0:
                v[o,n]= m.addVar(vtype=GRB.BINARY, name='v(%d,%d)' % (o,n))
                for s in S:
                    for e in E:
                        q[o,n,s,e]=m.addVar(vtype=GRB.CONTINUOUS, name='q(%d,%d,%d,%d)' % (o,n,s,e))

    w={}
    for d in D:
        for s in S:
            for e in E:
                if dry[s]==1:
                    w[d,s,e]=m.addVar(vtype=GRB.CONTINUOUS, name='w(%d,%d,%d)' % (d,s,e))

    r={}
    for o in O:
        for s in S:
            for e in E:
                if dry[s]==0:
                    r[o,s,e]=m.addVar(vtype=GRB.CONTINUOUS, name='r(%d,%d,%d)' % (o,s,e))

    #Create constraints
    #1. Guarantees that only one topology can be installed per feasible site
    for o in O:
        m.addConstr(quicksum(x[o,t] for t in T if max_ar[o,t]>0)<= 1,name='R1(%d)' % (o))

    #2. The selected area minimum requirement and maximum allowed
    for o in O:
        for t in T:
            if max_ar[o,t]>0:
                m.addConstr(y[o,t] <= max_ar[o,t]*x[o,t],name='R2_1(%d,%d)' % (o,t))
                m.addConstr(y[o,t] >= min_ar[o,t]*x[o,t],name='R2_2(%d,%d)' % (o,t))

    #3.Budget
    m.addConstr(quicksum(cu[t]*y[o,t] for t in T for o in O if max_ar[o,t]>0)+quicksum(tub[o,n]*v[o,n] for o in O for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)<=budget,name='R3')

    #4. volume between nodes is allowed only the conection exits
    for o in O:
        for n in N:
            if o!=n and dist[o,n]<=d_bound and tub[o,n]>0:
                m.addConstr(v[o,n]<=quicksum(x[o,t] for t in T if max_ar[o,t]>0),name='R4_1(%d,%d)' % (o,n))
                for s in S:
                    for e in E:
                        m.addConstr(q[o,n,s,e]<=bm1*v[o,n],name='R4_2(%d,%d,%d,%d)' % (o,n,s,e))

    for o in O:
        m.addConstr(quicksum(v[o,n] for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)<=max_conections,name='R4_4(%d)' % (o))
        for s in S:
            for e in E:
                m.addConstr(quicksum(q[o,n,s,e] for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)<=i[o,s,e],name='R4_3(%d,%d,%d)' % (o,s,e))

    #5 Inventory constraints
    #5.1 for the first wet period
    for o in O:
        for e in E:
            m.addConstr(i[o,1,e]==b[o,1,e]-r[o,1,e], name="R5_1(%d,1,%d)" % (o,e))

    for o in O:
        for s in S:
            for e in E:
                #5.2 for the wet periods
                if dry[s]==0 and s>1:
                    m.addConstr(i[o,s,e]==i[o,s-1,e]+b[o,s,e]-r[o,s,e],name="R5_2(%d,%d,%d)" % (o,s,e))
                    m.addConstr(b[o,s,e]>=r[o,s,e],name="R5_3(%d,%d,%d)" % (o,s,e))
                #5.3 for dry periods
                if dry[s]==1:
                    m.addConstr(i[o,s,e]==i[o,s-1,e]-quicksum(q[o,n,s,e] for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)+quicksum(q[n,o,s,e] for n in O if o!=n and dist[n,o]<=d_bound and tub[n,o]>0)-z[o,s,e],name="R5_4(%d,%d,%d)" %(o,s,e))
                #6. The inventory cannot exceed the capacity
                m.addConstr(i[o,s,e]<=quicksum(alpha[t]*y[o,t] for t in T if max_ar[o,t]>0),name="R6(%d,%d,%d)" %(o,s,e))
                #8 Linearization for evaporation
                m.addConstr(i[o,s,e]>=quicksum(y[o,t]*(evap[t]+inf[o]) for t in T if max_ar[o,t]>0)-bm2*teta[o,s,e],name="R8_1(%d,%d,%d)" %(o,s,e))
                m.addConstr(i[o,s,e]<=quicksum(y[o,t]*(evap[t]+inf[o]) for t in T if max_ar[o,t]>0)-sm1*teta[o,s,e]+(1-teta[o,s,e])*bm2,name="R8_2(%d,%d,%d)" %(o,s,e))
                if s>1:
                    m.addConstr(z[o,s,e]<=bm2*(1-teta[o,s-1,e]),name="R8_3(%d,%d,%d)" %(o,s,e))
                    m.addConstr(z[o,s,e]<=quicksum(y[o,t]*(evap[t]+inf[o]) for t in T if max_ar[o,t]>0),name="R8_4(%d,%d,%d)" %(o,s,e))
                    m.addConstr(z[o,s,e]>=quicksum(y[o,t]*(evap[t]+inf[o]) for t in T if max_ar[o,t]>0)-bm2*teta[o,s-1,e],name="R8_5(%d,%d,%d)" %(o,s,e))

    #7. The demand must be satisfied
    for d in D:
        for e in E:
            for s in S:
                if dry[s]==1:
                    m.addConstr(demand[d]==w[d,s,e]+quicksum(q[o,d,s,e] for o in O if o!=d and dist[o,d]<=d_bound and tub[o,d]>0),name="R7(%d,%d,%d)" %(d,s,e))


    def printSolution(which,tempo):
        if m.status == GRB.status.OPTIMAL:
            file2.write(str(which)+"\n")
            file2.write("Tiempo(s): "+str(tempo)+"\n")
            file2.write("Optimal solution\n")
            file2.write('Objective 1: '+str(sum(prob[e]*w[d,s,e].x for d in D for s in S for e in E if dry[s]==1))+"\n")
            file2.write('Objective 2: '+str(sum(prob[e]*r[o,s,e].x for o in O for s in S for e in E if dry[s]<1))+"\n")
            file2.write('Objective 3: '+str(sum(prob[e]*pot_wat*w[d,s,e].x for d in D for s in S for e in E if dry[s]==1)+sum(cu[t]*y[o,t].x for t in T for o in O if max_ar[o,t]>0)+sum(tub[o,n]*v[o,n].x for o in O for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0))+"\n")
            file2.write("Budget: "+str(sum(cu[t]*y[o,t].x for t in T for o in O if max_ar[o,t]>0)+sum(tub[o,n]*v[o,n].x for o in O for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0))+"\n")
            maximum=0
            where=-1
            for o in O:
                if sum(v[o,n].x for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)>maximum:
                    maximum=sum(v[o,n].x for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0)
                    where=o
            file2.write("Maximum conections: "+str(maximum) +" in: "+str(where)+"\n")
            file2.write("x \t y \n")
            for o in O:
                for t in T:
                    if max_ar[o,t]>0 and x[o,t].x>0.9:
                        file2.write(str(x[o,t].VarName)+": "+str(x[o,t].x)+"\t"+str(y[o,t].VarName)+": "+str(y[o,t].x)+"\n")
            file2.write("\n")
            file2.write("i \t teta \t z\n")
            for o in O:
                for s in S:
                    for e in E:
                        if i[o,s,e].x>0.00001 or z[o,s,e].x>0.00001 or teta[o,s,e].x>0.9:
                            file2.write(str(i[o,s,e].VarName)+": "+str(i[o,s,e].x)+"\t"+str(teta[o,s,e].VarName)+": "+str(teta[o,s,e].x)+"\t"+str(z[o,s,e].VarName)+": "+str(z[o,s,e].x)+"\n")
            file2.write("\n")
            file2.write("v \n")
            for o in O:
                for n in N:
                    if o!=n and dist[o,n]<=d_bound and tub[o,n]>0 and v[o,n].x>0.9:
                        file2.write(str(v[o,n].VarName)+": "+str(v[o,n].x)+"\n")
            file2.write("\n")
            file2.write("q \n")
            for o in O:
                for n in N:
                    if o!=n and dist[o,n]<=d_bound and tub[o,n]>0:
                        for s in S:
                            for e in E:
                                if q[o,n,s,e].x>0.00001:
                                    file2.write(str(q[o,n,s,e].VarName)+": "+str(q[o,n,s,e].x)+"\n")
            file2.write("\n")
            file2.write("w \n")
            for d in D:
                for s in S:
                    for e in E:
                        if dry[s]==1:
                            if w[d,s,e].x>0.00001:
                                file2.write(str(w[d,s,e].VarName)+": "+str(w[d,s,e].x)+"\n")

            file2.write("\n")
            file2.write("r \n")
            for o in O:
                for s in S:
                    for e in E:
                        if dry[s]==0:
                            if r[o,s,e].x>0.00001:
                                file2.write(str(r[o,s,e].VarName)+": "+str(r[o,s,e].x)+"\n")
            file2.write("----------------------------------------------------------------------------------------------------------\n")
            for e in E:
                file2.write("w: "+str(e)+"->"+str(sum(w[d,s,e].x for d in D for s in S if dry[s]==1)))

            for e in E:
                file2.write("r: "+str(e)+"->"+str(sum(r[o,s,e].x for o in O for s in S if dry[s]==0)))
        else:
            file2.write("No solution\n")
            file2.write("----------------------------------------------------------------------------------------------------------\n")
    #Objective function
    #Potable water
    print("Potable water")
    m.setObjective(quicksum(prob[e]*w[d,s,e] for d in D for s in S for e in E if dry[s]==1),GRB.MINIMIZE)
    m.update()
    creating_time=time.time()-creating_time
    run_time1=time.time()
    #m.write('John_SUDS_db_'+str(d_bound)+'_mc_'+str(max_conections)+'_beta_'+str(int(beta1*100))+'_gapallow_'+str(gapallow)+".lp")
    m.optimize()

    run_time1=time.time()-run_time1
    printSolution("minimizing potable water",run_time1)


    for e in E:
        m.addConstr(quicksum(w[d,s,e] for d in D for s in S if dry[s]==1)<=quicksum(w[d,s,e].x for d in D for s in S if dry[s]==1)*(1+beta1))

    m.setObjective(quicksum(prob[e]*r[o,s,e] for o in O for s in S for e in E if dry[s]<1),GRB.MINIMIZE)
    m.update()
    run_time2=time.time()
    print("Runoff")
    m.optimize()
    run_time2=time.time()-run_time2
    printSolution("minimizing runoff",run_time2)
    for e in E:
        m.addConstr(quicksum(r[o,s,e] for o in O for s in S if dry[s]<1)<=quicksum(r[o,s,e].x for o in O for s in S if dry[s]<1)*(1+beta2))

    m.setObjective(quicksum(prob[e]*pot_wat*w[d,s,e] for d in D for s in S for e in E if dry[s]==1)+quicksum(cu[t]*y[o,t] for t in T for o in O if max_ar[o,t]>0)+quicksum(tub[o,n]*v[o,n] for o in O for n in N if o!=n and dist[o,n]<=d_bound and tub[o,n]>0),GRB.MINIMIZE)
    m.update()
    run_time3=time.time()
    m.optimize()
    print("Budget")
    run_time3=time.time()-run_time3
    printSolution("minimizing budget",run_time3)
