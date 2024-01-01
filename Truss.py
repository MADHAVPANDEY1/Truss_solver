import numpy as np                          ##import the python library "Numpy" for use
r=int(input("Enter the number of nodes:"))  ##enter the number of nodes in the truss
single_row = []                             
for i in range (0,r*2):                     ##take input of the x and y coordinates of the n nodes
    k=int(i/2)+1                            ##node number whose coordinates are being taken as input  
    #K=k+str(1)
    if(i%2==0):
        a=float(input("Enter the x co-ordinate of node " +str(k)+ ": "))
    else:
        a=float(input("Enter the y co-ordinate of node " +str(k)+ ": "))
    single_row.append(a)                    ##add to the array containing the coordinates

matrix_node=np.array(single_row).reshape(r,2)   ##reshape the array as r*2 2-D matrix
print(matrix_node)                              ##display this matrix

single_row1 = []                           ##take input of the x and y components of the loads  
for i in range (r*2):
    k=int(i/2)+1                           ##node number whose loads are being taken as input     
    if(i%2==0):
        a=float(input("Enter the load in x in node"+str(k)+": "))
    else:
        a=float(input("Enter the load in y in node"+str(k)+": "))
    single_row1.append(a)                  ##add to the array

matrix_node_loads=np.array(single_row1).reshape(r,2)   ##reshape the array as r*2 2-D matrix
print(matrix_node_loads)                               ##display this matrix 

single_row2 = []
s_r=[]                          ##array to store the node number and constraint direction of a node having non-zero constraint value as pairs
for i in range (1,r*2+1):       ##take input of the constraints at each node in the x and y directions
    if(i%2==0):
        a=float(input("Enter the constraint in y : "))
    else:
        a=float(input("Enter the constraint in x : "))
    single_row2.append(a)
    if(a==1):
        if(i%2==0):
            s_r.append((int((i+1)/2),2))
        else:
            s_r.append((int((i+1)/2),1))
        

matrix_constraints=np.array(single_row2).reshape(r,2)             ##reshape the array as r*2 2-D matrix  
print(matrix_constraints)                                         ##display this matrix     
print(s_r)                                                        ##display this array 

single_row3 = []
connected=[]                                                      ##array to store the pair of nodes that are connected by a beam in between
cnt=0                                                             ##used to count the number of beams
for i in range (r):
    for j in range(r):    
        a=int(input("Connection between "+str(i+1)+ " and "+ str(j+1)+ ": ")) 
        single_row3.append(a)
        if(a==1):
            cnt+=1
            if(i<j):
                connected.append((i+1,j+1))                      ##if the two nodes are connected append it to the array "connected" 

matrix_adjacency=np.array(single_row3).reshape(r,r)              ##reshape the array into 2-D adjacency matrix
print(matrix_adjacency)                                          ##display the adjacency matrix
beams=int(cnt/2)                                                 ##gives the number of beams in the truss
reactions=len(s_r)                                               ##give the value of total number of support reactions to be calculated
def x_position(a,b):                                             ##function which checks whether the node a is to the left or right of node b
    if(matrix_node[a-1][0]<matrix_node[b-1][0]):                 ##return 1 if a is left of b 
        return 1
    else:
        return 2                                                 ##return 2 if a is right of b
def y_position(a,b):                                             ##function which checks whether the node a is below or above the node b
    if(matrix_node[a-1][1]<matrix_node[b-1][1]):                 ##return 1 if a is below b 
        return 1
    else:                                                        ##return 2 if a is above of b 
        return 2
def Ha(a,b):                                                     ##function which finds the horizontal angle between nodes a and b
    p=abs(matrix_node[b-1][1]-matrix_node[a-1][1])
    b=abs(matrix_node[b-1][0]-matrix_node[a-1][0])
    h=(p**2 + b**2)**(1/2)
    return b/h
def Va(a,b):                                                    ##function which finds the vertical angle between nodes a and b
    p=abs(matrix_node[b-1][1]-matrix_node[a-1][1])
    b=abs(matrix_node[b-1][0]-matrix_node[a-1][0])
    h=(p**2 + b**2)**(1/2)
    return p/h



main_arr=[]
for j in range (1,r+1):     ##traverse through each node and find the coefficients of all the unknown member forces and support reactions in the x and y eqlm equations
    Fx=[]                               
    Fy=[]
    for i in range(beams):         ##evaluates the coefficients of force in each beam in x and y directions in a particular node being considered    
        if(connected[i][0]==j):
            Fx.append(pow(-1,(x_position(connected[i][1],j))) * Ha(connected[i][1],j))    
            Fy.append(pow(-1,(y_position(connected[i][1],j))) * Va(connected[i][1],j))   
            
        elif(connected[i][1]==j):
            Fx.append(pow(-1,(x_position(connected[i][0],j))) * Ha(connected[i][0],j))
            Fy.append(pow(-1,(y_position(connected[i][0],j))) * Va(connected[i][0],j))
            
        else:
            Fx.append(0)
            Fy.append(0)
    for a in range(reactions):    ##evaluates the coefficients of force at each support in x and y directions in a particular node being considered 
        if(s_r[a][0]==j):
            if(s_r[a][1]==1):
                Fx.append(1)
                Fy.append(0)
            else:
                Fx.append(0)
                Fy.append(1)
        else:
            Fx.append(0)
            Fy.append(0)

    main_arr.append(Fx)
    main_arr.append(Fy)
np.resize(main_arr,(2*r,1))   ##resize the array into 2*r x 2*r matrix
 


B=[]
for i in range(r):           ##add the node load at each node which are the constant values
    B.append(-1*matrix_node_loads[i][0])
    B.append(-1*matrix_node_loads[i][1])
np.resize(B,(2*r,1))         ##resize the array into 2*r x 1 matrix


X=np.linalg.solve(main_arr,B)    ##Solves the system AX=B using the Cramer's rule

for j in range(beams):           ##Output the member forces
    print("Force in member between " ,connected[j][0] ,"and" ,connected[j][1] ,"is", abs(X[j]))
for i in range(beams,2*r):       ##Output the support reactions
    print("reaction at support" ,s_r[i-beams][0])
    if (s_r[i-beams][1]==2):
        print("in the y direction is", abs(X[i]))
    else:
        print("in the x direction is", abs(X[i]))





