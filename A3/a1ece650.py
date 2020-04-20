#!/usr/bin/env python
import re
import sys

class point():
    def __init__(self, data): 
        self.x = data[0]
        self.y = data[1]

class linesegment():   
    def __init__(self,point1,point2):
        self.x1 = float(point1[0])
        self.y1 = float(point1[1])
        self.x2 = float(point2[0])
        self.y2 = float(point2[1])


def findpoint(point,LS_list):
    Key = []
    for i in range(len(LS_list)):
        for j in range(len(LS_list[i])):
            x1 = LS_list[i][j].x1
            y1 = LS_list[i][j].y1
            x2 = LS_list[i][j].x2
            y2 = LS_list[i][j].y2
            if point == (x1,y1):
                Key.append(i)
            elif point == (x2,y2):
                Key.append(i)
    
    ################################ 
    for i in range(len(LS_list)):
        for j in range(len(LS_list[i])):
            for k in range(i+1,len(LS_list)):
                for n in range(len(LS_list[k])):
                    x2 = max(LS_list[i][j].x1,LS_list[i][j].x2)
                    x1 = min(LS_list[i][j].x1,LS_list[i][j].x2)
                    x4 = max(LS_list[k][n].x1,LS_list[k][n].x2)
                    x3 = min(LS_list[k][n].x1,LS_list[k][n].x2)
                    y2 = max(LS_list[i][j].y1,LS_list[i][j].y2)
                    y1 = min(LS_list[i][j].y1,LS_list[i][j].y2)
                    y4 = max(LS_list[k][n].y1,LS_list[k][n].y2)
                    y3 = min(LS_list[k][n].y1,LS_list[k][n].y2)
                    if (x2 < x3) or (x1 > x4) or (y2 < y3) or (y1 > y4):
                        # no intersection
                        tmp = []
                    else:
                        # use cramers rule (determinants)
                        a1 = LS_list[i][j].y1 - LS_list[i][j].y2
                        b1 = LS_list[i][j].x2 - LS_list[i][j].x1
                        c1 = LS_list[i][j].x2*LS_list[i][j].y1 - LS_list[i][j].x1*LS_list[i][j].y2
                        a2 = LS_list[k][n].y1 - LS_list[k][n].y2
                        b2 = LS_list[k][n].x2 - LS_list[k][n].x1
                        c2 = LS_list[k][n].x2*LS_list[k][n].y1 - LS_list[k][n].x1*LS_list[k][n].y2
                        
                        
                        # calculate determinants
                        Det = a1*b2 - b1*a2
                        Detx = c1*b2 - b1*c2
                        Dety = a1*c2 - c1*a2
                        
                        # check determinant and for overlapping line segments
                        Testp1 = (LS_list[k][n].x1-LS_list[i][j].x1)*(LS_list[i][j].y2-LS_list[k][n].y2) - (LS_list[i][j].x2-LS_list[k][n].x2)*(LS_list[k][n].y1-LS_list[i][j].y1)
                        Testp2 = (LS_list[k][n].x2-LS_list[i][j].x1)*(LS_list[k][n].y1-LS_list[i][j].y2) - (LS_list[k][n].x1-LS_list[i][j].x2)*(LS_list[k][n].y2-LS_list[i][j].y1) 
                        flagp2 = False
                        if (Testp1 == 0 and Testp2 == 0) and (Det == 0):
                            #there is parallel... calculate it again
                            flagp = True
                            # use cramers rule (determinants)
    
                            Testp3 = (LS_list[i][j].x2-LS_list[i][j].x1)*(LS_list[k][n].y1-LS_list[k][n].y2) - (LS_list[k][n].x1-LS_list[k][n].x2)*(LS_list[i][j].y2-LS_list[i][j].y1)
                            
                            Testp4 = (LS_list[k][n].x2-LS_list[i][j].x1)*(LS_list[i][j].y2-LS_list[k][n].y1) - (LS_list[i][j].x2-LS_list[k][n].x1)*(LS_list[k][n].y2-LS_list[i][j].y1) 
                            
                            
                            if Testp3 == 0 and Testp4 == 0:
                                flagp = False
                                IP4 = [(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)]
                                IP4 = sorted(IP4)
                                X = []
                                Y = []
                                X.append(IP4[1][0])
                                X.append(IP4[2][0])
                                Y.append(IP4[1][1])
                                Y.append(IP4[2][1])
                                flagp2 = True
                        else:
                            flagp = False 
                        if Det != 0 and not flagp:
                            if Det!= 0:
                                # calcuate coords
                                X = Detx / Det
                                Y = Dety / Det
                                # check if its in the interval
                                if X < max(x1, x3) or X > min(x2,x4) or Y < max(y1,y3) or Y > min(y2,y4):
                                    #print("No intersection.")
                                    tmp = []
                                else:
                                    if (X,Y) == point:
                                        Key.append(i)
                                        Key.append(k)
                        elif(not flagp and flagp2):
                            for q in range(len(X)):
                                if X[q] < max(x1, x3) or X[q] > min(x2,x4) or Y[q] < max(y1,y3) or Y[q] > min(y2,y4):
                                    #print("No intersection.")
                                    tmp = []
                                else:
                                    if (X[q],Y[q]) == point:
                                        Key.append(i)
                                        Key.append(k)
    Key = list(dict.fromkeys(Key))                
    return(Key)

def graph(StreetDatabase):
    Points = StreetDatabase.points
    LS_list = []
    # compute line segments
    for i in Points:
        tmp = []
        LS = []
        for j in i:
            if not tmp:
                tmp = j
            else:
                LS.append(linesegment(j,tmp))
                tmp = j
        LS_list.append(LS)
    I_list = []
    V = []
    Key = []
    I_Points = []
    # compute any intersections (4 forloops deep to compare every segment)
    for i in range(len(LS_list)):
        for j in range(len(LS_list[i])):
            for k in range(i+1,len(LS_list)):
                for n in range(len(LS_list[k])):
                    x2 = max(LS_list[i][j].x1,LS_list[i][j].x2)
                    x1 = min(LS_list[i][j].x1,LS_list[i][j].x2)
                    x4 = max(LS_list[k][n].x1,LS_list[k][n].x2)
                    x3 = min(LS_list[k][n].x1,LS_list[k][n].x2)
                    y2 = max(LS_list[i][j].y1,LS_list[i][j].y2)
                    y1 = min(LS_list[i][j].y1,LS_list[i][j].y2)
                    y4 = max(LS_list[k][n].y1,LS_list[k][n].y2)
                    y3 = min(LS_list[k][n].y1,LS_list[k][n].y2)
                    if (x2 < x3) or (x1 > x4) or (y2 < y3) or (y1 > y4):
                        #print('no intersection')
                        tmp = []
                    else:
                        #######################################################
                        # use cramers rule (determinants)
                        a1 = LS_list[i][j].y1 - LS_list[i][j].y2
                        b1 = LS_list[i][j].x2 - LS_list[i][j].x1
                        c1 = LS_list[i][j].x2*LS_list[i][j].y1 - LS_list[i][j].x1*LS_list[i][j].y2
                        a2 = LS_list[k][n].y1 - LS_list[k][n].y2
                        b2 = LS_list[k][n].x2 - LS_list[k][n].x1
                        c2 = LS_list[k][n].x2*LS_list[k][n].y1 - LS_list[k][n].x1*LS_list[k][n].y2
                                              
                        # calculate determinants
                        Det = a1*b2 - b1*a2
                        Detx = c1*b2 - b1*c2
                        Dety = a1*c2 - c1*a2
                            
                        # check determinant and for overlapping line segments
                        Testp1 = (LS_list[k][n].x1-LS_list[i][j].x1)*(LS_list[i][j].y2-LS_list[k][n].y2) - (LS_list[i][j].x2-LS_list[k][n].x2)*(LS_list[k][n].y1-LS_list[i][j].y1)
                        Testp2 = (LS_list[k][n].x2-LS_list[i][j].x1)*(LS_list[k][n].y1-LS_list[i][j].y2) - (LS_list[k][n].x1-LS_list[i][j].x2)*(LS_list[k][n].y2-LS_list[i][j].y1) 
                        flagp2 = False
                        if (Testp1 == 0 and Testp2 == 0) and (Det == 0):
                            #there is parallel... calculate it again
                            flagp = True
                            # use cramers rule (determinants)

                            Testp3 = (LS_list[i][j].x2-LS_list[i][j].x1)*(LS_list[k][n].y1-LS_list[k][n].y2) - (LS_list[k][n].x1-LS_list[k][n].x2)*(LS_list[i][j].y2-LS_list[i][j].y1)
                            Testp4 = (LS_list[k][n].x2-LS_list[i][j].x1)*(LS_list[i][j].y2-LS_list[k][n].y1) - (LS_list[i][j].x2-LS_list[k][n].x1)*(LS_list[k][n].y2-LS_list[i][j].y1) 
                            
                            if Testp3 == 0 and Testp4 == 0:
                                flagp = False
                                flagp2 = True
                                IP4 = [(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)]
                                IP4 = sorted(IP4)
                                X = []
                                Y = []
                                X.append(IP4[1][0])
                                X.append(IP4[2][0])
                                Y.append(IP4[1][1])
                                Y.append(IP4[2][1])
                        else:
                            flagp = False 
                        if Det != 0 and not flagp:
                            if Det!= 0:
                                # calcuate coords
                                X = Detx / Det
                                Y = Dety / Det
                            if X < max(x1, x3) or X > min(x2,x4) or Y < max(y1,y3) or Y > min(y2,y4):
                                #print("No intersection.")
                                tmp = []
                            else:
                                #print("Intersection at coordinates: " + "("+str(X)+","+str(Y)+")")
                                
                                # flags for matching the line points for verts.
                                flagi = False
                                flagv = False
                                flag1 = False
                                flag2 = False
                                flag3 = False
                                flag4 = False
                                
                                # check vertices for line segment points existence
                                for c in range(len(V)):
                                    if (LS_list[i][j].x1,LS_list[i][j].y1) == V[c]:
                                        flag1 = True
                                # add line segment point 1 to vertices if they aren't already
                                if not flag1:
                                    V.append((LS_list[i][j].x1,LS_list[i][j].y1))
                                 
                                for c in range(len(V)):
                                    if (LS_list[i][j].x2,LS_list[i][j].y2) == V[c]:
                                        flag2 = True     
                                # add line segment point 2 to vertices if they aren't already
                                if not flag2:
                                    V.append((LS_list[i][j].x2,LS_list[i][j].y2))
                                                                       
                                for c in range(len(V)):
                                    if (LS_list[k][n].x1,LS_list[k][n].y1) == V[c]:
                                        flag3 = True
                                # add line segment point 3 to vertices if they aren't already
                                if not flag3:
                                    V.append((LS_list[k][n].x1,LS_list[k][n].y1)) 
                                
                                for c in range(len(V)):
                                    if (LS_list[k][n].x2,LS_list[k][n].y2) == V[c]:
                                        flag4 = True
                                # add line segment point 3 to vertices if they aren't already
                                if not flag4:
                                    V.append((LS_list[k][n].x2,LS_list[k][n].y2))
                                
                                # check vertices and intersection lists for intersection existence
                                for c in range(len(V)):
                                    if (V[c] == (X,Y)):
                                        flagv = True
                                idxi = []
                                for c in range(len(I_list)):
                                    if I_list[c] == (X,Y):
                                        flagi = True
                                        idxi = c
                                
                                # add intersection if does not exist in the corresponding list
                                if not flagi and not flagv:
                                    I_list.append((X,Y))
                                    V.append((X,Y))
                                    I_Points.append([(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)])
                                elif not flagi and flagv:
                                    I_list.append((X,Y))  
                                    I_Points.append([(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)])
                                else:
                                    I_Points[idxi].append((LS_list[i][j].x1,LS_list[i][j].y1))
                                    I_Points[idxi].append((LS_list[i][j].x2,LS_list[i][j].y2))
                                    I_Points[idxi].append((LS_list[k][n].x1,LS_list[k][n].y1))
                                    I_Points[idxi].append((LS_list[k][n].x2,LS_list[k][n].y2))  
                        elif(not flagp and flagp2):
                            for q in range(len(X)):
                                if X[q] < max(x1, x3) or X[q] > min(x2,x4) or Y[q] < max(y1,y3) or Y[q] > min(y2,y4):
                                    #print("No intersection.")
                                    tmp = []
                                else:
                                    #print("Intersection at coordinates: " + "("+str(X)+","+str(Y)+")")
                                    
                                    # flags for matching the line points for verts.
                                    flagi = False
                                    flagv = False
                                    flag1 = False
                                    flag2 = False
                                    flag3 = False
                                    flag4 = False
                                    
                                    # check vertices for line segment points existence
                                    for c in range(len(V)):
                                        if (LS_list[i][j].x1,LS_list[i][j].y1) == V[c]:
                                            flag1 = True
                                    # add line segment point 1 to vertices if they aren't already
                                    if not flag1:
                                        V.append((LS_list[i][j].x1,LS_list[i][j].y1))
                                     
                                    for c in range(len(V)):
                                        if (LS_list[i][j].x2,LS_list[i][j].y2) == V[c]:
                                            flag2 = True     
                                    # add line segment point 2 to vertices if they aren't already
                                    if not flag2:
                                        V.append((LS_list[i][j].x2,LS_list[i][j].y2))
                                                                           
                                    for c in range(len(V)):
                                        if (LS_list[k][n].x1,LS_list[k][n].y1) == V[c]:
                                            flag3 = True
                                    # add line segment point 3 to vertices if they aren't already
                                    if not flag3:
                                        V.append((LS_list[k][n].x1,LS_list[k][n].y1)) 
                                    
                                    for c in range(len(V)):
                                        if (LS_list[k][n].x2,LS_list[k][n].y2) == V[c]:
                                            flag4 = True
                                    # add line segment point 3 to vertices if they aren't already
                                    if not flag4:
                                        V.append((LS_list[k][n].x2,LS_list[k][n].y2))
                                    
                                    # check vertices and intersection lists for intersection existence
                                    for c in range(len(V)):
                                        if (V[c] == (X[q],Y[q])):
                                            flagv = True
                                    idxi = []
                                    for c in range(len(I_list)):
                                        if I_list[c] == (X[q],Y[q]):
                                            flagi = True
                                            idxi = c
                                    
                                    # add intersection if does not exist in the corresponding list
                                    if not flagi and not flagv:
                                        I_list.append((X[q],Y[q]))
                                        V.append((X[q],Y[q]))
                                        I_Points.append([(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)])
                                    elif not flagi and flagv:
                                        I_list.append((X[q],Y[q]))  
                                        I_Points.append([(LS_list[i][j].x1,LS_list[i][j].y1),(LS_list[i][j].x2,LS_list[i][j].y2),(LS_list[k][n].x1,LS_list[k][n].y1),(LS_list[k][n].x2,LS_list[k][n].y2)])
                                    else:
                                        I_Points[idxi].append((LS_list[i][j].x1,LS_list[i][j].y1))
                                        I_Points[idxi].append((LS_list[i][j].x2,LS_list[i][j].y2))
                                        I_Points[idxi].append((LS_list[k][n].x1,LS_list[k][n].y1))
                                        I_Points[idxi].append((LS_list[k][n].x2,LS_list[k][n].y2)) 
    for i in range(len(Points)):
        Points[i] = list(dict.fromkeys(Points[i]))
    # Compute edges from list of intersections and vertices
    V = sorted(V)
    for i in range(len(V)):
        Key.append(findpoint(V[i],LS_list))
    
    V2 = [[]]*len(LS_list)
    for i in range(len(V)):
        for j in range(len(LS_list)):
            if Key[i].count(j)>0:
                if not V2[j]:
                    V2[j] = [V[i]]
                V2[j].append(V[i])
    # remove dupes
    for i in range(len(V2)):
        V2[i] =  list(dict.fromkeys(V2[i]))     
    for i in range(len(I_Points)):
        I_Points[i] = list(dict.fromkeys(I_Points[i]))
 
    # now we "sort" vsorted again by order of point connection. this is to facilitate edge calculation later..
    Vsorted = [[]]*len(V2)
    for i in range(len(V2)):
        if(len(V2[i])>len(Points[i])):         
            Vsorted[i] = [None]*len(V2[i])
        else:
            Vsorted[i] = [None]*len(Points[i])
        for j in range(len(V2[i])):
            flag = False
            for n in range(len(I_list)):
                if V2[i][j] == I_list[n]:
                    flag = True
            idx = []
            for k in range(len(Points[i])):
                if (float(Points[i][k][0]),float(Points[i][k][1])) == V2[i][j]:
                    idx = k
                    flag = False
            if not flag:
                Vsorted[i][idx] = V2[i][j]
        Vsorted[i] = list(filter(None,Vsorted[i]))
    
    # repeat above for the intersections now
    for i in range(len(V2)):
        for j in range(len(V2[i])):
            flag = False
            flag2 = False
            idx = []
            for n in range(len(I_list)):
                for m in range(len(Vsorted[i])):
                    if V2[i][j] == Vsorted[i][m]:
                        flag2 = True
                if V2[i][j] == I_list[n]:
                    flag = True
                    IP = I_Points[n]*1
                    idx = n
            if flag and not flag2:
                idx1 = []
                IP_V = []
                # match intersection points with the verts already in Vsorted
                for k in range(len(IP)):
                    for w in range(len(Vsorted[i])):
                        if (IP[k] == Vsorted[i][w]):
                            idx1.append(w)
                            IP_V.append(IP[k])
                if(max(idx1) == (min(idx1)+1)):
                    Vsorted[i].insert(max(idx1),V2[i][j])
                else:
                    for h in range(1,len(IP_V)):
                        x = IP_V[h-1][0]-IP_V[h][0]
                        y = IP_V[h-1][1]-IP_V[h][1]
                        for v in range(min(idx1), max(idx1)):
                            # sort by ascending x
                            if(x > 0):
                                if((Vsorted[i][v][0] < V2[i][j][0]) and (Vsorted[i][v+1][0] > V2[i][j][0])):
                                    Vsorted[i].insert(v+1,V2[i][j])
                            # sort by descending x
                            elif(x < 0):
                                if((Vsorted[i][v][0] > V2[i][j][0]) and (Vsorted[i][v+1][0] < V2[i][j][0])):
                                    Vsorted[i].insert(v+1,V2[i][j])
                            # sort by ascending y
                            elif(y > 0):     
                                if((Vsorted[i][v][1] < V2[i][j][1]) and (Vsorted[i][v+1][1] > V2[i][j][1])):
                                    Vsorted[i].insert(v+1,V2[i][j])                  
                            # sort by descending y
                            elif(y < 0):
                                if((Vsorted[i][v][1] > V2[i][j][1]) and (Vsorted[i][v+1][1] < V2[i][j][1])):
                                    Vsorted[i].insert(v+1,V2[i][j])
                                
    # after all that we can finally calculate which segments are edges woo
    Edge = []
    for i in range(len(Vsorted)):
        for j in range(len(Vsorted[i])-1):
            flag1 = False
            flag2 = False
            for m in range(len(I_list)):
                if Vsorted[i][j] == I_list[m]:
                    flag1 = True
                elif Vsorted[i][j+1] == I_list[m]:
                    flag2 = True
            if flag1 or flag2:
                Edge.append(Vsorted[i][j])
                Edge.append(Vsorted[i][j+1])
     
    E = []       
    # get rid of duplicate edges and finish list            
    for i in range(0,len(Edge)-1,2):
        idx1 = []
        idx2 = []
        for m in range(len(V)):
            if V[m] == Edge[i]:
                idx1 = m
            if V[m] == Edge[i+1]:
                idx2 = m
        if not E:
            E.append((idx1,idx2))
        else:
            flag1 = False
            for j in range(0,len(E)):
                if (E[j][0] == idx1) and (E[j][1] == idx2):
                    flag1 = True
                elif (E[j][0] == idx2) and (E[j][1] == idx1):
                    flag1 = True
            if not flag1:
                E.append((idx1,idx2))                   
    E = list(dict.fromkeys(E))
                  
    return(V, E)

class streetdatabase():
    def __init__(self):
        self.names = []
        self.points = []
    def addstreet(self, street):
        # use reg expressions to determine if street is here
        flag = False
        if street[0]:
            j = 0
            for i in self.names:
                if(i.replace("'","").replace("\"","").lower().find(street[0].replace("'","").replace("\"","").lower()) != -1):
                    flag = True
                j = j+1
            if (flag == True):   
                print("Error: 'a' specified for a street that already exists! Keeping old entry")
            else:
                self.names.append(street[0])
                self.points.append(street[1])
        else:
            print('Error: You did not enter a street name or entered an incorrect name. Try again!')
    def removestreet(self, street):
        # use reg expressions to determine if street is here
        flag = False
        if street[0]:
            j = 0
            for i in self.names:
                if(i.replace("'","").replace("\"","").lower().find(street[0].replace("'","").replace("\"","").lower()) != -1):
                    flag = True
                    idx = j
                j = j+1
            if (flag == True):
                # remove data
                self.names.pop(idx)
                self.points.pop(idx)
            else:
                print("Error: 'r' specified for a street that does not exist.")
        else:
            print('Error: You did not enter a street name or entered an incorrect name. Try again!')
    def modifystreet(self,street):
        # use reg expressions to determine if street is here
        flag = False
        if street[0]:
            j = 0
            for i in self.names:
                if(i.replace("'","").replace("\"","").lower().find(street[0].replace("'","").replace("\"","").lower()) != -1):
                    flag = True
                    idx = j
                j = j + 1
            if (flag == True):
                # change street data
                self.points[idx] = street[1]
            else:
                print("Error: 'c' specified for a street that does not exist.")
        else:
            print('Error: You did not enter a street name. Try again!')
def main():
    Data = streetdatabase()
    try:
        while True:
            Command = input("")
            # Reg. expression to capture input command + data
            c = re.compile('^\s*\w\s+')
            g = re.compile('^\s*g\s*$')
            s = re.compile('["][a-zA-Z\s]+["]\s+')
            r = re.compile('["][a-zA-Z\s]+["]\s*')
            p = re.compile('(\s*[-+]?[\d.]+\s*),(\s*[-+]?[\d.]+\s*)')
            
            flagpo = True
            # determine whether the round brackets exist
            flagb = True
            if "(" not in Command or ")" not in Command:
                flagb = False
            brackets = 0
            for i in Command:
                if i == "(":
                    brackets = brackets + 1
                elif i == ")":
                    if brackets == 0:
                        flagb = False
                    brackets = brackets - 1
            flagb = brackets == 0
            if not flagb:
                print('Error: Brackets are missing!')
            # Check if command exists
            try:
                Com2 = g.findall(Command)
                Com = c.search(Command).group()
                # pull out just the command character
                Com = re.search('\w',Com).group()
                Flag = True
            except:
                if(Com2):
                    if(Command=='g'):
                        Flag = True
                        Com = 'g'
                    else:
                        Flag = False
                        Com = ''
                else:
                    print("Error: You did not enter any recognizable command or command sequence. Try again.")
                    Flag = False
                    Com = ''
            if ((Com in {'a','c','r','g'})) and (Flag == True) and flagb:
                streetname = s.findall(Command)
                if streetname:
                    streetname = re.findall('[\'"][a-zA-Z\s]+[\'"]',streetname[0])
                elif(Com in {'r'}):
                    streetname = r.findall(Command)
                elif(Com in {'g'}):
                    streetname = r.findall(Command)
                else:
                    print("Error: Make sure to leave proper whitespacing in your command!")
                points = p.findall(Command)
                points = list(dict.fromkeys(points))
                if len(points)>1 or Com in {'r','g'}:
                    if streetname or (not streetname and Com == "g"):
                        if Com == "a":
                            if points:
                                Data.addstreet([streetname[0],points])
                                #print(Data.names)
                                #print(Data.points)
                            else:
                                print("Error: 'a' specified with no points given.")
                        elif Com == "c":
                            if points:
                                Data.modifystreet([streetname[0],points])
                                #print(Data.names)
                                #print(Data.points)
                            else:
                                print("Error: 'c' specified with no points given.")
                        elif Com == "g":
                            # Depending on reg exp. detection, will change array size
                            if points and streetname:
                                print("Error: Ignoring street and point entry with command 'g'.")
                            elif points:
                                print("Error: Ignoring point entry with command 'g'.")
                            elif streetname:
                                print("Error: Ignoring street entry with command 'g'.")
                            elif not Com2:
                                print("Error: Invalid command entered with command 'g'.")
                            else:    
                                [V,E] = graph(Data)
                                edge = 'E {'
                                vert = "V = {"
                                for i in range(len(V)):
                                    vert = vert + '\n  ' + str(i+1) + ':' + '  (' + str(round(V[i][0],2)) + ',' + str(round(V[i][1],2)) + ')'
                                vert = vert + '\n}'
                                for i in range(len(E)-1):
                                    edge = edge + '<' + str(round(E[i][0],2)) + ',' + str(round(E[i][1],2)) + '>' + ','
                                if E:
                                    edge = edge + '<' + str(round(E[len(E)-1][0],2)) + ',' + str(round(E[len(E)-1][1],2)) + '>'
                                edge = edge + '}'
                                Vprint = "V ";
                                Vprint = Vprint + str(len(V))
                                print(Vprint)
                                print(edge)
                        else:
                            Data.removestreet([streetname[0],points])
                    else:
                        print("Error: You did not enter a streetname! Try again!")
                else:
                    print("Error: You've not entered enough points for your street.")
                    flagpo = False
            elif((Flag == True) and not(Com in {'a','c','r','g'}) and flagb and flagpo):
                print(Com+" is not a correct command. Try again.")
            else:
                print('Error: Command sequence not recognized. Try again!')
    except EOFError:
        return
            
#    ### sample code to read from stdin.
#    ### make sure to remove all spurious print statements as required
#    ### by the assignment
#    while True:
#        line = sys.stdin.readline()
#        if line == '':
#            break
#        print('read a line:', line)
#
#    print('Finished reading input')
#    # return exit code 0 on successful termination
#    sys.exit(0)

if __name__ == '__main__':
    main()

