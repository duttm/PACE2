import math

#BFS traverse to create list of all clusters, with the largest cluster being the one we consider    
def main(data,dimensions,adjacencyList):
    clusterList = []

    #Mark all the vertices as not visited 
    visited = [False] * (len(data)) 

    #BFS traverse
    for i in range(len(data)):
        if visited[i] == False:
            curClusterList = []
            curClusterList = BFS(data,i,visited,adjacencyList)
            clusterList.append(curClusterList)

    #Find largest cluster
    largestClusterIndex = -1
    maxClusterLength = 0

    for i in range(len(clusterList)):
        if len(clusterList[i]) > maxClusterLength:
            maxClusterLength = len(clusterList)
            largestClusterIndex = i

    #Return largest cluster
    return clusterList[largestClusterIndex]

def BFS(data,s,visited,adjacencyList): 

    #Create a queue for BFS 
    queue = [] 

    #Create cluster list
    cluster = []

    # Mark the source node as visited and enqueue it 
    queue.append(s) 
    visited[s] = True

    #Keep track of previous 
    prev = s

    while len(queue) > 0: 

        # Dequeue a vertex from  
        # queue and print it 
        s = queue.pop(0) 
        #print (s, end = " ") 
        cluster.append(s)

        # Get all adjacent vertices of the 
        # dequeued vertex s. If a adjacent 
        # has not been visited, then mark it 
        # visited and enqueue it 
        for i in (adjacencyList[s]): 
            if visited[i] == False: 
                queue.append(i) 
                visited[i] = True

    return cluster