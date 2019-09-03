from sys import argv
from time import time
from math import ceil, sqrt
from heapq import heappop, heappush
from collections import defaultdict
from pprint import pprint



def octileDist(a, b):
    dx = abs(a.x - b.x)
    dy = abs(a.y - b.y)
    return max(dx, dy) + ((sqrt(2)-1)*min(dx, dy))


class Cell:
    def __init__(self, x, y, c, block):
        self.x = x
        self.y = y
        self.empty = c in ['.', 'S']
        self.block = block
        self.neighboursin = []
        self.neighboursout = []
        self.g = float("inf")
    def coord(self):
        return (self.x, self.y)
    def reset(self):
        self.g = float('inf')
    def __str__(self):
        return '<Cell ({},{})>'.format(self.x, self.y)
    def __repr__(self):
        return str(self)
    def __lt__(self, other):
        return self.coord() < other.coord()


class Block:
    def __init__(self, x, y, grid):
        self.x = x
        self.y = y
        self.grid = grid
        self.ingress = set()
        self.egress = []
        self.heapval = float('inf')
    def coord(self):
        return (self.x, self.y)
    def reset(self):
        self.ingress = set()
        self.heapval = float('inf')
    def __str__(self):
        return '<Block ({},{})>'.format(self.x, self.y)
    def __repr__(self):
        return str(self)
    def __lt__(self, other):
        return self.coord() < other.coord()


class Grid:
    def __init__(self, filename, K):
        with open('maps/{}.map'.format(filename)) as file:
            lines = file.readlines()            
            self.direction = lines[0].split()[1]
            self.height = H = int(lines[1].split()[1])
            self.width = W = int(lines[2].split()[1])
            self.K = K
            self.block = [[Block(i, j, self) 
                            for j in range(ceil(W/K))] 
                            for i in range(ceil(H/K))]
            self.cell = [[Cell(i, j, lines[i+4][j], self.block[i//K][j//K])
                            for j in range(W)]
                            for i in range(H)]
        for i in range(H):
            for j in range(W):
                currCell = self.cell[i][j]
                if currCell.empty:
                    # Build cell's neighbours tile direction
                    for (di, dj) in [(0,1), (0,-1), (1,0), (-1,0)]:
                        ii = i+di
                        jj = j+dj
                        if (0<=ii<H) and (0<=jj<W):
                            nextCell = self.cell[ii][jj]
                            if nextCell.empty:
                                distance = octileDist(currCell, nextCell)
                                neighbour = (nextCell,distance)
                                if currCell.block == nextCell.block:
                                    currCell.neighboursin.append(neighbour)
                                else:
                                    currCell.neighboursout.append(neighbour)
                    # Build cell's neighbours octile direction
                    for (di, dj) in [(1,1), (1,-1), (-1,1), (-1,-1)]:
                        ii = i+di
                        jj = j+dj
                        if (0<=ii<H) and (0<=jj<W):
                            nextCell = self.cell[ii][jj]
                            sideCell1 = self.cell[i][jj]
                            sideCell2 = self.cell[ii][j]
                            if nextCell.empty and sideCell1.empty and sideCell2.empty:
                                distance = octileDist(currCell, nextCell)
                                neighbour = (nextCell,distance)
                                if currCell.block == nextCell.block:
                                    currCell.neighboursin.append(neighbour)
                                else:
                                    currCell.neighboursout.append(neighbour)
                    # Build block's egress cells
                    if ((i%K==0) or (j%K==0) or ((i+1)%K==0) or ((j+1)%K==0)) and currCell.neighboursout:
                        currCell.block.egress.append(currCell)
    def cell_at(self, coord):
        return self.cell[coord[0]][coord[1]]
    def reset(self):
        for blockrow in self.block:
            for block in blockrow:
                block.reset()
        for cellrow in self.cell:
            for cell in cellrow:
                cell.reset()


def dijkstra(startCell):
    dist = defaultdict(lambda: float('inf'))
    dist[startCell] = 0.0
    pq = [(dist[startCell], startCell)]
    while pq:
        pathLength, currCell = heappop(pq)
        if pathLength > dist[currCell]:
            continue
        for nextCell, weight in currCell.neighboursin:
            if dist[nextCell] > dist[currCell] + weight:
                dist[nextCell] = dist[currCell] + weight
                heappush(pq, (dist[nextCell], nextCell))
    return dist


class LDDB:
    def __init__(self, grid):
        self.grid = grid
        self.dist = dict()
        for cellrow in grid.cell:
            for cell in cellrow:
                self.dist[cell] = dijkstra(cell)


def initBlock(lddb, cell, tipe):
    block = cell.block
    if tipe == 'start':
        for c in block.egress:
            c.g = lddb.dist[cell][c]
            if (c.g != float('inf')):
                block.ingress.add(c)
    return block


def blockAStar(grid, lddb, start, goal):
    startBlock = initBlock(lddb, start, 'start')
    goalBlock = initBlock(lddb, goal, 'goal')
    if startBlock == goalBlock:
        length = lddb.dist[start][goal]
    else:
        length = float('inf')
    pq = [(octileDist(start,goal), startBlock)]
    while pq:
        currheapval, currBlock = heappop(pq)
        if currheapval >= length:
            break
        if currheapval > currBlock.heapval:
            continue
        if currBlock == goalBlock:
            for y in currBlock.ingress:
                length = min(length, y.g + lddb.dist[goal][y])
        newheapval = defaultdict(lambda: float('inf'))
        for xin in currBlock.egress:
            for y in currBlock.ingress:
                xin.g = min(xin.g, y.g + lddb.dist[xin][y])
            for xout,weight in xin.neighboursout:
                if xout.g > xin.g + weight:
                    xout.g = xin.g + weight
                    xout.block.ingress.add(xout)
                    if newheapval[xout.block] > xout.g + octileDist(xout, goal):
                        newheapval[xout.block] = xout.g + octileDist(xout, goal)
        for block,heapval in newheapval.items():
            block.heapval = heapval
            heappush(pq, (heapval, block))
        currBlock.ingress = set()
    return length


def test_single():
    grid = Grid(argv[1], int(argv[2]))
    start = grid.cell_at((int(argv[3]), int(argv[4])))
    goal = grid.cell_at((int(argv[5]), int(argv[6])))
    lddb = LDDB(grid)
    tstart = time()
    length = blockAStar(grid, lddb, start, goal)
    tend = time()
    duration = (tend - tstart)*1000
    print(length)
    print(duration, 's')


def test_all():
    filename = argv[1]
    lowK = argv[2]
    higK = argv[3]
    with open('maps/{}.map.scen'.format(filename)) as file:
        lines = file.readlines()            
        testcases = [row.split('\t') for row in lines[1:]]
    for K in range(int(lowK),int(higK)+1):
        grid = Grid(filename, K)
        lddb = LDDB(grid)
        failed = 0
        totalDuration = 0.0
        for tc in testcases:
            tstartreset = time()
            grid.reset()
            tendreset = time()
            start = grid.cell_at((int(tc[5]), int(tc[4])))
            goal = grid.cell_at((int(tc[7]), int(tc[6])))
            expected = float(tc[8])
            tstart = time()
            length = blockAStar(grid, lddb, start, goal)
            tend = time()
            duration = tend - tstart + tendreset - tstartreset
            totalDuration += duration
            # print('{:.3f}   K={}   x1=({},{})   x2=({},{})   {:.3f}'.format(abs(length-expected), K, tc[5], tc[4], tc[7], tc[6], length))
            if not (abs(length-expected)<1e-6):
                failed += 1
        # print()
        # print('Failed={}   K={}   duration={:.3f} s'.format(failed, K, totalDuration))
        print(totalDuration/len(testcases))
        # print()


if __name__ == '__main__':
    test_all()
    # test_single()