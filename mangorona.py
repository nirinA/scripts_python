'''game of mangorona.

goal:
    keep more pawn on the board than your opponent.

movement:
    move your pawn to an unoccupied place.

pick:
    fill or create an empty place beetwen your pawn and
    your opponent's, and pick all opponent pawn in the
    same line of movement.
'''
import sys
import random
import time
import profile
import traceback

class IllegalMove(Exception):
    pass

class NoMoreMove(Exception):
    pass

class Init(object):
    def __init__(self, dimension, players, lattice):
        self.dimension = dimension
        self.xmax, self.ymax = dimension
        self.player1, self.player2, self.blank = players
        self.lattice = lattice
        self.all = [(x,y) for x in range(self.xmax) for y in range(self.ymax)]
        self.gain = {self.player1:0, self.player2:0}
        
class Position(Init):
    '''get all positions around one point'''
    def __init__(self, p, dimension, lattice):
        Init.__init__(self, dimension, ('','',''), lattice)
        self.xi, self.yi = p
        ##'''pawn can move only horizontally '''
        self.p1 = self.xi+1, self.yi
        self.p2 = self.xi-1, self.yi
        ##'''pawn can move only verticaly'''
        self.p3 = self.xi, self.yi+1
        self.p4 = self.xi, self.yi-1
        ##'''pawn can also move diagonaly'''
        self.p5 = self.xi-1, self.yi-1
        self.p6 = self.xi-1, self.yi+1
        self.p7 = self.xi+1, self.yi-1
        self.p8 = self.xi+1, self.yi+1
        
        if lattice is None:
            if sum(p)%2:
                self.around = self.p1,self.p2,self.p3,self.p4
            else:
                self.around = self.p1,self.p2,self.p3,self.p4,\
                              self.p5,self.p6,self.p7,self.p8
        elif lattice == 'star':
            if sum(p)%2:
                self.around = self.p1,self.p2,self.p3,self.p4
            else:
                self.around = self.p1,self.p2,self.p3,self.p4,\
                              self.p5,self.p6,self.p7,self.p8
        elif lattice == 'diamond':
            if sum(p)%2:
                self.around = self.p1,self.p2,self.p3,self.p4,\
                              self.p5,self.p6,self.p7,self.p8
            else:
                self.around = self.p1,self.p2,self.p3,self.p4
        elif lattice == 'cubic':
            self.around = self.p1,self.p2,self.p3,self.p4
        elif lattice == 'web':
            self.around = self.p1,self.p2,self.p3,self.p4,\
                          self.p5,self.p6,self.p7,self.p8
        elif lattice == 'X':
            self.around = self.p5,self.p6,self.p7,self.p8

    def Movable(self):
        return [p for p in self.around if p in self.all]

    def Deletable(self, final):
        xf, yf = final
        deltax = xf - self.xi
        deltay = yf - self.yi

        removeup = []
        removedown = []
        xu = xd = self.xi
        yu = yd = self.yi
        while (0<=xu<=self.xmax) and (0<=yu<=self.ymax):
            xu += deltax
            yu += deltay
            removeup.append((xu,yu))
        removeup.remove((xf, yf))
        while (0<=xd<=self.xmax) and (0<=yd<=self.ymax):
            xd -= deltax
            yd -= deltay
            removedown.append((xd,yd))
        return [xy for xy in removeup if xy in self.all],\
               [xy for xy in removedown if xy in self.all]

class Mangorona(Init):
    def __init__(self, players, lattice, dimension, matrix=None):
        '''set matrix to None to create an initial board with self.Create'''
        if matrix is None:
            self.matrix = self.Create(dimension, players)
        else:
            self.matrix = matrix
        Init.__init__(self, (len(self.matrix), len(self.matrix[0])), players, lattice)

    def Create(self, dimension, players):
        xmax, ymax = dimension
        player1, player2, blank = players
        m =[[None for i in range(ymax)] for j in range(xmax)]
        for x in range(xmax):
            for y in range(ymax):
                if (x < int(xmax/2)):
                    m[x][y]=player1
                elif (x == int(xmax/2)):
                    if (y < int(ymax/2)):
                        if y%2 != 0:
                            m[x][y]=player2
                        else:
                            m[x][y]=player1
                    elif (y == int(ymax/2)):
                        m[x][y]=blank
                    else:
                        if y%2 != 0:
                            m[x][y]=player1
                        else:
                            m[x][y]=player2
                else:
                    m[x][y]=player2
        return m

    def Zero(self):
        '''return the position(s) of blank'''
        w = []
        for i in range(self.xmax):
            c = self.matrix[i].count(self.blank)
            s = 0
            while c > 0:
                n = self.matrix[i].index(self.blank, s)
                w.append((i, n))
                s = n + 1
                c -= 1
        return w

    def Pawn(self, position, turn):
        x, y = position
        if self.matrix[x][y] == turn:
            return True
        return False

    def MovablePawn(self, turn):
        movable = []
        wherezero = self.Zero()
        for p in wherezero:
            pos = Position(p, self.dimension, self.lattice)
            turnmovable = [i for i in pos.Movable() if self.Pawn(i,turn)]
            movable.extend(turnmovable)
        return movable

    def ChangePawn(self, turn, initial, final):
        xi, yi = initial
        xf, yf = final
        self.matrix[xi][yi]=self.blank
        self.matrix[xf][yf]=turn
        todelete = Position(initial, self.dimension, self.lattice).Deletable(final)
        for t in todelete:
            for p in t:
                x,y = p
                if (not self.Pawn(p, turn) and self.matrix[x][y] != self.blank):
                    self.matrix[x][y] = self.blank
                    self.gain[turn] += 1
                else:
                    break

    def Move(self, turn, initial, final):
        if initial == final:
            raise IllegalMove("you don't move")
        if not self.Pawn(initial, turn):
            raise IllegalMove('not your pawn')
        if final not in self.Zero():
            raise IllegalMove('destination must be empty')
        if initial not in self.MovablePawn(turn):
            raise IllegalMove('this pawn cannot move')
        if final not in Position(initial, self.dimension, self.lattice).around:
            raise IllegalMove('not allowable move')
        self.ChangePawn(turn, initial, final)

    def Winner(self):
        if self.gain[self.player1]<self.gain[self.player2]:
            return self.player2
        elif self.gain[self.player1]>self.gain[self.player2]:
            return self.player1
        else:
            return self.blank

class AllowableMovement(object):
    def __init__(self, m, turn):
        self.m = m.matrix
        self.blank = m.blank
        self.mZero = m.Zero()
        self.mMovablePawn = m.MovablePawn(turn)
        self.mdimension = m.dimension
        self.player = turn
        self.mlattice = m.lattice

    def Move(self, maximum=False, getall=False):
        '''check if the player can move, and used as machine player'''
        move = {}
        for i in self.mMovablePawn:
            pos = Position(i, self.mdimension, self.mlattice)
            listf = [f for f in pos.around if f in self.mZero]
            for f in listf:
                if getall:
                    move.update({(i,f):0})
                else:
                    moveup , movedown = pos.Deletable(f)
                    up = [self.m[x][y] for (x,y) in moveup]
                    down = [self.m[x][y] for (x,y) in movedown]
                    if self.blank in up:
                        up = up[:up.index(self.blank)]
                    if self.player in up:
                        up = up[:up.index(self.player)]
                    if self.blank in down:
                        down = down[:down.index(self.blank)]
                    if self.player in down:
                        down = down[:down.index(self.player)]
                    get = len(up+down)
                    if get>0:
                        move.update({(i,f):get})
        if move:
            if maximum:
                getmax = max(move.values())
                for k in list(move.keys()):
                    if move[k]<getmax:
                        move.pop(k)
            return list(move.keys())
        else:
            raise NoMoreMove('%s cannot move anymore'%self.player)

class Board(object):
    '''displaying the game in command line mode'''
    def __init__(self, m):
        self.m = m.matrix
        self.x = m.xmax
        self.y = m.ymax
        self.evenline = [chr(92), '/']
        self.oddline = ['/', chr(92)]
        if m.lattice == 'diamond':
            self.evenline.reverse()
            self.oddline.reverse()
        if m.lattice == 'cubic':
            self.evenline = [' ', ' ']
            self.oddline = [' ', ' ']
        if m.lattice == 'web':
            self.evenline = ['x', 'x']
            self.oddline = ['x', 'x']
        
    def WidthLine(self, listline):
        if self.y%2==0:
            return '  |%s|'%'|'.join(listline*int(self.y/2))[:-2]
        return '  |%s|'%'|'.join(listline*int(self.y/2))
    def Inline(self, i):
        if i%2==0:
            return self.WidthLine(self.evenline)
        if i%2!=0:
            return self.WidthLine(self.oddline)

    def Display(self):
        d = '  '+' '.join([str(j) for j in range(self.y)])+'\n'
        for i in range(self.x):
            d += str(i)+' '
            d += '-'.join([str(self.m[i][j]) for j in range(self.y)])
            d += ' '+str(i)+'\n'
            if i!=self.x-1:
                d += self.Inline(i)+'\n'
        return d+'  '+' '.join([str(j) for j in range(self.y)])+'\n'
        
def MachineMachine():
    LATTICE = 'star' ##, 'diamond'
    DIMENSION = (5,11)
    PLAYERS = 'a', 'b', ' '
    ##mc = Mangorona(PLAYERS,'cubic', DIMENSION, None)
    ##maximum=True
    ##getall=True
    mc = Mangorona(PLAYERS,'diamond', (7,11), None)
    maximum=True
    getall=False
    t = PLAYERS[:2]
    tab = 0
    print(Board(mc).Display())
    while True:
        try:
            turn = t[tab%2]
            movable = AllowableMovement(mc, turn).Move(maximum=maximum, getall=getall)
            machine = random.choice(movable)
            print(turn, 'move:', machine[0], machine[1])
            mc.Move(turn, machine[0], machine[1])
            print(Board(mc).Display())
            print(mc.gain['a'], mc.gain['b']) ##, t1-t0
            print()
            tab += 1
        except IllegalMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            print(exc)
        except NoMoreMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            print(exc)
            print('winner:', mc.Winner())
            break

def TestvsMachine():
    LATTICE = 'star'
    DIMENSION = 5, 9
    PLAYERS = 'a', 'b', ' '
    machineplayer = PLAYERS[0]
    mc = Mangorona(PLAYERS,LATTICE, DIMENSION, None)
    maximum=True
    getall=False
    t = PLAYERS[:2]
    tab = 0
    print(Board(mc).Display())
    while True:
        try:
            turn = t[tab%2]
            movable = AllowableMovement(mc, turn).Move(maximum=maximum, getall=getall)
            if turn == machineplayer:
                machine = random.choice(movable)
                print(turn, 'move:', machine[0], machine[1])
                mc.Move(turn, machine[0], machine[1])
                print(Board(mc).Display())
                print(mc.gain['a'], mc.gain['b']) ##, t1-t0
                print()
                tab += 1
            else:
                h = input("type:'?' for movable, 'z' for Zero, 'h' for rules\nyour move - :")
                if h == '?':
                    print(mc.MovablePawn(turn))
                elif h == 'z':
                    print(mc.Zero())
                elif h == 'h':
                    print(__doc__)
                else:
                    human = eval(h)
                    if human not in movable:
                        raise IllegalMove('not allowable move')
                    mc.Move(turn, human[0], human[1])
                    print(Board(mc).Display())
                    tab += 1
        except IllegalMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            print(exc)
        except NoMoreMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            print(exc)
            print('winner:', mc.Winner())
            break
        except KeyboardInterrupt:
            raise SystemExit
        except:
            traceback.print_exc()

__version__ = '3k-0.0.0'
__author__ = 'nirinA'
__date__ = 'Sat May 10 21:52:15 2008'
