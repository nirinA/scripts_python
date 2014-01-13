'''bataille navale

. les bateaux doivent être espacé d'une case.
. taille maxi est de (26, 26)
. pour plus de bateaux definir differents noms comme:
    NAVIRES = {'porteavion1':5,
               'porteavion2':4
               'croiseur': 3,
               'corvette':2,
               'sousmarin1':1,
               'sousmarin2':2}
. si plus de navires ou plus grande taille de navire
  vérifier que la taille du jeu peut supporter ces nombres de navires.
. machine_start = False ou machine_start = 0 pour commencer,
  sinon, c'est la machine qui commence

'''
import sys
import random

machine_starts = True
TAILLE = (10,10) ##(12, 12)
NAVIRES = {'porteavion':4,
           'croiseur': 3,
           'corvette':2,
           'sousmarin':1}

class Error(Exception):
    pass

class GameOver(Exception):
    pass
    
class Setup(object):
    '''machine setup grid'''
    def __init__(self, size=(10,10)):
        try:
            self.size = self.check_size(size)
        except Warning as w:
            print(w)
            self.size = 10, 10
        self.mesh = []
        for i in range(self.size[0]):
            for j in range(self.size[1]):
                self.mesh.append((i,j))
        self.player_target = self.mesh[:]
        self.player_position = {i:list(range(NAVIRES[i])) for i in NAVIRES}
        self.machine_target = self.mesh[:]
        self.machine_ships = self.mesh[:]
        self.machine_position = {} 
        self.set_ships(NAVIRES)
        self.player_missed = []
        self.player_touched = []
        self.player_sank = []
        self.machine_missed = []
        self.machine_touched = []
        self.where = ''
        for i in self.machine_position:
            self.where += '%s :'%i
            for x,y in self.machine_position[i]:
                self.where += '%s%d, '%(chr(x+65), y)
            self.where += '\n'
        #print(self.machine_position)
        for i in self.machine_position.values():
            self.machine_touched.extend(i)
        self.machine_sank = []

    def check_size(self, size):
        if (not (8 < size[0] <= 26)) or (not (8 < size[1] <= 26)):
            raise Warning('trop petite ou rop grande taille, le définir à 10x10')
        else:
            return size

    def check_around_position(self, position, touched=False):
        x, y = position
        if touched:
            return [(x+1, y), (x, y-1),
                    (x, y+1), (x-1, y)]
        return [(x+1, y+1), (x+1, y),
                (x+1, y-1), (x, y-1),
                (x, y+1), (x-1, y-1),
                (x-1, y), (x-1, y+1)]
    
    def check_around(self, position, ship_length):
        x, y = position
        p1 = [(x, y+i) for i in range(1, ship_length)]
        p2 = [(x+i, y) for i in range(1, ship_length)]
        p3 = [(x-i, y) for i in range(1, ship_length)]
        p4 = [(x, y-i) for i in range(1, ship_length)]
        return p1, p2, p3, p4

    def valid_position(self, ship):
        ship_length = NAVIRES[ship]
        s = random.choice(self.machine_ships)
        if ship_length == 1:
            return [s]
        add_pos = []
        for j in range(1, ship_length):
            for p in self.check_around(s, ship_length):
                if set(p).issubset(set(self.machine_ships)):
                    add_pos.append(p)
        if add_pos:
            a = random.choice(add_pos)
            a.append(s)
            return a
        else:
            return None

    def remove_around(self, position):
        pos = []
        for p in position:
            pos.extend(self.check_around_position(p))
        return [i for i in list(set(pos)) if i in self.machine_ships]

    def set_ships(self, ships):
        for i in ships:
            ship_position = None
            while ship_position is None:
                ship_position = self.valid_position(i)
            around_position = self.remove_around(ship_position)
            self.machine_position[i] = ship_position
            for p in around_position:
                self.machine_ships.pop(self.machine_ships.index(p))                

    def parse_shoot(self, shoot):
        if len(shoot) == 0:
            raise Error('pas dans la grille, réessayez!')
        x = shoot[0]
        if not x.isalpha():
            raise Error('pas dans la grille, réessayez!')
        else:
            x = ord(x.upper()) - 65
            if not (0 <= x < self.size[0]):
                raise Error('pas dans la grille, réessayez!')
        y = shoot[1:]
        if not y.isnumeric():
            raise Error('pas dans la grille, réessayez!')
        else:
            y = int(y) - 1
            if not (0 <= y < self.size[1]):
                raise Error('pas dans la grille, réessayez!')
        return x, y
        
    def shoot_target(self, shoot):
        shoot = self.parse_shoot(shoot)
        if shoot in self.machine_touched:
            self.machine_sank.append(shoot)
            if len(self.machine_sank) == len(self.machine_touched):
                raise GameOver('argh! toutes mes navires sont touchées et coulées.')
            for i in self.machine_position:
                if shoot in self.machine_position[i]:
                    self.machine_position[i].pop(self.machine_position[i].index(shoot))
                    if len(self.machine_position[i]) == 0:
                        print('le tir  %s%s a coulé mon %s!'%(chr(shoot[0]+65), str(shoot[1]+1), i))
                    else:
                        print('mon navire est touchée à %s%s!'%(chr(shoot[0]+65), str(shoot[1]+1)))
        else:
            print('Raté!')

    def check_target(self, check, shoot):
        mst = {'R':'Ratée', 'T':'Touchée', 'C':'Coulée'}
        if (check.upper() not in 'RTC') or (check.upper() == ''):
            raise Error('(R)até, (T)ouché or (C)oulé :')
        else:
            confirm = input('confirmer si votre navire est %s! ou 0-zero pour refaire: '%mst[check.upper()])
            if (confirm == '0'):
                raise Error('confirmer si(R)até, (T)ouché or (C)oulé: ')
            else:
                self.player_target.pop(self.player_target.index(shoot))
                if check.upper() == 'R': ## should not be necessary
                    self.player_missed.append(shoot)
                elif check.upper() =='T':
                    self.player_touched.append(shoot)
                elif check.upper() == 'C':
                    self.player_touched.append(shoot)
                    self.player_sank.extend(self.player_touched)
                    win = sum(NAVIRES.values())
                    if len(self.player_sank) == win:
                        raise GameOver('uh! toutes vos navires sont coulées.')
                    a,b,c,d = self.check_around(shoot, max(NAVIRES.values()))
                    p = [i for i in a+b+c+d if i in self.player_touched]
                    p.append(shoot)
                    for ps in p:
                        ap = self.check_around_position(ps)
                        for a in ap:
                            if a in self.player_target:
                                self.player_target.pop(self.player_target.index(a))
                    self.player_touched = []

    def pick_player_target(self):
        if self.player_touched:
            s = []
            for p in self.player_touched:
                ap = self.check_around_position(p, touched=True)
                for a in ap:
                    if a in self.player_target:
                        s.append(a)
            if len(self.player_touched) == 1:
                return random.choice(s)
            else:
                rs = []
                t0 = self.player_touched[0]
                t1 = self.player_touched[1]
                if t0[0] == t1[0]:
                    for i in s:
                        if i[0] == t0[0]:
                            rs.append(i)
                        else:
                            pass
                elif t0[1] == t1[1]:
                    for i in s:
                        if i[1] == t0[1]:
                            rs.append(i)
                        else:
                            pass
                return random.choice(rs)
        else:
            return random.choice(self.player_target)

if __name__ == '__main__':
    game = Setup(size=TAILLE)
    print('la bataille commence ...')
    #print(game.machine_touched)
    if machine_starts:
        turn = 0
    else:
        turn = 1
    while True:
        try:
            if (turn % 2) == 0:
                print('-- la machine joue --')
                shoot = game.pick_player_target() #random.choice(game.player_target)
                go = False
                while not go:
                    try:
                        check = input('je tire sur %s%s! (R)até, (T)ouché or (C)oulé? '%(chr(shoot[0]+65),str(shoot[1]+1)))
                        game.check_target(check, shoot)
                        go = True
                    except Error as e:
                        print(e)
                turn += 1                
            else:
                print('-- vous jouez --')
                go = False
                while not go:
                    try:
                        shoot = input('entrez votre tir: ')
                        game.shoot_target(shoot)
                        go = True
                    except Error as e:
                        print(e)
                turn += 1
        except GameOver as e:
            print(e)
            print('la position de mes navires:')
            print(game.where)
            break
