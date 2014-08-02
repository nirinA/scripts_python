'''<h1>Jeu de Mangorona</h1>
<h2>nouveau jeu , new game</h2>
<h3>configuration</h3>
on peut configurer les &eacute;l&eacute;ments suivants
pour un nouveau jeu:
<h4>lattice: r&eacute;seau</h4>
<ul>
<li>star: le r&eacute;seau habituel du jeu
<li>diamond: idem &agrave; star mais invers&eacute;.
mais moins int&eacute;ressant car les coins sont &quot;morts&quot;
<li>cubic: on joue seulement sur les lignes horizontales et
verticales.<!-- int&eacute;ressant pour une grande dimension et
avec les configurations &quot;maximun&quot; sur &quot;no&quot; et &quot;all&quot;
sur &quot;yes&quote&. -->
<li>web: plus de lignes sur le plateau mais il y a des endroits o&ugrave; 
on ne peut pas se placer (les intersections &agrave; l'int&eacute;rieur
des cases)
<li>X: idem &agrave; web mais sans les lignes horizontales et verticales
</ul>
les deux derniers <!-- ne sont pas vraiment jouables; mais --> ont &eacute;t&eacute;
rajout&eacute; car ce sont seulement 2 lignes de codes dans le programme entier.
<h4>dimension</h4>
comme son nom l'indique.<br>
on peut rajouter des dimensions, mais cela peut ne pas &ecirc;tre
jouable si la dimension est grande et l'&eacute;cran petit.
<br>pour ajouter les dimensions: dans le code, localiser les
lignes suivantes:
<pre>
    list_dimension = ['5 x 7',
                      '5 x 9',
                      '7 x 9',
                      '3 x 5',
                      '3 x 7',
                      '5 x 11',
                      '7 x 11',
                      '9 x 11']

</pre>
rajouter les dimension voulues, par ex:
<pre>
    list_dimension = ['5 x 7',
                      '5 x 9',
                      '7 x 9',
                      '3 x 5',
                      '3 x 7',
                      '5 x 11',
                      '7 x 11',
                      '9 x 11',
                      '3 x 9',
                      '3 x 11',
                      '5 x 5',]

</pre>
ici on a ajout&eacute;:
<pre>
                      '3 x 9',
                      '3 x 11',
                      '5 x 5',
</pre>
enregister et jouer...
<br>
<b>note: les dimensions doivent &ecirc;tre des nombres impairs.</b>
<h4>maximum</h4>
AI de la machine.<br>
si sur &quot;yes&quot; la machine essaie de prendre
 le maximum de pi&egrave;ces. sinon, elle joue un peu n'importe comment.
<h4>all</h4>
si sur &quot;yes&quot; on peut jouer tant qu'il y ait des pi&egrave;ces
sur le plateau, tant qu'on peut se mouvoir. sinon, le jeu s'arr&ecirc;tte
 d&egrave;s que l'on ne peut plus prendre de pi&egrave;ces.
<h4>who ?</h4>
qui commence? <br>les rouges.
<hr>
les configurations par d&eacute;faut:<br>
<ul>
<li>lattice:star
<li>dimension:5 x 9
<li>maximum:yes
<li>all:no
<li>who:machine
</ul>
<h2>charger jeu , load game</h2>
on peut charger les jeux enregistr&eacute;s avec les 2 programmes,
qu'importe le nom, mais avec les extensions .mang pour pouvoir les
rechercher facilement.
<h2>enregister jeu , save game</h2>
idem pour charger, le nom n'a pas d'importance
<h2>quit</h2>
no comment
<h2>aide</h2>
cette page
'''

import tkinter
import time
import random
import traceback
import webbrowser

import tkinter.scrolledtext as ScrolledText
from tkinter.filedialog import asksaveasfilename, askopenfilename
import pickle

from mangorona import *
from mangorona import __doc__ as maindoc

global root, _trace, w_frame

class TkMangorona(Mangorona, tkinter.Frame):
    def __init__(self, master, players, lattice, dimension):
        tkinter.Frame.__init__(self, master)
        Mangorona.__init__(self, players, lattice, dimension)
        self.frame = tkinter.Frame(self)
        self.frame.pack()
        self.framelog = tkinter.Frame(self)
        self.framelog.pack(side=tkinter.BOTTOM)
        self.textlog = ScrolledText.ScrolledText(self.framelog, height=5)
        self.textlog.pack()

        self.canvas = tkinter.Canvas(self.frame, bg='white')
        self.canvas.pack()
        self.player = {self.player1:'red', self.player2:'blue', self.blank:None}
        self.height = self.xmax*50
        self.width = self.ymax*50
        self.canvas['height'] = self.height
        self.canvas['width'] = self.width+50
        self.coords = [(x*50+10,y*50+10,x*50+26,y*50+26) for (y,x) in self.all]
        self.tag = ['(%i,%i)'%(x,y) for (x,y) in self.all]
        for (p,w,h) in [(self.player1,self.width+10, 50),\
                        (self.player2,self.width+10, self.height-50)]:
            self.DisplayGain(p,w,h)
        self.Board()

    def DisplayGain(self, p, w, h):
        self.canvas.create_text(w, h, text='0', fill=self.player[p],
                                font=('arial', 40, 'bold'), tags=p)
    def Board(self):
        if self.lattice != 'X':
            for i in range(self.xmax):
                self.canvas.create_line(18,18+i*50,self.width-32,18+i*50)
            for i in range(self.ymax):
                self.canvas.create_line(18+i*50,18,18+i*50,self.height-32)
        for (x,y) in self.all:
            a = [i for i in Position((x,y), self.dimension, self.lattice).around if i in self.all]
            for (z,t) in a:
                self.canvas.create_line(list(map(lambda i:18+i*50, (y,x)+(t,z))))

    def Table(self, m):
        for i in self.tag:
            self.canvas.delete(i)
        ##self.canvas.delete('all')
        self.pawns = [self.player[m.matrix[x][y]] for (x,y) in self.all]
        self.oval = [self.canvas.create_oval(xy,fill=c,outline=c,tags=t)
                     for (xy,c,t) in zip(self.coords,self.pawns,self.tag)
                     if c is not None]
        for p in [self.player1, self.player2]:
            self.canvas.itemconfig(p, text=str(m.gain[p]))
        for i in [self.canvas.find_withtag(t) for t in self.oval]:
            self.canvas.tag_bind(i, '<Button-1>', self.Start)
            self.canvas.tag_bind(i, '<B1-Motion>', self.Move)
            self.canvas.tag_bind(i, '<ButtonRelease-1>', self.Destination)
        self.frame.update()

    moving = False
    initial = None
    final = None

    def Start(self, event):
        self.finishmoving()
        self.nowtag = self.canvas.find_withtag(tkinter.CURRENT)
        self.taginitial = self.canvas.gettags(self.nowtag)
        self.initial = eval(self.taginitial[0])
        self.tagcolor = self.canvas.itemcget(self.nowtag, 'fill')
        self.canvas.itemconfig(self.nowtag, fill='green')
        self.startmoving(event)

    def Move(self, event):
        self.keepmoving(event)

    def Destination(self, event):
        self.canvas.itemconfig(self.nowtag, fill=self.tagcolor)
        self.keepmoving(event)
        self.finishmoving()

    def startmoving(self, event):
        self.moving = False
        self.ix = self.lastx = event.x
        self.iy = self.lasty = event.y
        self.moving = True

    def keepmoving(self, event):
        if not self.moving:
            return
        self.canvas.move(self.nowtag, event.x - self.lastx, event.y - self.lasty)
        self.lastx = event.x
        self.lasty = event.y
        xix, yiy = self.initial
        xfx = xix + int(self.lasty/50)-int(self.iy/50)
        yfy = yiy + int(self.lastx/50)-int(self.ix/50)
        self.final = (xfx,yfy)

    def finishmoving(self):
        self.moving = False
        move = (self.initial, self.final)
        if all(move):
            return move

def TkGamevsMachine(root, players, lattice, dimension, maximum, getall, matrix,
                    whostart=0, tab=0, gain=None):
    global _trace
    machineplayer = players[whostart]
    mc = Mangorona(players, lattice, dimension, matrix=matrix)
    mat = TkMangorona(root, players, lattice, dimension)
    mat.pack()
    mat.Table(mc)
    if gain:
        mc.gain = gain
##    tab = tab
    _trace = [(players, lattice, dimension, maximum, getall, mc.matrix,
               whostart, tab, mc.gain)]
    while True:
        try:
            turn = players[tab%2]
            movable = AllowableMovement(mc, turn).Move(maximum=maximum, getall=getall)
            if turn == machineplayer:
                machine = random.choice(movable)
                mat.textlog.insert(tkinter.END,
                                   '%s move:%s to %s\n'%(mat.player[turn], str(machine[0]), str(machine[1])))
                mat.textlog.see(tkinter.END)
                mc.Move(turn, machine[0], machine[1])
                tab += 1
                mat.Table(mc)
                _trace.append((players, lattice, dimension, maximum, getall, mc.matrix,
                               whostart, tab, mc.gain))
            else:
                while not mat.moving:
                    mat.Table(mc)
                    time.sleep(1)
                    if mat.moving:
                        human = mat.finishmoving()
                        break
                mc.Move(turn, human[0], human[1])
                mat.textlog.insert(tkinter.END,
                                   '%s move:%s to %s\n'%(mat.player[turn], str(human[0]), str(human[1])))
                mat.textlog.see(tkinter.END)
                mat.Table(mc)
                tab += 1
                _trace.append((players, lattice, dimension, maximum, getall, mc.matrix,
                               whostart, tab, mc.gain))
        except NoMoreMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            mat.textlog.insert(tkinter.END, exc)
            mat.textlog.insert(tkinter.END, 'winner:%s'%mat.player[mc.Winner()])
            mat.textlog.see(tkinter.END)
            break
        except IllegalMove:
            exc = traceback.format_exception(*sys.exc_info())[-1]
            mat.textlog.insert(tkinter.END, exc)
        except KeyboardInterrupt:
            raise SystemExit
        except:
            traceback.print_exc()
            raise SystemExit

def TkTest():
    global w, root, _trace
    root = tkinter.Tk()
    root.title('MANGORONA')
    menubar = tkinter.Menu(root)

    gamemenu =  tkinter.Menu(menubar)
    gamemenu.add_command(label="new game", command=new_game)
    gamemenu.add_command(label="load game", command=load_game)
    gamemenu.add_command(label="save game", command=save_game)
    gamemenu.add_separator()
    gamemenu.add_command(label="quit", command=root.destroy)

    helpmenu = tkinter.Menu(menubar, name='help')
    helpmenu.add_command(label="aide...blah blah", command=aide)

    menubar.add_cascade(label="game", menu=gamemenu)
    menubar.add_cascade(label="Aide", menu=helpmenu)

    root['menu']=menubar
    root.mainloop()

def new_game():
    global root, w_frame
    w_frame = tkinter.Frame(root)
    w_frame.pack()
    tkinter.Label(w_frame, text='configuration').pack(side=tkinter.TOP)
    bp = tkinter.Button(w_frame, text="play")
    bp.pack(side=tkinter.BOTTOM)

    ## lattice
    menu_lattice = tkinter.Frame(w_frame)
    menu_lattice.pack(side=tkinter.LEFT)
    tkinter.Label(menu_lattice, text='lattice',
                  bg='green', fg='blue').pack(anchor=tkinter.N)
    scrollbar_l = tkinter.Scrollbar(menu_lattice,
                                    orient=tkinter.VERTICAL)
    listbox_l = tkinter.Listbox(menu_lattice,height=3, width=10,
                                yscrollcommand=scrollbar_l.set,
                                selectmode=tkinter.SINGLE,
                                exportselection=0)
    scrollbar_l.config(command=listbox_l.yview)
    scrollbar_l.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    listbox_l.pack(side=tkinter.LEFT,
                   fill=tkinter.BOTH, expand=tkinter.YES)
    list_lattice = ['star', 'diamond', 'cubic', 'web', 'X']
    for i in list_lattice:
        listbox_l.insert(tkinter.END, i)
    listbox_l.select_set(0)

    ## dimension
    menu_dimension = tkinter.Frame(w_frame)
    menu_dimension.pack(side=tkinter.LEFT)
    tkinter.Label(menu_dimension, text='dimension',
                  bg='green', fg='blue').pack(anchor=tkinter.N)
    scrollbar = tkinter.Scrollbar(menu_dimension,
                                  orient=tkinter.VERTICAL)
    listbox = tkinter.Listbox(menu_dimension,height=3, width=10,
                              yscrollcommand=scrollbar.set,
                              selectmode=tkinter.SINGLE,
                              exportselection=0)
    scrollbar.config(command=listbox.yview)
    scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    listbox.pack(side=tkinter.LEFT,
                 fill=tkinter.BOTH, expand=tkinter.YES)
    list_dimension = ['5 x 7',
                      '5 x 9',
                      '7 x 9',
                      '3 x 5',
                      '3 x 7',
                      '5 x 11',
                      '7 x 11',
                      '9 x 11']
    for i in list_dimension:
        listbox.insert(tkinter.END, i)
    listbox.select_set(1)
    
    ## machine rule
    menu_machine_rule = tkinter.Frame(w_frame)
    menu_machine_rule.pack(side=tkinter.LEFT)
    tkinter.Label(menu_machine_rule, text='maximum',
                  bg='green', fg='blue').pack(anchor=tkinter.N)
    machine_rule = tkinter.StringVar()
    machine_rule.set('yes')
    for t in ['yes', 'no']:
        tkinter.Radiobutton(menu_machine_rule, text=t,
                            variable=machine_rule,
                            value=t).pack(anchor=tkinter.W)
    
    ## game rule
    menu_game_rule = tkinter.Frame(w_frame)
    menu_game_rule.pack(side=tkinter.LEFT)
    tkinter.Label(menu_game_rule, text='all',
                  bg='green', fg='blue').pack(anchor=tkinter.N)
    game_rule = tkinter.StringVar()
    game_rule.set('no')
    for t in ['yes', 'no']:
        tkinter.Radiobutton(menu_game_rule, text=t,
                            variable=game_rule,
                            value=t).pack(anchor=tkinter.W)

    ## who starts
    menu_starts_rule = tkinter.Frame(w_frame)
    menu_starts_rule.pack(side=tkinter.LEFT)
    tkinter.Label(menu_starts_rule, text='who ?',
                  bg='green', fg='blue').pack(anchor=tkinter.N)
    starts_rule = tkinter.StringVar()
    starts_rule.set('machine')
    for t in ['machine', 'human']:
        tkinter.Radiobutton(menu_starts_rule, text=t,
                            variable=starts_rule,
                            value=t).pack(anchor=tkinter.W)

    TF = {'yes':True,'no':False}
    MH = {'machine':0, 'human':1}
    
    bp['command']=lambda :run(list_lattice[int(listbox_l.curselection()[0])],
                              get_dimension(list_dimension[int(listbox.curselection()[0])]),
                              TF[machine_rule.get()],
                              TF[game_rule.get()],
                              MH[starts_rule.get()])

def get_dimension(d):
    l = d.split('x')
    return int(l[0]), int(l[1])

def load_game():
    global root, _trace
    try:
        filename = askopenfilename(parent=None,
                                     filetypes=[("mangorona", "*.mang"),
                                                ("tous", "*")],
                                     initialfile='save.mang'
                                     )
        with open(filename, 'rb') as o:
            mmt = pickle.load(o)
        players, lattice, dimension, maximum, getall, matrix, whostart, tab, gain = mmt[-1]
        TkGamevsMachine(root, players, lattice, dimension, maximum, getall, matrix,
                        whostart, tab, gain)
    except:
        pass

def run(lattice, dimension, m_r, g_r, w):
    global root, _trace, w_frame
    w_frame.destroy()
    print(lattice, dimension, m_r, g_r, w)
    PLAYERS = 'a','b',' '
    TkGamevsMachine(root, PLAYERS, lattice, dimension, m_r, g_r, None, whostart=w)

def save_game():
    global _trace
    try:
        filename = asksaveasfilename(parent=None,
                                     filetypes=[("mangorona", ".mang"),
                                                ("tous", "*")],
                                     initialfile='save.mang'
                                     )
        with open(filename, 'wb') as s:
            pickle.dump(_trace, s)
    except:
        pass

def aide():
    top = tkinter.Toplevel()
    txt = tkinter.Text(top)
    txt.insert(tkinter.END, maindoc)
    txt.pack()
##    with open('./mangorona.html', 'w') as wbr:
##        wbr.write(__doc__)
##    webbrowser.open('./mangorona.html', new=0)
##    #print(__doc__)

__version__ = '1.0.0'
__author__ = 'nirinA'
__date__ = 'Sun Dec 18 23:01:41 2011'

if __name__ == '__main__':
    TkTest()
