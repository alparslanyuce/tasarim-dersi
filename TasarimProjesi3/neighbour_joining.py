from neighbour_joining_tree import NeighbourJoiningTree
from sequence import Sequence
from ete3 import Tree
from tkinter import *
import tkinter as tk

from tkinter import filedialog, Text
import tkinter.messagebox
import os
import matplotlib
import copy
import numpy as np

from Bio import Phylo


def dosyaAc3():
    os.system("sekanslar.txt")




def hakkinda():
    tkinter.messagebox.showinfo("Hakkında","Merhaba! Bu demo programı, Tasarım dersi için oluşturulmuştur.")

def cikis():

        quit(0)



root = tk.Tk()


root.title("Proje")

menu = Menu(root)
root.config(menu=menu)

subm1 = Menu(menu)
menu.add_cascade(label="Dosya", menu=subm1)
subm1.add_command(label="Çıkış", command=cikis)


subm2 = Menu(menu)
menu.add_cascade(label="Yardım", menu=subm2)
subm2.add_command(label="Hakkında", command = hakkinda)


canvas = tk.Canvas(root, height=500, width=500, bg="#776ec5")
canvas.pack()

frame = tk.Frame(root, bg = "white")
frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)









dosyaAcmaEtiketi = tk.Button(root, text="Metin dosyasini ac", padx=10, pady=5, fg="white", bg="red", command=dosyaAc3).place(x=190, y=300)
uyarix = Label(frame, text = "Dikkat! Metin dosyasındaki değişikliklerin işlenmesi için ", fg = 'black',  font = ('arial', 9, 'bold')).place(x=45, y=285)
uyari2x = Label(frame, text = "programın yeniden başlatılması gerekmektedir.", fg = 'black',  font = ('arial', 9, 'bold')).place(x=60, y=305)

baslik = Label(frame, text = "Neighbor Joining Hesaplama", fg = 'green',  font = ('arial', 18, 'bold'))
baslik.pack()
bosluk = Label(frame, text = "                               ", fg = 'magenta',  font = ('arial', 18, 'bold'))
bosluk.pack()
dosya = open('sekanslar.txt')
yazdirma = dosya.read()


yaz = Label(frame, text = yazdirma, fg = 'black',  font = ('arial', 7, 'normal'))
yaz.pack()

uyari3x = Label(frame, text = "Dikkat! Yukarıdakı sekansların hepsi ekrana sığmayabilir. ", fg = 'black',  font = ('arial', 9, 'bold')).place(x=45, y=200)
uyari3x = Label(frame, text = "Tam halini görmek için metin belgesini açınız.", fg = 'black',  font = ('arial', 9, 'bold')).place(x=45, y=220)



class NeighbourJoining:
    def __init__(self, sekanslar):
        self._tree = NeighbourJoiningTree()
        self._tree.normalize(sekanslar)
        self._tree.calculate_matrix_distances()
        self._nodes = []

    def _size(self):
        return len(self._tree.distances[0])

    def Hesaplama(self):
        n = 0
        mapping = [i for i in range(0, self._size())]
        while n <= self._size() - 2:
            r = []
            m = [[0] * self._size() for i in range(self._size())]
            nodes = []

            for i in range(self._size()):
                result = 0
                for j in range(i):
                    result += self._tree.distances[i][j]

                for k in range(i + 1, self._size()):
                    result += self._tree.distances[k][i]

                r.append(result)

            for i in range(1, self._size()):
                for j in range(i):
                    m[i][j] = self._tree.distances[i][j] - (r[i] + r[j]) / (self._size() - 2)

            minimum = m[0][0]
            minimum_index = [0, 0]

            for i in range(1, self._size()):
                for j in range(i):
                    if m[i][j] < minimum:
                        minimum = m[i][j]
                        minimum_index[0] = i
                        minimum_index[1] = j

            p1 = self._tree.distances[minimum_index[0]][minimum_index[1]]

            sum_u1 = p1 / 2 + (r[0] - r[1]) / (2 * (self._size() - 2))
            sum_u2 = p1 - sum_u1

            uzakliklar = [sum_u1, sum_u2]
            pozisyonlar = [minimum_index[0], minimum_index[1]]
            node = (uzakliklar, pozisyonlar)

            distance_u = []
            for i in range(self._size()):
                if i != minimum_index[0] and i != minimum_index[1]:
                    d1 = self._tree.distances[i][minimum_index[1]]
                    d2 = self._tree.distances[i][minimum_index[0]]
                    distance_u.append((d1 + d2 - self._tree.distances[minimum_index[0]][minimum_index[1]]) / 2)

            new_distances = [[0] * (self._size() - 1) for i in range(self._size() - 1)]

            column_u = [[0] * (self._size() - 1) for i in range(1)]

            for i in range(len(distance_u)):
                new_distances[i + 1][0] = distance_u[i]

            if len(distance_u) > 2:
                for i in range(1, len(new_distances)):
                    column_u[i][0] = distance_u[i - 1]

            new_distances = []

            for i in range(1, self._size()):
                dist = list()
                dist.append(column_u[0][i - 1])
                for j in range(self._size()):
                    if j != minimum_index[0] and j != minimum_index[1]:
                        dist.append(self._tree.distances[i][j])
                new_distances.append(dist)

            self._nodes.append([mapping[minimum_index[1]], mapping[minimum_index[0]]])

            n += 1

            new_mapping = [n * (-1), ]
            for i in range(self._size()):
                if i != minimum_index[0] and i != minimum_index[1]:
                    new_mapping.append(mapping[i])

            mapping = new_mapping
            self._tree.distances = new_distances

        self._nodes.append(mapping)

        return self._nodes
   
  
def main():
    window = Toplevel()
    
    
    window.title('Sonuçlar')
    window.geometry('500x500')
    output_file = 'output.png'
    sequence_filename = 'sekanslar.txt'
    (names, sekanslar) = Sequence.read(sequence_filename)

    neighbour_joining = NeighbourJoining(sekanslar)
    nodes = neighbour_joining.Hesaplama()

    neighbour_joining_tree = []

    for node in nodes:
        neighbour_joining_tree.append([])
        for index in node:
            if index >= 0:
                neighbour_joining_tree[-1].append(names[index])
            else:
                neighbour_joining_tree[-1].append(index)

    tree_newick_str = []

    for i in range(len(neighbour_joining_tree)):
        if nodes[i][0] < 0 and nodes[i][1] < 0:
            tree_newick_str.append("({},{})".format(tree_newick_str[neighbour_joining_tree[i][0] * (-1) - 1],
                                             tree_newick_str[neighbour_joining_tree[i][1] * (-1) - 1]))
        elif nodes[i][1] < 0:
            tree_newick_str.append("({},{})".format(neighbour_joining_tree[i][0],
                                             tree_newick_str[neighbour_joining_tree[i][1] * (-1) - 1]))
        elif nodes[i][0] < 0:
            tree_newick_str.append("({},{})".format(tree_newick_str[neighbour_joining_tree[i][0] * (-1) - 1],
                                             neighbour_joining_tree[i][1]))
        else:
            tree_newick_str.append("({},{})".format(neighbour_joining_tree[i][0], neighbour_joining_tree[i][1]))

    tree_newick_str = "{};".format(tree_newick_str[len(tree_newick_str) - 1])

    tree = Tree(tree_newick_str)
    tree.render(output_file, dpi=600)
    
    canvas = Canvas(window, width=300, height=300)
    canvas.pack()
    image = PhotoImage(file=output_file)
    canvas.create_image(20, 20, anchor=NW, image=image)
    window.mainloop()
    
    
  



        
calistir = tk.Button(root, text="Neighbor Joining sonucunu görmek için tıklayınız", padx=10, pady=5, fg="white", bg="#734a12", command=main)
calistir.pack()
 
        
cikis = tk.Button(root, text="Çıkış yap", padx=10, pady=5, fg="white", bg="#434343", command=cikis)
cikis.pack()

root.mainloop()

