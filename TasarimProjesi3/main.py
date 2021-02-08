import tkinter as tk
from tkinter import filedialog, Text
from tkinter import *
import tkinter.messagebox
import os
import random
import matplotlib
import copy
import numpy as np

from Bio import Phylo
from io import StringIO

from neighbour_joining_tree import NeighbourJoiningTree
from sequence import Sequence
from ete3 import Tree

import sys
import statistics
from binary_tree import Node, pprint, stringify, _bst_insert, _new_node, _validate_tree, _build_tree, _weight_of, \
    _add_left, _add_right, _build_list, _right_of, _left_of, _value_of, _null

match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)
pairwise_alignment_score = 0
sequence = []
names = []
finalSequence = []
score_sequence = []
deltaMatrix = []
fitchIndex = []
filename = "sekanslar3.txt"


def hakkinda():
    tkinter.messagebox.showinfo("Hakkında", "Merhaba! Bu demo programı, Tasarım dersi için oluşturulmuştur.")


def cikis():
    quit(0)


def cikis2():
    quit(0)


def dosyaAc():
    os.system("text.txt")


def dosyaAc2():
    os.system("Sekanslar.txt")


def dosyaAc3():
    os.system("sekanslar2.txt")


def dosyaAc4():
    os.system("sekanslar3.txt")


def UPGMAPenceresi():
    window = Tk()
    window.title("UPGMA")
    menu = Menu(window)
    window.config(menu=menu)

    def geriDon():
        window.destroy()

    subm1 = Menu(menu)
    menu.add_cascade(label="Dosya", menu=subm1)
    subm1.add_command(label="Çıkış", command=cikis2)

    subm2 = Menu(menu)
    menu.add_cascade(label="Yardım", menu=subm2)
    subm2.add_command(label="Hakkında", command=hakkinda)

    canvas = tk.Canvas(window, height=500, width=500, bg="#263D42")
    canvas.pack()

    frame = tk.Frame(window, bg="white")
    frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

    hesaplama = tk.Button(window, text="Hesaplamaları Gör", padx=10, pady=5, fg="white", bg="#263D42",
                          command=Hesaplamalar)
    hesaplama.pack()

    calistir = tk.Button(window, text="UPGMA sonucunu görmek için tıkla", padx=10, pady=5, fg="white", bg="#263D42",
                         command=UPGMAHesabi)
    calistir.pack()

    dosyaAcmaEtiketi = tk.Button(window, text="Metin dosyasini ac", padx=10, pady=5, fg="white", bg="red",
                                 command=dosyaAc).place(x=190, y=290)
    uyarix = Label(frame, text="Dikkat! Metin dosyasındaki değişikliklerin işlenmesi için ", fg='black',
                   font=('arial', 9, 'bold')).place(x=45, y=300)
    uyari2x = Label(frame, text="programın yeniden başlatılması gerekmektedir.", fg='black',
                    font=('arial', 9, 'bold')).place(x=60, y=320)

    baslik = Label(frame, text="UPGMA Hesaplama", fg='magenta', font=('arial', 18, 'bold'))
    baslik.pack()
    cik1 = tk.Button(window, text="Çık", padx=10, pady=5, fg="white", bg="#263D42", command=geriDon)
    cik1.pack()
    bosluk = Label(frame, text="                               ", fg='magenta', font=('arial', 18, 'bold'))
    bosluk.pack()
    label2 = tk.Label(frame, text="Dizinler:", bg="white")
    label2.pack()

    label3 = tk.Label(frame, text="A: " + A, bg="white")
    label3.pack()
    label4 = tk.Label(frame, text="B: " + B, bg="white")
    label4.pack()
    label5 = tk.Label(frame, text="C: " + C, bg="white")
    label5.pack()
    label6 = tk.Label(frame, text="D: " + D, bg="white")
    label6.pack()
    label7 = tk.Label(frame, text="E: " + E, bg="white")
    label7.pack()
    label8 = tk.Label(frame, text="F: " + F, bg="white")
    label8.pack()


def UPGMAHesabi():
    window = Tk()
    window.title("UPGMA")
    window.geometry('350x250')

    yazi = "UPGMA: ", UPGMA(M, M_etiketler)

    print(yazi)
    label = Label(window, text=yazi, fg='red', bg='yellow', relief="solid", font=('arial', 12, 'bold')).place(x=60,
                                                                                                              y=95)
    Phylo.draw(agac)


def sekansOku(filename):
    sonucListesi = []
    infile = open(filename, 'r')

    satir = infile.readline()
    baslik = satir.rstrip()
    etiket = baslik[1:]

    sekans = ''

    for satir in infile:
        satir = satir.rstrip()

        # Boş satırları yoksay.
        if satir != '':

            if satir[0] == '>':
                sonucListesi.append([etiket, sekans])
                baslik = satir.rstrip()
                etiket = baslik[1:]
                sekans = ''

            else:
                sekans += satir

    infile.close()
    sonucListesi.append([etiket, sekans])
    return sonucListesi


sekansBilgisi = sekansOku('Sekanslar.txt')


def main():
    window = Tk()
    window.title("Sonuçlar")
    window.geometry('600x600')
    labelx = Label(window, text="Hesaplama yapılıyor... ", fg='red', font=('arial', 13, 'bold')).place(x=90, y=100)
    # Dosyadan bilgi oku
    sekansBilgisi = sekansOku('Sekanslar.txt')

    sekansListesi = []
    etiketListesi = []
    # Dosyadaki dizilerinin ve etiketlerinin ayrı bir listesi oluşturulur.
    for i in range(0, len(sekansBilgisi)):
        sekans = sekansBilgisi[i][1]
        sekansListesi.append(sekans)
        etiket = sekansBilgisi[i][0]
        etiketListesi.append(etiket)

    # Dizileri bilgi verici (informative) sitelere atayın.
    sekansListesi = bilgiGoster(sekansListesi)
    n = len(sekansListesi[0])  # Bilgi verici sitelerin sayısı

    # A, B, C ve D sekanslarını ayarlama
    cladeCifti1 = 0  # A-B sınıflarının (clade) sayısını (ve dolayısıyla C-D sınıflarını da) sayar.
    cladeCifti2 = 0  # A-C sınıflarının (clade) sayısını (ve dolayısıyla B-D sınıflarını da) sayar.
    cladeCifti3 = 0  # A-D sınıflarının (clade) sayısını (ve dolayısıyla B-C sınıflarını da) sayar.

    # Önyüklemeyi (bootstrapping) 1000 kez çalıştır.
    denemeSayisi = 1000
    for deneme in range(denemeSayisi):
        ornekSiraListesi = []
        for sekans2 in sekansListesi:  # Her sekans için.
            ornek = ''
            for x in range(n):  # Her bilgi verici site için.
                site = random.randrange(n)
                ornek += sekans2[site]
            ornekSiraListesi.append(ornek)  # Burada çoğaltmadan sonra dizilerin bir listesini yapılır.

        ilkAgac = rastgeleAgacOlustur(ornekSiraListesi)  # Rastgele bir başlangıç ağacı seçer.
        uzunluk = uzunlukHesapla(ilkAgac)  # Toplam dal uzunluğunu bul.
        enIyiAgac = ilkAgac

        # En iyi ağaçta değişiklik olmaksızın ağaçların kaç kez karşılaştırıldığını sayın.
        degisiklikYok = 0
        # 10 keyfi olarak seçilen bir eşiktir, ancak en iyi ağaç ince değerini yansıtmalıdır.
        while (degisiklikYok < 10):
            agacKiyasla = rastgeleAgacOlustur(ornekSiraListesi)  # Yeni bir rastgele ağacı oluştur.
            uzunlukKiyasla = uzunlukHesapla(agacKiyasla)  # Bu ağacın dal uzunluğunu bulun.
            if (uzunlukKiyasla < uzunluk):  # Eğer bu şimdiye kadar bulunan en iyi dal uzunluğuysa,
                uzunluk = uzunlukKiyasla  # bunu karşılaştırılacak yeni değer hasatir getirin.
                enIyiAgac = agacKiyasla  # Ardından bunu şimdiye kadarki en iyi değer olarak sakla.
            else:
                degisiklikYok += 1

        # İlk üç deneme için bu topolojileri ekrana yazdırın.
        if (deneme == 0 or deneme == 1 or deneme == 2):

            agacYazdir = ((enIyiAgac[0][0][0], enIyiAgac[0][1][0]),
                          (enIyiAgac[1][0][0], enIyiAgac[1][1][0]))
            # Denemeye göre baskı ön ekini değiştirin.
            if (deneme == 0):
                basKisim = 'Birinci'
            elif (deneme == 1):
                basKisim = 'İkinci'
            elif (deneme == 2):
                basKisim = 'Ucuncu'
            print(basKisim + ' bootstrap (ön yükleme) hesaplaması aşağıdaki topolojiyi ortaya koymaktadır: ')
            print(agacYazdir)
            print('')

            # Birinci sınıftaki (first clade) dizilerin etiketlerini belirleyin.
        # (Örnek: 0 ve 1 değerleri, birinci dizinin ve ikinci dizinin bir sınıfını (a clade) temsil eder.)
        etiket1 = enIyiAgac[0][0][1]
        etiket2 = enIyiAgac[0][1][1]
        # Bir örnek topolojide sınıflar (clades) gözlemlendiğinde sınıf sayaçlarını artırın.
        if ((etiket1 == 0 and etiket2 == 1) or
                (etiket1 == 1 and etiket2 == 0)):
            cladeCifti1 += 1
        elif ((etiket1 == 0 and etiket2 == 2) or
              (etiket1 == 2 and etiket2 == 0)):
            cladeCifti2 += 1
        elif ((etiket1 == 0 and etiket2 == 3) or
              (etiket1 == 3 and etiket2 == 0)):
            cladeCifti3 += 1
        elif ((etiket1 == 1 and etiket2 == 2) or
              (etiket1 == 2 and etiket2 == 1)):
            cladeCifti3 += 1
        elif ((etiket1 == 1 and etiket2 == 3) or
              (etiket1 == 3 and etiket2 == 1)):
            cladeCifti2 += 1
        elif ((etiket1 == 2 and etiket2 == 3) or
              (etiket1 == 3 and etiket2 == 2)):
            cladeCifti1 += 1

    # En çok hangi sınıf (clade) çiftinin göründüğünü belirleyin.
    consensus = max(cladeCifti1, cladeCifti2, cladeCifti3)
    # consensusAgaci'na yazdırılmak üzere kaldırılan etiketler.
    # clade (sınıf) sayaçları, bootstrap (önyükleme) değerini elde etmek için toplam deneme sayısına bölünür.
    if (consensus == cladeCifti1):
        consensusAgaci = [[etiketListesi[0], etiketListesi[1]],
                          [etiketListesi[2], etiketListesi[3]]]
        bootstrap = cladeCifti1 / denemeSayisi
    elif (consensus == cladeCifti2):
        consensusAgaci = [[etiketListesi[0], etiketListesi[2]],
                          [etiketListesi[1], etiketListesi[3]]]
        bootstrap = cladeCifti2 / denemeSayisi
    elif (consensus == cladeCifti3):
        consensusAgaci = [[etiketListesi[0], etiketListesi[3]],
                          [etiketListesi[1], etiketListesi[2]]]
        bootstrap = cladeCifti3 / denemeSayisi
    # Consensus ağacını ve önyükleme değerini yazdırma:
    print('Verilen diziler için consensus ağacı şu şekilde bulundu: ')
    labely = Label(window, text="Verilen diziler için consensus ağacı\n     şu şekilde bulundu: ", fg='red',
                   font=('arial', 13, 'bold')).place(x=90, y=100)
    print('')

    print(str(consensusAgaci[0][0]) + '         ' + str(consensusAgaci[1][0]))
    labelz = Label(window, text=str(consensusAgaci[0][0]) + "         " + str(consensusAgaci[1][0]), fg='magenta',
                   font=('arial', 13, 'bold')).place(x=90, y=155)
    print('                 \________/')
    labelt = Label(window, text="                 \________/", fg='magenta', font=('arial', 13, 'bold')).place(x=90,
                                                                                                               y=180)
    print('                 /        \\')
    labela = Label(window, text="                 /               \\", fg='magenta', font=('arial', 13, 'bold')).place(
        x=90, y=200)
    print(str(consensusAgaci[0][1]) + '         ' + str(consensusAgaci[1][1]))
    labelb = Label(window, text=str(consensusAgaci[0][1]) + "         " + str(consensusAgaci[1][1]), fg='magenta',
                   font=('arial', 13, 'bold')).place(x=90, y=220)
    print('')
    print('bootstrap (önyükleme) değeri: ' + str(bootstrap) + '.')
    labelc = Label(window, text="bootstrap (önyükleme) değeri: " + str(bootstrap) + ".", fg='blue',
                   font=('arial', 13, 'bold')).place(x=90, y=300)


# bilgiGoster: Bir dizi listesindeki bilgilendirici siteleri kontrol eden bir fonksiyondur.
def bilgiGoster(sekansListesi):
    sekansUzunlugu = len(sekansListesi[0])
    bilgiTakip = []  # Dizilerdeki bilgi verici sitelerin yerlerini takip eder.
    for x in range(sekansUzunlugu):
        sayacA = 0
        sayacC = 0
        sayacG = 0
        sayacT = 0
        # Dizileri yineleyin, ardından her bir nükleotid için görünüm sayısını hesaplayın.
        for sekans in sekansListesi:
            if (sekans[x] == 'A'):
                sayacA += 1
            elif (sekans[x] == 'C'):
                sayacC += 1
            elif (sekans[x] == 'G'):
                sayacG += 1
            elif (sekans[x] == 'T'):
                sayacT += 1
        # Bilgi verici bir site, en az iki farklı nükleotidin iki veya daha fazla görüntüsüne sahiptir.
        if ((sayacA >= 2 and sayacC >= 2) or
                (sayacA >= 2 and sayacG >= 2) or
                (sayacA >= 2 and sayacT >= 2) or
                (sayacC >= 2 and sayacG >= 2) or
                (sayacC >= 2 and sayacT >= 2) or
                (sayacG >= 2 and sayacT >= 2)):
            bilgiTakip.append(1)  # "1" ile izleyicide bilgi verici siteleri temsil eder.
        else:
            bilgiTakip.append(0)  # bilgi verici olmayan siteleri 0 ile temsil eder.

    sekansBilgiListe = []

    for sekans in sekansListesi:
        sekansBilgi = ''
        for x in range(sekansUzunlugu):
            if (bilgiTakip[x] == 1):
                sekansBilgi += sekans[x]
        sekansBilgiListe.append(sekansBilgi)

    return sekansBilgiListe


# rastgeleAgacOlustur - Dört diziden oluşan bir listeden rastgele bir ağaç oluşturan bir fonksiyondur.

def rastgeleAgacOlustur(sekansListesi):
    agac = []
    seqnum = 0
    ek = []
    # Sınıfların (clades) daha sonra belirlenmesine yardımcı olmak için dizilerin her birine bir etiket atayın.
    for sekans in sekansListesi:
        etiket = seqnum
        seqnum += 1
        ek.append([sekans, etiket])
    # Diziler için rastgele bir sıra seçin.
    a = random.choice(ek)
    ek.remove(a)
    b = random.choice(ek)
    ek.remove(b)
    c = random.choice(ek)
    ek.remove(c)
    d = random.choice(ek)
    # Rastgele düzeni yansıtan bir ağaç oluşturun.
    agac.append([a, b])
    agac.append([c, d])
    return agac


# uzunlukHesapla - Belirli bir ağacın toplam dal uzunluğunu hesaplayan bir fonksiyondur.

def uzunlukHesapla(agac):
    # Etiketlere değil, sadece dizilere bakılır.
    a = agac[0][0][0]
    b = agac[0][1][0]
    c = agac[1][0][0]
    d = agac[1][1][0]

    # Uzunluk, verilen iki dizi arasındaki ikame sayısı ile belirlenir.
    ab_uzunluk = subsay(a, b)
    cd_uzunluk = subsay(c, d)
    ac_uzunluk = subsay(a, c)
    bd_uzunluk = subsay(b, d)
    # Toplam dal uzunluğunun üst üste binen kısımlarından hesaplanması gerekir.
    l = (ac_uzunluk + bd_uzunluk + ab_uzunluk + cd_uzunluk) / 2

    return l


# subsay - İki dizge arasındaki ikame sayısını (ve dolayısıyla iki dizi için dal uzunluğunu) sayan bir fonksiyondur.

def subsay(seq1, seq2):
    subsayisi = 0
    for x in range(len(seq1)):
        # Bir indeksteki nükleotidler eşleşmediğinde bir ikame meydana gelir.
        if (seq1[x] != seq2[x]):
            subsayisi += 1
    return subsayisi


# sekansOku - FASTA'ya benzer bir biçime sahip dosyaları okuyan bir fonksiyon.


def sekansOku(filename):
    sonucListesi = []
    infile = open(filename, 'r')

    satir = infile.readline()
    baslik = satir.rstrip()
    etiket = baslik[1:]

    sekans = ''

    for satir in infile:
        satir = satir.rstrip()

        # Boş satırları yoksay.
        if satir != '':

            if satir[0] == '>':
                sonucListesi.append([etiket, sekans])
                baslik = satir.rstrip()
                etiket = baslik[1:]
                sekans = ''

            else:
                sekans += satir

    infile.close()
    sonucListesi.append([etiket, sekans])
    return sonucListesi


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


def main2():
    window = Toplevel()

    window.title('Sonuçlar')
    window.geometry('500x500')
    output_file = 'output.png'
    sequence_filename = 'sekanslar2.txt'
    (names, sekanslar) = Sequence.read(sequence_filename)

    def geriDon():
        window.destroy()

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
    cik4 = tk.Button(window, text="Çık", padx=30, pady=8, fg="white", bg="#263D42", command=geriDon)
    cik4.pack()
    window.mainloop()


def MaxParsimonyPenceresi():
    window2 = Tk()
    window2.title("Maksimum Parsimony")
    menu = Menu(window2)
    window2.config(menu=menu)

    def geriDon():
        window2.destroy()

    subm1 = Menu(menu)
    menu.add_cascade(label="Dosya", menu=subm1)
    subm1.add_command(label="Çıkış", command=cikis2)

    subm2 = Menu(menu)
    menu.add_cascade(label="Yardım", menu=subm2)
    subm2.add_command(label="Hakkında", command=hakkinda)

    canvas = tk.Canvas(window2, height=500, width=700, bg="#00a778")
    canvas.pack()

    frame = tk.Frame(window2, bg="white")
    frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)
    label__0 = Label(window2, text="Maksimum Parsimony Hesaplama", fg='magenta', font=('arial', 18, 'bold')).place(
        x=170, y=90)
    label1 = Label(window2, text="Sekansları metin dosyasında görmek için tıklayınız.", fg='red',
                   font=('arial', 13, 'bold')).place(x=170, y=160)

    dosyaAcmaEtiketi = tk.Button(window2, text="Metin dosyasini ac", padx=10, pady=5, fg="white", bg="red",
                                 command=dosyaAc2).place(x=290, y=190)
    uyarix = Label(frame, text="Dikkat! Metin dosyasındaki değişikliklerin işlenmesi için ", fg='black',
                   font=('arial', 9, 'bold')).place(x=125, y=200)
    uyari2x = Label(frame, text="programın yeniden başlatılması gerekmektedir.", fg='black',
                    font=('arial', 9, 'bold')).place(x=140, y=220)

    hesaplama = tk.Button(window2, text="Hesapla", padx=10, pady=5, fg="white", bg="green", command=main).place(x=315,
                                                                                                                y=320)
    uyari = Label(window2, text="Dikkat! Hesaplama süreci 5 ila 10 sn sürmektedir.", fg='black',
                  font=('arial', 9, 'bold')).place(x=220, y=355)

    cik2 = tk.Button(window2, text="Çık", padx=10, pady=5, fg="white", bg="#263D42", command=geriDon)
    cik2.pack()


def NeighborJoiningPenceresi():
    window3 = Tk()
    window3.title("Neighbor Joining")
    menu = Menu(window3)
    window3.config(menu=menu)

    def geriDon():
        window3.destroy()

    subm1 = Menu(menu)
    menu.add_cascade(label="Dosya", menu=subm1)
    subm1.add_command(label="Çıkış", command=cikis2)

    subm2 = Menu(menu)
    menu.add_cascade(label="Yardım", menu=subm2)
    subm2.add_command(label="Hakkında", command=hakkinda)

    canvas = tk.Canvas(window3, height=500, width=500, bg="pink")
    canvas.pack()

    frame = tk.Frame(window3, bg="white")
    frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

    dosyaAcmaEtiketi = tk.Button(window3, text="Metin dosyasini ac", padx=10, pady=5, fg="white", bg="red",
                                 command=dosyaAc3).place(x=190, y=300)
    uyarix = Label(frame, text="Dikkat! Metin dosyasındaki değişikliklerin işlenmesi için ", fg='black',
                   font=('arial', 9, 'bold')).place(x=45, y=285)
    uyari2x = Label(frame, text="programın yeniden başlatılması gerekmektedir.", fg='black',
                    font=('arial', 9, 'bold')).place(x=60, y=305)

    baslik = Label(frame, text="Neighbor Joining Hesaplama", fg='green', font=('arial', 18, 'bold'))
    baslik.pack()
    bosluk = Label(frame, text="                               ", fg='magenta', font=('arial', 18, 'bold'))
    bosluk.pack()
    dosya = open('sekanslar2.txt')
    yazdirma = dosya.read()

    yaz = Label(frame, text=yazdirma, fg='black', font=('arial', 7, 'normal'))
    yaz.pack()

    uyari3x = Label(frame, text="Dikkat! Yukarıdakı sekansların hepsi ekrana sığmayabilir. ", fg='black',
                    font=('arial', 9, 'bold')).place(x=45, y=200)
    uyari3x = Label(frame, text="Tam halini görmek için metin belgesini açınız.", fg='black',
                    font=('arial', 9, 'bold')).place(x=45, y=220)
    calistir = tk.Button(window3, text="Neighbor Joining sonucunu görmek için tıklayınız", padx=10, pady=5, fg="white",
                         bg="#734a12", command=main2)
    calistir.pack()

    cik3 = tk.Button(window3, text="Çık", padx=10, pady=5, fg="white", bg="brown", command=geriDon)
    cik3.pack()


def FitchPenceresi():
    window4 = Tk()
    window4.title("Neighbor Joining")
    menu = Menu(window4)
    window4.config(menu=menu)

    def geriDon():
        window4.destroy()

    subm1 = Menu(menu)
    menu.add_cascade(label="Dosya", menu=subm1)
    subm1.add_command(label="Çıkış", command=cikis2)

    subm2 = Menu(menu)
    menu.add_cascade(label="Yardım", menu=subm2)
    subm2.add_command(label="Hakkında", command=hakkinda)

    canvas = tk.Canvas(window4, height=500, width=600, bg="#ba2b1a")
    canvas.pack()

    frame = tk.Frame(window4, bg="white")
    frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)
    dosyaAcmaEtiketi = tk.Button(window4, text="Metin dosyasini ac", padx=10, pady=5, fg="white", bg="red",
                                 command=dosyaAc4).place(x=240, y=300)
    uyarix = Label(frame, text="Dikkat! Metin dosyasındaki değişikliklerin işlenmesi için ", fg='black',
                   font=('arial', 9, 'bold')).place(x=75, y=285)
    uyari2x = Label(frame, text="programın yeniden başlatılması gerekmektedir.", fg='black',
                    font=('arial', 9, 'bold')).place(x=90, y=305)

    baslik = Label(frame, text="Fitch Algoritması Hesaplama", fg='blue', font=('arial', 18, 'bold'))
    baslik.pack()
    bosluk = Label(frame, text="                               ", fg='black', font=('arial', 18, 'bold'))
    bosluk.pack()
    dosya = open('sekanslar3.txt')
    yazdirma = dosya.read()

    yaz = Label(frame, text=yazdirma, fg='black', font=('arial', 7, 'normal'))
    yaz.pack()

    uyari3x = Label(frame, text="Dikkat! Yukarıdakı sekansların hepsi ekrana sığmayabilir. ", fg='black',
                    font=('arial', 9, 'bold')).place(x=75, y=200)
    uyari3x = Label(frame, text="Tam halini görmek için metin belgesini açınız.", fg='black',
                    font=('arial', 9, 'bold')).place(x=75, y=220)

    calistir = tk.Button(window4, text="Fitch Hesapla", padx=10, pady=5, fg="white", bg="#734a12", command=main3).place(
        x=255, y=390)
    uyari = Label(window4, text="Dikkat! Hesaplama süreci 5 ila 10 sn sürmektedir.", fg='black',
                  font=('arial', 9, 'bold')).place(x=150, y=430)

    cik4 = tk.Button(window4, text="Çıkış yap", padx=10, pady=5, fg="white", bg="#434343", command=geriDon)
    cik4.pack()


def main3():
    # Sequence input structure
    sequence, distance_matrix = fileReader(filename)
    for i in range(1, len(sequence), 2):
        names.append(sequence[i - 1][1:].strip())
        finalSequence.append(sequence[i].strip())

    # Initialize distance_matrix
    distance_matrix = pairwiseDistanceMatrix(finalSequence, distance_matrix)
    score_sequence, min_sequence = scoreSequence(distance_matrix)

    new_root = buildTree(names, score_sequence, finalSequence)
    fitchsIndexCreation(new_root)

    window = Tk()
    window.title("Sonuçlar")
    window.geometry('600x600')

    text_area = Text(window)
    text_area.pack()

    text_area.insert(INSERT, stringify(new_root))

    window.mainloop()


def fileReader(filename):
    # open file with extended ascii
    with open(filename, 'r') as newFile:
        data = newFile.readlines()
        for line in data:
            line.rstrip('\n')
            line = line.replace(' ', " ")  # strip out any whitespace
            sequence.append(line)
    distance_matrix = [
        [0 for col in range(len(sequence) // 2)] for row in range(len(sequence) // 2)]
    return sequence, distance_matrix


def delta(leftSequence, rightSequence):
    deltaScore = 0

    if (leftSequence == rightSequence):
        deltaScore = 0
    else:
        deltaScore = 1
    print(deltaScore)
    # print("Deltas")
    deltaMatrix.append(deltaScore)


def buildTree(names, weights, values):
    medianValues = list(weights)

    medianValues.sort()

    if len(medianValues) % 2 == 0:
        root_value = medianValues[int(len(medianValues) / 2)]
    else:
        root_value = statistics.median(medianValues)

    for x in range(len(weights)):
        if weights[x] == root_value:
            root_index = x
        else:
            continue
    new_root = _new_node(names[root_index], values[root_index], weights[root_index])
    for x in range(len(weights)):
        if x != root_index:
            _bst_insert(new_root, names[x], values[x], weights[x])
        else:
            continue
    return (new_root)


def fitchsIndexCreation(root_node):
    global fitchIndex
    node = root_node
    if root_node == _null:
        return
    if _left_of(node) == _null and _right_of(node) == _null:
        return
    try:
        test_1 = set(_value_of(_left_of(node)))
    except:
        test_1 = set()
    try:
        test_2 = set(_value_of(_right_of(node)))
    except:
        test_2 = set()
    inSet = set.intersection(test_1, test_2)
    if len(inSet) == 0:
        node.value = set.union(test_1, test_2)
    else:
        node.value = inSet

    fitchsIndexCreation(_left_of(node))
    fitchsIndexCreation(_right_of(node))


def compare(a, b):
    notEqual = []
    for x, y in zip(a, b):
        if x == y:
            print('equal')
            notEqual.append((x, y))
        else:
            print('not equal')
            notEqual.append((x, y))
    if len(notEqual) < 1:
        return a
    else:
        return notEqual

        ############################################
        ## Calculates the pairwise distance score ##
        ## for all sequences in S                 ##
        ############################################


def pairwiseDistanceMatrix(sequence_list, distance_matrix):
    global seq1
    global seq2
    global rows
    global cols
    for x in range(len(sequence_list)):
        for y in range(len(sequence_list)):
            seq1 = sequence_list[x]
            seq2 = sequence_list[y]

            rows = len(seq1) + 1
            cols = len(seq2) + 1

            score_matrix, start_pos = createScoreMatrix(rows, cols)
            distance_score = traceback(score_matrix, start_pos)

            distance_matrix[x][y] = distance_score
        # offset += 1
    return distance_matrix


##############################################
## Creates the sequence of alignment scores ##
##############################################

def scoreSequence(d_matrix):
    score_min = 10000000
    score_value = 0
    sc_sequence = []
    for x in range(len(d_matrix)):
        for y in range(len(d_matrix)):
            score_value = score_value + d_matrix[x][y]
        if score_value < score_min:
            score_min_index = x
            score_min = score_value
        sc_sequence.append(score_value)
        score_value = 0
    return sc_sequence, score_min_index


###############################################
## Creates scoring matrix from input strings ##
###############################################


def createScoreMatrix(rows, cols):
    global maxScore
    maxPosition = (0, 0)
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            similarity = match if seq1[i - 1] == seq2[j - 1] else other
            diag_score = score_matrix[i - 1][j - 1] + similarity
            up_score = score_matrix[i - 1][j] + other
            left_score = score_matrix[i][j - 1] + other
            curMax = max(0, diag_score, up_score, left_score)
            if curMax > maxScore:
                maxScore = curMax
                maxPosition = (i, j)
            score_matrix[i][j] = curMax
    maxScore = 0
    return score_matrix, maxPosition


#######################################
## Creates best fit alignment string ##
## based on scoring matrix           ##
#######################################


def traceback(score_matrix, start_pos):
    pairwise_alignment_score = 0
    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = start_pos
    move = nextMove(score_matrix, x, y)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
            pairwise_alignment_score += 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1
            pairwise_alignment_score += 1

        move = nextMove(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq1[y - 1])

    return pairwise_alignment_score


################################################
## Determines the next move for the traceback ##
################################################


def nextMove(score_matrix, x, y):
    # Assign the diagonal score
    diag = score_matrix[x - 1][y - 1]

    # Assign insertion/deletion scores
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]

    # Check all three cases to find next character/insertion/deletion
    if diag >= up and diag >= left:
        return 1 if diag != 0 else 0
    elif up > diag and up >= left:
        return 2 if up != 0 else 0
    elif left > diag and left > up:
        return 3 if left != 0 else 0

    # Error detection
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')


###############################################
## Creates the alignment string for printing ##
###############################################


def alignment_string(aligned_seq1, aligned_seq2):
    # Sets initial values
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []

    # Runs through both strings
    for base1, base2 in zip(aligned_seq1, aligned_seq2):

        # Checks for match
        if base1 == base2:
            alignment_string.append('|')
            idents += 1

        # Checks for insertion/deletion
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1

        # If neither of the above, it's mismatch
        else:
            alignment_string.append(':')
            mismatches += 1

    # Returns the "alignment" string and the alignment characteristics
    return ''.join(alignment_string), idents, gaps, mismatches


def Hesaplamalar():
    window = Tk()
    window.title("Sonuçlar")
    window.geometry('450x350')
    yazi = "\nA ile B uzaklığı: ", uzaklikAB, "\nA ile C uzaklığı: ", uzaklikAC, "\nA ile D uzaklığı: ", uzaklikAD, "\nA ile E uzaklığı: ", uzaklikAE, "\nA ile F uzaklığı", uzaklikAF, "\nB ile C uzaklığı: ", uzaklikBC, "\nB ile D uzaklığı: ", uzaklikBD, "\nB ile E uzaklığı: ", uzaklikBE, "\nB ile F uzaklığı: ", uzaklikBF, "\nC ile D uzaklığı: ", uzaklikCD, "\nC ile E uzaklığı: ", uzaklikCE, "\nC ile F uzaklığı: ", uzaklikCF, "\nD ile E uzaklığı: ", uzaklikDE, "\nD ile F uzaklığı: ", uzaklikDF, "\nE ile F uzaklığı: ", uzaklikEF
    yazi2 = "Sekans başına harf sayısı: ", harfSayisi, "\n"
    label = Label(window, text=yazi2, fg='red', font=('arial', 13, 'bold')).place(x=70, y=10)
    label2 = Label(window, text=yazi, fg='blue', font=('arial', 10, 'normal')).place(x=120, y=35)


def Hesaplamalar2():
    window = Tk()
    window.title("Sonuçlar")
    window.geometry('450x350')
    yazi = "\nA ile B uzaklığı: ", uzaklikAB2, "\nA ile C uzaklığı: ", uzaklikAC2, "\nA ile D uzaklığı: ", uzaklikAD2, "\nA ile E uzaklığı: ", uzaklikAE2, "\nB ile C uzaklığı: ", uzaklikBC2, "\nB ile D uzaklığı: ", uzaklikBD2, "\nB ile E uzaklığı: ", uzaklikBE2, "\nC ile D uzaklığı: ", uzaklikCD2, "\nC ile E uzaklığı: ", uzaklikCE2, "\nD ile E uzaklığı: ", uzaklikDE2
    yazi2 = "Sekans başına harf sayısı: ", harfSayisi2, "\n"
    label = Label(window, text=yazi2, fg='#082567', font=('arial', 13, 'bold')).place(x=70, y=10)
    label2 = Label(window, text=yazi, fg='#321414', font=('arial', 10, 'normal')).place(x=120, y=35)


root = tk.Tk()

root.title("Proje")

menu = Menu(root)
root.config(menu=menu)

subm1 = Menu(menu)
menu.add_cascade(label="Dosya", menu=subm1)
subm1.add_command(label="Çıkış", command=cikis)

subm2 = Menu(menu)
menu.add_cascade(label="Yardım", menu=subm2)
subm2.add_command(label="Hakkında", command=hakkinda)

canvas = tk.Canvas(root, height=500, width=500, bg="green")
canvas.pack()

frame = tk.Frame(root, bg="white")
frame.place(relwidth=0.8, relheight=0.8, relx=0.1, rely=0.1)

# hesaplama = tk.Button(window2, text="Hesapla", padx=10, pady=5, fg="white", bg="green", command=main).place(x=315, y=320)


labelBas = Label(root, text="Tasarım Projesi", fg='#3d4caa', font=('arial', 18, 'bold')).place(x=160, y=90)

upgmaHesabi = tk.Button(root, text="UPGMA", padx=10, pady=5, fg="white", bg="#263D42", font=('arial', 13, 'bold'),
                        command=UPGMAPenceresi).place(x=200, y=180)

maxParsimonyHesabi = tk.Button(root, text="Maksimum Parsimony", padx=10, pady=5, fg="white", bg="#00a778",
                               font=('arial', 13, 'bold'), command=MaxParsimonyPenceresi).place(x=145, y=240)

neighborjoiningHesabi = tk.Button(root, text="Neighbor Joining", padx=10, pady=5, fg="white", bg="pink",
                                  font=('arial', 13, 'bold'), command=NeighborJoiningPenceresi).place(x=165, y=300)

fitchHesabi = tk.Button(root, text="Fitch Algoritması", padx=10, pady=5, fg="white", bg="#ba2b1a",
                        font=('arial', 13, 'bold'), command=FitchPenceresi).place(x=165, y=360)

cikis = tk.Button(root, text="Çıkış yap", padx=10, pady=5, fg="white", bg="red", font=('arial', 9, 'bold'),
                  command=cikis)
cikis.pack()

A = ""
B = ""
C = ""
D = ""
E = ""
F = ""
bayrak1 = 1;
bayrak2 = 0;
bayrak3 = 0;
bayrak4 = 0;
bayrak5 = 0;
bayrak6 = 0;
satir = '';

# Dosya okuma işlemi
dosya = open('text.txt', 'r')
while True:
    char = dosya.read(1)

    # Sırasıyla okumalar başlıyor. Eğer '\n' sembolü gelirse duracaktır, ardından diğer dizilim için işlemler
    # gerçekleştirilecektir.

    if bayrak6 == 1:
        if (char != '\n'):
            F += char
        if (char == ''):
            bayrak6 = 0
            break

    if bayrak5 == 1:
        if (char != '\n'):
            E += char
        if (char == '\n'):
            bayrak5 = 0
            bayrak6 = 1

    if bayrak4 == 1:
        if (char != '\n'):
            D += char
        if (char == '\n'):
            bayrak4 = 0
            bayrak5 = 1

    if bayrak3 == 1:
        if (char != '\n'):
            C += char
        if (char == '\n'):
            bayrak3 = 0
            bayrak4 = 1

    if bayrak2 == 1:
        if (char != '\n'):
            B += char
        if (char == '\n'):
            bayrak2 = 0
            bayrak3 = 1

    if bayrak1 == 1:
        if (char != '\n'):
            A += char

        if (char == '\n'):
            bayrak1 = 0
            bayrak2 = 1

print("Dizilimler: \n")

print(A)
print(B)
print(C)
print(D)
print(E)
print(F)
print("")

harfSayisi = len(A)  # Karakter uzunluğu belirleme işlemi

print("Harf sayisi: ", harfSayisi, "\n")

sayi = 0
sayac = 1
uzaklikAB = 0
uzaklikAC = 0
uzaklikAD = 0
uzaklikAE = 0
uzaklikAF = 0
uzaklikBC = 0
uzaklikBD = 0
uzaklikBE = 0
uzaklikBF = 0
uzaklikCD = 0
uzaklikCE = 0
uzaklikCF = 0
uzaklikDE = 0
uzaklikDF = 0
uzaklikEF = 0
yazi1 = ""
yazi2 = ""
uzaklik = 0

# Uzaklık kontrolleri için satır ve sütunların seçimleri.

while (sayac <= 15):
    if (sayac == 1):
        yazi1 = A
        yazi2 = B
    if (sayac == 2):
        yazi1 = A
        yazi2 = C
    if (sayac == 3):
        yazi1 = A
        yazi2 = D
    if (sayac == 4):
        yazi1 = A
        yazi2 = E
    if (sayac == 5):
        yazi1 = A
        yazi2 = F
    if (sayac == 6):
        yazi1 = B
        yazi2 = C
    if (sayac == 7):
        yazi1 = B
        yazi2 = D
    if (sayac == 8):
        yazi1 = B
        yazi2 = E
    if (sayac == 9):
        yazi1 = B
        yazi2 = F
    if (sayac == 10):
        yazi1 = C
        yazi2 = D
    if (sayac == 11):
        yazi1 = C
        yazi2 = E
    if (sayac == 12):
        yazi1 = C
        yazi2 = F
    if (sayac == 13):
        yazi1 = D
        yazi2 = E
    if (sayac == 14):
        yazi1 = D
        yazi2 = F
    if (sayac == 15):
        yazi1 = E
        yazi2 = F

    # Hangi karakterler farklı ise o halde uzaklık değişkeni birer arttırılacaktır.

    while (sayi < harfSayisi):
        if (yazi1[sayi] != yazi2[sayi]):
            uzaklik += 1
        sayi += 1

        # sayac değişkenine göre satır ve sütun değişkenlerinin seçimi.
    if (sayi == harfSayisi):
        if (sayac == 1):
            uzaklikAB = uzaklik
            print("A ile B uzakligi: ", uzaklikAB)
        if (sayac == 2):
            uzaklikAC = uzaklik
            print("A ile C uzakligi: ", uzaklikAC)
        if (sayac == 3):
            uzaklikAD = uzaklik
            print("A ile D uzakligi: ", uzaklikAD)
        if (sayac == 4):
            uzaklikAE = uzaklik
            print("A ile E uzakligi: ", uzaklikAE)
        if (sayac == 5):
            uzaklikAF = uzaklik
            print("A ile F uzakligi: ", uzaklikAF)
        if (sayac == 6):
            uzaklikBC = uzaklik
            print("B ile C uzakligi: ", uzaklikBC)
        if (sayac == 7):
            uzaklikBD = uzaklik
            print("B ile D uzakligi: ", uzaklikBD)
        if (sayac == 8):
            uzaklikBE = uzaklik
            print("B ile E uzakligi: ", uzaklikBE)
        if (sayac == 9):
            uzaklikBF = uzaklik
            print("B ile F uzakligi: ", uzaklikBF)
        if (sayac == 10):
            uzaklikCD = uzaklik
            print("C ile D uzakligi: ", uzaklikCD)
        if (sayac == 11):
            uzaklikCE = uzaklik
            print("C ile E uzakligi: ", uzaklikCE)
        if (sayac == 12):
            uzaklikCF = uzaklik
            print("C ile F uzakligi: ", uzaklikCF)
        if (sayac == 13):
            uzaklikDE = uzaklik
            print("D ile E uzakligi: ", uzaklikDE)
        if (sayac == 14):
            uzaklikDF = uzaklik
            print("D ile F uzakligi: ", uzaklikDF)
        if (sayac == 15):
            uzaklikEF = uzaklik
            print("E ile F uzakligi: ", uzaklikEF)

        uzaklik = 0
        sayi = 0
        sayac += 1


# Hızlı Bir UPGMA Uygulaması (Aritmetik Ortalama ile Ağırlıksız Çift Grup Yöntemi)

# enKucuk_Hucre:
#   Tablodaki en küçük hücreyi bulur.
def enKucuk_Hucre(tablo):
    # Varsayılanı sonsuz olarak ayarla.
    min_hucre = float("inf")
    x, y = -1, -1

    # Her hücreye gidip, en küçük olanı ara.
    for i in range(len(tablo)):
        for j in range(len(tablo[i])):
            if tablo[i][j] < min_hucre:
                min_hucre = tablo[i][j]
                x, y = i, j

    # Hücrenin x ve y koordinatlarını döndür.
    return x, y


# etiketleri_Ekle:
#   Bir etiket listesinde iki etiketi birleştirir.
def etiketleri_Ekle(etiketler, a, b):
    # İndeksler sıralanmamışsa değiştirin.
    if b < a:
        a, b = b, a

    # İlk dizindeki etiketleri birleştirin.
    etiketler[a] = "(" + etiketler[a] + "," + etiketler[b] + ")"

    # İkinci dizindeki (artık gereksiz olan) etiketi kaldırın.
    del etiketler[b]


# tabloyu_Ekle:
#   Veri girişlerinin ortalamasını alarak hücredeki bir tablonun girdilerini (a, b) birleştirir.
def tabloyu_Ekle(tablo, a, b):
    # İndeksler sıralanmamışsa değiştirin.
    if b < a:
        a, b = b, a

    # Daha düşük indekslerde, i < A için; (A, i) olan tüm satırları yeniden oluşturun.
    satir = []
    for i in range(0, a):
        satir.append((tablo[a][i] + tablo[b][i]) / 2)
    tablo[a] = satir

    # Ardından, i > A için; tüm (i, A) sütunlarını yeniden oluşturun.
    #   Not: Matris daha küçük bir üçgen olduğundan dolayı, b satırı yalnızca b'den küçük indisleri için değerler taşır.
    for i in range(a + 1, b):
        tablo[i][a] = (tablo[i][a] + tablo[b][i]) / 2

    #   Değerlerin geri kalanını i satırından alıyoruz.
    for i in range(b + 1, len(tablo)):
        tablo[i][a] = (tablo[i][a] + tablo[i][b]) / 2
        # (Artık gereksiz olan) ikinci dizin sütunu girişini kaldırın.
        del tablo[i][b]

    # (Artık gereksiz olan) ikinci dizin satırını kaldırın.
    del tablo[b]


# UPGMA:
#   UPGMA algoritmasını etiketli bir tabloda çalıştırır.
def UPGMA(tablo, etiketler):
    # Tüm etiketler birleştirilene kadar işlemleri yap.
    while len(etiketler) > 1:
        # Tablodaki en küçük hücreyi bulun.
        x, y = enKucuk_Hucre(tablo)

        # Hücre koordinatlarındaki tabloya ekleme yapın.
        tabloyu_Ekle(tablo, x, y)

        # Etiketleri uygun şekilde güncelleyin.
        etiketleri_Ekle(etiketler, x, y)

    # Son etiketi döndürün.
    return etiketler[0]


## Burada, şu sitedeki UPGMA tablosu test edildi: https://www.youtube.com/watch?v=09eD4A_HxVQ&t=303s

# alfa_etiketleri:
#   Başlangıç harfinden bitiş harfine kadar etiketler yapar.
def alfa_etiketleri(baslangic, bitis):
    etiketler = []
    for i in range(ord(baslangic), ord(bitis) + 1):
        etiketler.append(chr(i))
    return etiketler


# Test tablosu verileri ve ilgili etiketler
M_etiketler = alfa_etiketleri("A", "F")  # A ile F
M = [
    [],  # A
    [uzaklikAB],  # B
    [uzaklikAC, uzaklikBC],  # C
    [uzaklikAD, uzaklikBD, uzaklikCD],  # D
    [uzaklikAE, uzaklikBE, uzaklikCE, uzaklikDE],  # E
    [uzaklikAF, uzaklikBF, uzaklikCF, uzaklikDF, uzaklikEF]  # F

]

print()
print("UPGMA: ", UPGMA(M, M_etiketler))
print()
agacVerisi = UPGMA(M, M_etiketler)
print(agacVerisi)
handle = StringIO(agacVerisi)
agac = Phylo.read(handle, "newick")

print(agac)

# UPGMA(M, M_etiketler) şu şekilde çıktı vermelidir: '(((A,B),E),(C,D))'


# ----------------------------------------------------------------------------------------


root.mainloop()


