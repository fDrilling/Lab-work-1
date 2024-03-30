import math

#Відкриття файла і зчитування данних. 
xyz_file = open("1_r9_xyz_dftV3#r9_001_AaabC.xyz")



# претворення в 2d list формату
# [[atom №, x, y, z]]
line_data_list = [] 
numberC = 0 
numberO = 0 
numberH = 0 
numberN = 0  
                                 
for d in xyz_file:
    sep_list = d.split()
    if sep_list[0] == 'C':
        line_data_list.append([f'{sep_list[0]}{numberC}', 
                            float(sep_list[1]), 
                            float(sep_list[2]), 
                            float(sep_list[3])])
        numberC += 1
    elif sep_list[0] == 'N':
        line_data_list.append([f'{sep_list[0]}{numberN}', 
                            float(sep_list[1]), 
                            float(sep_list[2]), 
                            float(sep_list[3])])
        numberN += 1
    elif sep_list[0] == 'H':
        line_data_list.append([f'{sep_list[0]}{numberH}', 
                            float(sep_list[1]), 
                            float(sep_list[2]), 
                            float(sep_list[3])])
        numberH += 1
    elif sep_list[0] == 'O':
        line_data_list.append([f'{sep_list[0]}{numberO}', 
                            float(sep_list[1]), 
                            float(sep_list[2]), 
                            float(sep_list[3])])
        numberO += 1

atom_number = numberC + numberO + numberH + numberN
xyz_file.close()
    

# Функція для обрахунку дистанції подавати дані атомів у вмгляді [atom, x, y, z] (line_data_list) видає результат у вигляді 10^-1 А(ангстреми)

def distance_calculation (atom1, atom2):
    sum_x = atom1[1] - atom2[1]
    sum_y = atom1[2] - atom2[2]
    sum_z = atom1[3] - atom2[3]
    result = math.hypot(sum_x, sum_y, sum_z)
    return result

# 2d матриця яка показує прари атомів, які зв'язані між собою bond_atoms_matrix

bond_atoms_matrix = []
for i in range (atom_number):
    for k in range (atom_number):
        distance = distance_calculation(line_data_list[i], line_data_list[k])
        if 4 < int(distance*10) < 17:
            bond_atoms_matrix.append([line_data_list[i][0],line_data_list[k][0]])


# CN = 13 - 15 
# CH = 6 - 8
# CO = 12 - 14
# CC = 14 - 16
# C2N = 11 - 13
# C2O = 13 - 15
            

# DA atoms number - 31
# DG atoms number - 32
# DT - 31
# DC - 29

if atom_number == 32:
    nucleoside_is = 'DG'
elif atom_number == 29:
    nucleoside_is = 'DC'
elif atom_number == 31:
    if numberN == 2:
        nucleoside_is = 'DT'
    else:
        nucleoside_is = 'DA'


# Шукаємо особливий N що звёязаний з 3 атомами С, остримуємо n_atom_3C = [C, C, C, N]

n_atom_3C = []
n = 0
for w in bond_atoms_matrix:
    if w[0][0] == 'N':
        n_atom_3C.append(w[1])
        for a in n_atom_3C:
            if a[0] == 'C':
                n += 1
        if n > 3:                     ####
            n_atom_3C.append(w[0])
            break
   
    
#Створюємо result list, де дані у вигляді [[atom_fileindex, [x,y,z], atom_num]], і додаємо туди N
line_data_dict = {}
for i in line_data_list:
    line_data_dict. update({i[0] :[i[1], i[2], i[3]]})

result_list = []
result_list.append([n_atom_3C[3], line_data_dict[n_atom_3C[3]]])
if nucleoside_is == 'DT' or 'DC':
    result_list[0].append('N6')
else:
    result_list[0].append('N9')




#Функція бере А№ i видає список атомів поряд
def near_atoms(atom = str):                
    sub_atoms = []
    for i in bond_atoms_matrix:
        if i[0] == atom:
            sub_atoms.append(i[1][0])
    return sub_atoms 

#Функція бере А№ i видає список атомів поряд з файловими індексами
def near_atomswindex(atom = str):                
    sub_atomsindex =[]
    for i in bond_atoms_matrix:
        if i[0] == atom:
            sub_atomsindex.append(i[1])
    return sub_atomsindex

# Функція бере А№ i A2№ знаходить спсисок атомів поряд з A№ і видалає з ноьго A2№ - уже опрацьований
def path_division(atom = str, already_checked = str):
    sub_atoms =[]
    for i in bond_atoms_matrix:
        if i[0] == atom:
            sub_atoms.append(i[1])

    sub_atoms.remove(already_checked)
    return sub_atoms

#функція отримує значення С№ і значення номеру атому в молекулі, і додає їх в результат
def update_result(atom_fileindex, atom_num): 

        result_list.append([atom_fileindex, line_data_dict[atom_fileindex]])
        result_list[-1].append(atom_num)
        





for i in n_atom_3C[0:3]:            # desoxiribose part
    ls = near_atoms(i) 

    if ls.count('O') == 1 and ls.count('C') == 1 and ls.count('H') == 1 and ls.count('N') == 1:  # По С1' > O
        update_result(i, 'C1`')
        divide = path_division(i, n_atom_3C[3])
        for s in divide:
            
            if s[0] == 'H':
                update_result(s, 'H1`')
            elif s[0] == 'O':
                new_atom = s
                update_result(s, 'O4`')
                
                for k in path_division(new_atom, i):                     # По О > C4`                                                                              
                    update_result(k, 'C4`')
                for w in path_division(k, s):                             #C4 >
                    if w[0] == 'H':
                        update_result(w, 'H4`')
                    elif w[0] == 'C':
                            
                            for l in path_division(w, k):           #C4`` > O5` or C2`
                                lst = near_atoms(l)
                                
                                if lst.count('C') == 2 and lst.count('H') == 2:    # C2`
                                    update_result(l, 'C2`')
                                    for h in near_atomswindex(l):
                                        if h[0] == 'H':
                                            update_result(h, 'H2`')
                                
                                elif  lst.count('C') == 1 and lst.count('H') == 1:         # O5` or O3`
                                    for g in near_atomswindex(l):
                                        if g[0] == 'C':
                                            jk = near_atoms(g)
                                            if jk.count('C') == 2:                 # C3` and O3` -> result
                                                update_result(g, 'C3`')
                                                update_result(l, 'O3`')
                                                O3_found = [l, 'O3']
                                                C3_found = [g, 'C3']
                                            elif jk.count('C') == 1:                  # C5` and O5` -> result
                                                update_result(g, 'C5`')
                                                update_result(l, 'O5`')
                                                O5_found = [l, 'O5']
                                                C5_found = [g, 'C5']


    if nucleoside_is == 'DT' or 'DC':
        ls = near_atoms(i)
        if ls.count('O') == 1  and ls.count('N') == 2:  # По С1' > O
            update_result(i, 'C2')
            for s in path_division(i, n_atom_3C[3]):
                if s[0] == 'O':
                    update_result(s, 'O2')
                
                elif s[0] == 'N':
                    update_result(s, 'N3')
                    for pl in path_division(s, i):
                        if pl[0] == 'H':
                            update_result(pl, 'H3')
                        
                        elif pl[0] == 'C':
                            update_result(pl, 'C4')
                           
                            for nk in path_division(pl, s):
                                if nk[0] == 'O':
                                    update_result(nk, 'O4')
                               
                                elif nk[0] == 'N':
                                    update_result(nk, 'N4')
                                    for nns in path_division(nk, pl):
                                        update_result (nns, 'H4')
                                
                                elif nk[0] == 'C':
                                    update_result(nk, 'C5')
                                   
                                    for c6 in path_division(nk, pl):
                                        if c6[0] == 'H':
                                            update_result(pl, 'H5')

                                        elif c6[0] == 'C':
                                            if near_atoms(c6).count('H') == 1:
                                                update_result(c6, 'C6')                                                
                                                for h6 in path_division(c6, nk):
                                                    if h6[0] == 'H':
                                                        update_result(h6, 'H6')

                                            elif near_atoms(c6).count('H') == 3:
                                                update_result(c6, 'C55')
                                               
                                                for hh in path_division(c6, nk):
                                                    update_result(hh, 'H55')


    if nucleoside_is == 'DG' or 'DA':
        ls = near_atoms(i)
        if ls.count('C') == 1  and ls.count('N') == 2:  # По С1' > O
            update_result(i, 'C4')
            for s in path_division(i, n_atom_3C[3]):
                if s[0] == 'N':
                    update_result(s, 'N3')
                    
                    for c2 in path_division(s, i):
                        update_result(c2, 'C2')
                        
                        for n1 in path_division(c2, s):
                            if n1[0] == 'H':
                                update_result(n1, 'H2')
                            
                            elif n1[0] == 'N':
                                if near_atoms(n1).count('H') == 2:
                                    update_result(n1, 'N2')                                    
                                    for h2 in path_division(n1, c2):
                                        if h2[0] == 'H':
                                            update_result(h2, 'H2')
                                
                                x = near_atoms(n1).count('H')                                                     
                                if x == 1 or x == 0:
                                    update_result(n1, 'N1')
                                    
                                    for c6 in path_division(n1, c2):
                                        update_result(c6, 'C6')
                                        
                                        for c5 in path_division(c6, n1):
                                            if c5[0] == 'O':
                                                update_result(c5, 'O6')
                                            
                                            elif c5[0] == 'N':
                                                update_result(c5, 'N6')
                                                for n6 in path_division(c5, c6):
                                                    update_result(n6, 'H6')
                                            
                                            elif c5[0] == 'C':
                                                update_result(c5, 'C5')
                                                
                                                for n7 in path_division(c5, c6):
                                                    if n7[0] == 'N':
                                                        update_result(n7, 'N7')
                                                        
                                                        for c8 in path_division(n7, c5):
                                                            update_result(c8, 'C8')
                                                            for h8 in path_division(c8, n7):
                                                                if h8[0] == 'H':
                                                                    update_result(h8, 'H8')



# Додаэмо в results залишкові атоми Н
def H_atoms_sugar(atom_found):
    H_dict = {'O3':'H3``', 'C3' : 'H3`', 'O5':'H5``', 'C5' : 'H5`'}
    for n in(near_atomswindex(atom_found[0])):
        if n[0] == 'H':
            var = H_dict[atom_found[1]]
            update_result(n, var)


atom_found_ls = [O3_found, C3_found, O5_found, C5_found]

if O3_found is None:
    pass
else:
    for h in atom_found_ls:
        H_atoms_sugar(h)



len(result_list)

with open('result.pdb', 'w') as file1:
    for i in range(atom_number):                                # !!!!!!!!range(atom_number-1) маэ бути

        rsp = 5 - len(result_list[i][2])
        
        x_val = str(result_list[i][1][0])[0:6]
        y_val = str(result_list[i][1][1])[0:6]
        z_val = str(result_list[i][1][2])[0:6]
        if len(x_val) < 6:
            x_val = f'{x_val}{'0'*(6 - len(x_val))}'
        if len(y_val) < 6:
            y_val = f'{y_val}{'0'*(6 - len(y_val))}'
        if len(z_val) < 6:
            z_val = f'{z_val}{'0'*(6 - len(z_val))}'

        atom = result_list[i][0][0]
        atom_ind = result_list[i][2]

        if  i+1 <= 9:
            space = ' '
        else:
            space = ''

        file1.write(f'ATOM    {i+1}{space}   {atom_ind}{' '*rsp}{nucleoside_is}   1      {x_val}   {y_val}   {z_val}   1.00 0.00           {atom}\n')                                                                


                                                    
                                         
print(line_data_list)
