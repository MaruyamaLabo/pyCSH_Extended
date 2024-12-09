import numpy as np
import os



def get_lammps_cshff(filename, entries_crystal, entries_bonds_cshff, entries_angle_cshff, supercell):
    """
    CSH_FF用のLAMMPS .dataファイルを生成する関数
    
    原子タイプの定義:
    1: Si
    2: Ca (非Interlayer)
    3: Cab (非Interlayer bridging site Ca ion)
    4: Cw (Interlayer Ca ion)
    5: O (非OH, 非水分子)
    6: Ob (bridging site O)
    7: Ooh (OH基のO)
    8: Oohw (Interlayer OH基のO)
    9: Ow (水分子のO)
    10: Hoh (非Interlayer, 非水分子)
    11: Hohw (Interlayer OH基のH)
    12: Hw (水分子のH)
    """
    atom_id_cshff = [sublist[-3] for sublist in entries_crystal]
    N_atom = max(atom_id_cshff)
    N_bond = len(entries_bonds_cshff)
    N_angle = len(entries_angle_cshff)

    with open(filename, 'w') as f:
        f.write( "Generated with Brickcode \n\n" )
        f.write( "{: 8d} atoms \n".format(N_atom) )
        f.write( "{: 8d} bonds \n".format(N_bond) )
        f.write( "{: 8d} angles \n".format(N_angle) )
        f.write( "{: 8d} atom types \n".format(12) )
        f.write( "{: 8d} bond types \n".format(2) )
        f.write( "{: 8d} angle types \n".format(1) )
        f.write( " \n" )
        f.write( "{: 12.6f} {: 12.6f} xlo xhi \n".format(0.0, supercell[0,0]) )
        f.write( "{: 12.6f} {: 12.6f} ylo yhi \n".format(0.0, supercell[1,1]) )
        f.write( "{: 12.6f} {: 12.6f} zlo zhi \n".format(0.0, supercell[2,2]) )
        f.write( "{: 12.6f} {: 12.6f} {: 12.6f} xy xz yz \n".format( supercell[1,0], supercell[2,0], supercell[2,1] ) )
        f.write( " \n" )
        
        # Masses セクション
        f.write("Masses\n\n")
        f.write("1 28.0855   # Si\n")
        f.write("2 40.078    # Ca\n")
        f.write("3 40.078    # Cab (Ca ion of Non-Interlayer bridging site)\n")
        f.write("4 40.078    # Cw (Ca ion of Interlayer)\n")
        f.write("5 15.9994   # O\n")
        f.write("6 15.9994   # Ob (bridging site O)")
        f.write("7 15.9994   # Ooh\n")
        f.write("8 15.9994   # Oohw (O of Interlayer OH)\n")
        f.write("9 15.9994   # Ow\n")
        f.write("10 1.00794  # Hoh\n")
        f.write("11 1.00794  # Hohw (H of Interlayer OH)\n")
        f.write("12 1.00794  # Hw\n\n")
        
        # Atoms セクション
        f.write("Atoms \n\n")
        
        fmt = "{: 8d} {: 8d} {: 8d} {: 8.3f} {: 12.6f} {: 12.6f} {: 12.6f}\n"
        
        for i, entry in enumerate(entries_crystal, 1):
            # 原子情報の取得
            atom_id = entry[7]
            atom_type_original = entry[1]  # 元々のタイプ情報
            x, y, z = entry[3:6]
            brick_id = entry[6]
            
            # 原子タイプの変換
            if atom_type_original == 2: #Si
                atom_type = 1
                charge = 1.72

            elif atom_type_original == 1: #Ca
                if brick_id in ['CII', 'CIU', 'CID', 'XU', 'XD']: #Cw (Interlayer Ca ion)
                    atom_type = 4
                    charge = 1.70
                elif brick_id in ['CU', 'CD']: #Cab (Non-Interlayer bridging site Ca ion)
                    atom_type = 3
                    charge = 1.70
                else: #Non-Interlayer Ca
                    atom_type = 2
                    charge = 1.43

            elif atom_type_original == 3: #'O' その他の酸素 ＊atom_type_original == 4は入れてはいけない

                if entry[6] in ['!L', '!R'] and entry[9] == "Up-bridging": #[!L,SU or SUo, !R]のとき
                    if entry[6] == '!L' and (entry[8] == 4 or entry[8] == 5): #O13と11 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    elif entry[6] == '!R' and entry[8] == 4: #O25 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    else:
                        atom_type = 5
                        charge = -1.14
                elif entry[6] in ['!L', '!Lo'] and entry[9] == "non-bridging": #[!L, !R] or [!Lo, CU, !Ro]のとき
                    if entry[6] == '!L' and entry[8] == 4: #O13 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    elif entry[6] == '!Lo' and entry[8] == 4: #O13 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    else:
                        atom_type = 5
                        charge = -1.14
                elif entry[6] in ['@L', '@R'] and entry[9] == "Down-bridging": #[!L,SU or SUo, !R]のとき
                    if entry[6] == '@R' and (entry[8] == 4 or entry[8] == 5): #O14と26 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    elif entry[6] == '@L' and entry[8] == 4: #O12 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    else:
                        atom_type = 5
                        charge = -1.14
                elif entry[6] in ['@R', '@Ro'] and entry[9] == "non-bridging": #[!L, !R] or [!Lo, CU, !Ro]のとき
                    if entry[6] == '@R' and entry[8] == 4: #O14 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    elif entry[6] == '@Ro' and entry[8] == 4: #O14 Ob (Si-bridging site O)
                        atom_type = 6
                        charge = -1.14
                    else:
                        atom_type = 5
                        charge = -1.14
                else:
                    atom_type = 5
                    charge = -1.14

            elif atom_type_original == 5: #'Ow' 水分子の酸素
                atom_type = 9
                charge = -0.82
            elif atom_type_original == 6:
                if brick_id in ['oDL', 'oDR', 'oUL', 'oUR', 'oXU', 'oXD', 'oMDL', 'oMDR', 'oMUL', 'oMUR']: #Oohw (Interlayer OH)
                    atom_type = 8
                    charge = 1.70
                else: #'Ooh' Non-INterlayer OH基の酸素
                    atom_type = 7
                    charge = -1.00

            elif atom_type_original == 7: #'Hw' 水分子の水素
                atom_type = 12
                charge = 0.41
            elif atom_type_original == 8:
                if brick_id in ["oDL","oDR","oUL","oUR","oXU","oXD","oMDL","oMDR","oMUL","oMUR"]: #Hohw (Interlayer OH)
                    atom_type = 11
                    charge = 1.70
                else: #'Hoh' その他の水素
                    atom_type = 10
                    charge = 0.29
                    
            # すべての原子に対してmol_ID = 1を使用, 元Oshellは除外
            if atom_type_original != 4:
                f.write(fmt.format(atom_id, 1, atom_type, charge, x, y, z))
        
        # Bonds セクション
        if entries_bonds_cshff:
            f.write("\nBonds\n\n")
            for i, entry in enumerate(entries_bonds_cshff, 1):
                if entry[1] != 3: # 3 is O(S)-O
                     f.write(f"{i} {entry[1]} {entry[2]} {entry[3]}\n")
        
        # Angles セクション
        if entries_angle_cshff:
            f.write("\nAngles\n\n")
            for i, entry in enumerate(entries_angle_cshff, 1):
                f.write(f"{i} {entry[1]} {entry[2]} {entry[3]} {entry[4]}\n")

def get_lammps_input_Erica_umeki(input_file, entries_crystal, entries_bonds, entries_angle, supercell, write_lammps_erica):
    
    N_atom = len(entries_crystal)
    N_bond = len(entries_bonds)
    N_angle = len(entries_angle)

    with open(input_file, "w") as f:
        f.write( "Generated with Brickcode \n\n" )
        f.write( "{: 8d} atoms \n".format(N_atom) )
        f.write( "{: 8d} bonds \n".format(N_bond) )
        f.write( "{: 8d} angles \n".format(N_angle) )
        f.write( "{: 8d} atom types \n".format(12) )
        f.write( "{: 8d} bond types \n".format(3) )
        f.write( "{: 8d} angle types \n".format(3) )
        f.write( " \n" )
        f.write( "{: 12.6f} {: 12.6f} xlo xhi \n".format(0.0, supercell[0,0]) )
        f.write( "{: 12.6f} {: 12.6f} ylo yhi \n".format(0.0, supercell[1,1]) )
        f.write( "{: 12.6f} {: 12.6f} zlo zhi \n".format(0.0, supercell[2,2]) )
        f.write( "{: 12.6f} {: 12.6f} {: 12.6f} xy xz yz \n".format( supercell[1,0], supercell[2,0], supercell[2,1] ) )
        f.write( " \n" )
        f.write( "Masses \n" )
        f.write( " \n" )
        f.write( "1 40.08  #Ca  \n" )
        f.write( "2 28.10  #Si \n" )
        f.write( "3 15.59  #O \n" )
        f.write( "4 0.40   #O(S) \n" )
        f.write( "5 16.00  #Ow \n" )
        f.write( "6 16.00  #Oh \n" )
        f.write( "7 1.00   #Hw \n" )
        f.write( "8 1.00   #H \n" )
        f.write( "9 40.08  #Cw \n" )
        f.write( "10 40.08  #Cab \n" )
        f.write( "11 16.00  #Ohw \n" )
        f.write( "12 1.00  #How \n" )
        f.write( " \n" )
        f.write( "Atoms \n" )
        f.write( " \n" )
        fmt = "{: 8d} {: 8d} {: 8d} {: 8.3f} {: 12.6f} {: 12.6f} {: 12.6f}\n"
        molID = 2
        CSID = 2
        CS_info = []
        H2O_H1 = True
        for i in entries_crystal:
            # 3,4(O,O(S))は連続する
            if i[1] == 3: #O
                f.write( fmt.format(i[0], molID, *i[1:6]) )
                CS_info.append( [i[0], CSID] )
            elif i[1] == 4: #O(S)
                f.write( fmt.format(i[0], molID, *i[1:6]) )
                CS_info.append( [i[0], CSID] )
                molID += 1
                CSID += 1
            elif i[1] == 1: #Ca
                if i[6] in ['CII', 'CIU', 'CID', 'XU', 'XD']: #Cw (Interlayer Ca ion)
                    f.write( fmt.format(i[0], 1, 9, *i[2:6]) )
                    CS_info.append( [i[0], 1] )
                elif i[6] in ['CU', 'CD']: #Cab (Non-Interlayer bridging site Ca ion)
                    f.write( fmt.format(i[0], 1, 10, *i[2:6]) )
                    CS_info.append( [i[0], 1] )
                else:
                    f.write( fmt.format(i[0], 1, *i[1:6]) )
                    CS_info.append( [i[0], 1] )
            # 5,7, 7(Ow,Hw,Hw)は連続する
            elif i[1] == 5: #Ow
                f.write( fmt.format(i[0], molID, *i[1:6]) )
                CS_info.append( [i[0], 1] )
            elif i[1] == 7: #Hw
                if H2O_H1 == True:
                    H2O_H1 = False
                    f.write( fmt.format(i[0], molID, *i[1:6]) )
                    CS_info.append( [i[0], 1] )
                else: #H2O_H1 == False
                    H2O_H1 = True
                    f.write( fmt.format(i[0], molID, *i[1:6]) )
                    CS_info.append( [i[0], 1] )
                    molID += 1
            # 6,8(Oh,H)は連続する, 11,12(Ohw,How)は連続する
            elif i[1] == 6: #Oh
                if i[6] in ["oDL","oDR","oUL","oUR","oXU","oXD","oMDL","oMDR","oMUL","oMUR"]: #Ohw (Interlayer O of OH)
                    f.write( fmt.format(i[0], molID, 11, *i[2:6]) )
                    CS_info.append( [i[0], 1] )
                else: #oh 
                    f.write( fmt.format(i[0], molID, *i[1:6]) )
                    CS_info.append( [i[0], 1] )
            elif i[1] == 8: #H
                if i[6] in ["oDL","oDR","oUL","oUR","oXU","oXD","oMDL","oMDR","oMUL","oMUR"]: #How (Interlayer H of OH)
                    f.write( fmt.format(i[0], molID, 12, *i[2:6]) )
                    CS_info.append( [i[0], 1] )
                    molID += 1
                else: #H
                    f.write( fmt.format(i[0], molID, *i[1:6]) )
                    CS_info.append( [i[0], 1] )
                    molID += 1
            else:
                f.write( fmt.format(i[0], 1, *i[1:6]) )
                CS_info.append( [i[0], 1] )
        f.write( " \n" )

        f.write( "Bonds \n" )
        f.write( " \n" )
        fmt = "{: 8d} {: 8d} {: 8d} {: 8d} \n"
        for i in entries_bonds:
            f.write( fmt.format(*i) )
        f.write( " \n" )

        f.write( "Angles \n" )
        f.write( " \n" )
        fmt = "{: 8d} {: 8d} {: 8d} {: 8d} {: 8d} \n"
        for i in entries_angle:
            f.write( fmt.format(*i) )
        f.write( " \n" )

        f.write( " \n" )
        fmt = "{: 8d} {: 8d} \n"
        if write_lammps_erica:
            f.write( "CS-Info \n" )
            f.write( "\n" )
            for i in CS_info:
                f.write( fmt.format(*i))