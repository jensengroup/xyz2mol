from collections import defaultdict

__ATOM_LIST__ = \
    ['h',  'he',
     'li', 'be', 'b',  'c',  'n',  'o',  'f',  'ne',
     'na', 'mg', 'al', 'si', 'p',  's',  'cl', 'ar',
     'k',  'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu',
     'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',
     'rb', 'sr', 'y',  'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag',
     'cd', 'in', 'sn', 'sb', 'te', 'i',  'xe',
     'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy',
     'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w',  're', 'os', 'ir', 'pt',
     'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn',
     'fr', 'ra', 'ac', 'th', 'pa', 'u',  'np', 'pu']


atomic_valence = defaultdict(list)
atomic_valence[1] = [1]
atomic_valence[5] = [3,4]
atomic_valence[6] = [4]
atomic_valence[7] = [3,4]
atomic_valence[8] = [2,1,3]
atomic_valence[9] = [1]
atomic_valence[11] = [1] #not really needed
atomic_valence[14] = [4]
atomic_valence[15] = [5,3] #[5,4,3]
atomic_valence[16] = [6,3,2] #[6,4,2]
atomic_valence[17] = [1]
atomic_valence[21] = [4] #not really needed
atomic_valence[22] = [4] #not really needed
atomic_valence[23] = [4] #not really needed
atomic_valence[24] = [4] #not really needed
atomic_valence[25] = [4] #not really needed
atomic_valence[26] = [8] #not really needed
atomic_valence[27] = [4] #not really needed
atomic_valence[28] = [4] #not really needed
atomic_valence[29] = [4] #not really needed
atomic_valence[30] = [4] #not really needed
atomic_valence[32] = [4]
atomic_valence[35] = [1]
atomic_valence[39] = [4] #not really needed
atomic_valence[40] = [4] #not really needed
atomic_valence[41] = [4] #not really needed
atomic_valence[42] = [4] #not really needed
atomic_valence[43] = [4] #not really needed
atomic_valence[44] = [4] #not really needed
atomic_valence[45] = [4] #not really needed
atomic_valence[46] = [4] #not really needed
atomic_valence[47] = [4] #not really needed
atomic_valence[48] = [4] #not really needed
atomic_valence[53] = [1]
atomic_valence[57] = [4] #not really needed
atomic_valence[72] = [4] #not really needed
atomic_valence[73] = [4] #not really needed
atomic_valence[74] = [4] #not really needed
atomic_valence[75] = [4] #not really needed
atomic_valence[76] = [4] #not really needed
atomic_valence[77] = [4] #not really needed
atomic_valence[78] = [4] #not really needed
atomic_valence[79] = [4] #not really needed
atomic_valence[80] = [4] #not really needed

TM_charge_list = {}
TM_charge_list[11] = [1]
TM_charge_list[21] = [3]
TM_charge_list[22] = [4, 2, 3]
TM_charge_list[23] = [5, 4, 3, 2]
TM_charge_list[24] = [3, 6, 2, 5]
TM_charge_list[25] = [2, 4, 7, 3, 5, 6]
TM_charge_list[26] = [2, 3, 4, 6]
TM_charge_list[27] = [2, 3]
TM_charge_list[28] = [2, 3, 4]
TM_charge_list[29] = [2, 1, 3]
TM_charge_list[30] = [2]

TM_charge_list[39] = [3]
TM_charge_list[40] = [4, 3, 2]
TM_charge_list[41] = [5, 4, 3, 2]
TM_charge_list[42] = [3, 6, 2, 4, 5]
TM_charge_list[43] = [4, 7, 6, 5]
TM_charge_list[44] = [2, 3, 4, 6, 8]
TM_charge_list[45] = [2, 3, 4]
TM_charge_list[46] = [2, 4, 0]
TM_charge_list[47] = [1, 2]
TM_charge_list[48] = [2]

TM_charge_list[57] = [3]
TM_charge_list[72] = [4, 3, 2]
TM_charge_list[73] = [5, 4, 3, 2]
TM_charge_list[74] = [3, 6, 2, 4, 5]
TM_charge_list[75] = [1, 4, 7, 6, 5, 2, 3]
TM_charge_list[76] = [4, 6, 8, 2, 3]
TM_charge_list[77] = [3, 4, 1, 2]
TM_charge_list[78] = [2, 4, 0]
TM_charge_list[79] = [1, 3, 2]
TM_charge_list[80] = [2, 1]

atomic_valence_electrons = {}
atomic_valence_electrons[1] = 1
atomic_valence_electrons[5] = 3
atomic_valence_electrons[6] = 4
atomic_valence_electrons[7] = 5
atomic_valence_electrons[8] = 6
atomic_valence_electrons[9] = 7
atomic_valence_electrons[11] = 1
atomic_valence_electrons[11] = 0
atomic_valence_electrons[14] = 4
atomic_valence_electrons[15] = 5
atomic_valence_electrons[16] = 6
atomic_valence_electrons[17] = 7
atomic_valence_electrons[21] = 3
atomic_valence_electrons[22] = 4
atomic_valence_electrons[23] = 5
atomic_valence_electrons[24] = 6
atomic_valence_electrons[25] = 7
atomic_valence_electrons[26] = 8
atomic_valence_electrons[27] = 9
atomic_valence_electrons[28] = 10
atomic_valence_electrons[29] = 11
atomic_valence_electrons[30] = 2
atomic_valence_electrons[32] = 4
atomic_valence_electrons[35] = 7
atomic_valence_electrons[39] = 3
atomic_valence_electrons[40] = 4
atomic_valence_electrons[41] = 5
atomic_valence_electrons[42] = 6
atomic_valence_electrons[43] = 7
atomic_valence_electrons[44] = 8
atomic_valence_electrons[45] = 9
atomic_valence_electrons[46] = 10
atomic_valence_electrons[47] = 11
atomic_valence_electrons[48] = 2
atomic_valence_electrons[53] = 7
atomic_valence_electrons[57] = 3
atomic_valence_electrons[72] = 4
atomic_valence_electrons[73] = 5
atomic_valence_electrons[74] = 6
atomic_valence_electrons[75] = 7
atomic_valence_electrons[76] = 8
atomic_valence_electrons[77] = 9
atomic_valence_electrons[78] = 10
atomic_valence_electrons[79] = 11
atomic_valence_electrons[80] = 2

atom_aos = 2*[1] +\
           8*[4] +\
           8*[4] +\
           2*[4] + 9*[9] + 7*[4] + \
           2*[4] + 9*[9] + 7*[4] + \
           2*[4] + 1*[9] + 14*[0] + 8*[9] + 7*[4]   

TMs = {11,
       21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
       39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
       57, 72, 73, 74, 75, 76, 77, 78, 79, 80}