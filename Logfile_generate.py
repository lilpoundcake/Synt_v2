#!/usr/bin/python3

import math  # для округления длин фрагментов

log_title = ""  # для шапки файлов, подаваемых на вход DNAWorks

project_name = input("Введите название проекта\n")

react_temp = input("Введите температуру проведения реакции\n")
log_title += ("melting low " + str(react_temp) + "\n")

primer_lenght = input("Введите желаемую длину праймеров\n")
log_title += ("length low " + str(primer_lenght) + "\n")

log_title += "tbio\n"  # параметр для расчета олигов от центра к краям

sys_expr = int(input("Укажите систему экспрессии:\n1 - E.coli\n2 - CHO\n3 - Human\n4 - Sf9\n"))

if sys_expr == 1:
    log_title += "codon E.coli\n"
elif sys_expr == 3:
    log_title += "codon H. sapiens\n"
elif sys_expr == 2:
    CHO = open("CHO.txt", "r")
    for i in CHO:
        log_title += (i + "\n")
    CHO.close()
elif sys_expr == 4:
    Sf9 = open("Sf9.txt", "r")
    for i in Sf9:
        log_title += (i + "\n")
    Sf9.close()

if int(input("Последовательность аминокислотная или нуклеотидная?\n1 - аминокислотная\n2 - нуклеотидная\n")) == 1:
    seq_prot = input("Введите последовательность белка\n")
    seq_nt = 0
else:
    seq_nt = input("Введите нуклеотидную последовательность\n")
    seq_prot = 0

if seq_nt == 0:
    codon_freq = input("Введите пороговое значение частоты встречаемости кодонов\n")
    log_title += ("frequency threshold " + codon_freq + "\n")

fragment_list = []  # для хранения отдельных фрагментов гена. 1 фрагмент - 1 реакция синтеза

if seq_prot == 0:
    seq_len = len(seq_nt)
else:
    seq_len = len(seq_prot)

if seq_len >= 500 and seq_prot == 0:  # разбивание последовательности на фрагменты длиной не более 500 нуклеотидов
    fragment_num = math.ceil(len(seq_nt) / 500)
    fragment_len = math.ceil(seq_len / fragment_num)

    for i in range(fragment_num):
        if i == int(fragment_num - 1):
            fragment_list.append(seq_nt[(fragment_len * i):])
        else:
            fragment_list.append(seq_nt[(fragment_len * i):(fragment_len * (i + 1))])

elif seq_len < 500 and seq_prot == 0:
    fragment_num = 1
    fragment_list.append(seq_nt)

elif seq_len >= 166 and seq_nt == 0:  # разбивание последовательности на фрагменты не более 166 а/к
    fragment_num = int(math.ceil(seq_len / 166))
    fragment_len = math.ceil(seq_len / fragment_num)

    for i in range(fragment_num):
        if i == int(fragment_num - 1):
            fragment_list.append(seq_prot[(fragment_len * i):])
        else:
            fragment_list.append(seq_prot[fragment_len * i:(fragment_len * (i + 1))])

elif seq_len < 166 and seq_nt == 0:
    fragment_num = 1
    fragment_list.append(seq_prot)

# запись всех input файлов для DNAWorks

for i in range(len(fragment_list)):
    output_i = open(project_name + "_" + str(i + 1) + ".inp", "w")

    if seq_prot == 0:
        output_i.write(log_title + "\nlogfile " + project_name + "_" + str(i + 1) + ".txt\n" + "\nnucleotide\n")

        for j in range(int(math.ceil((len(fragment_list[i])) / 100))):
            if j == (int(len(fragment_list[i]) / 100)):
                output_i.write(fragment_list[i][j * 100:] + "\n")
            else:
                output_i.write(fragment_list[i][j * 100:(j + 1) * 100] + "\n")

        output_i.write("//\n" + "\n#" + str(len(fragment_list[i])))

    else:
        output_i.write(log_title + "\nlogfile " + project_name + "_" + str(i + 1) + ".txt\n" + "\nprotein\n" +
                       fragment_list[i] + "\n//\n" + "\n#" + str(len(fragment_list[i])))
    output_i.close()

# в этом файле хранятся параметры проведения реакции, название проекта, число и последовательности фрагментов

data_log = open("data_log.txt", "w")
data_log.write("project_name " + project_name + "\nfragment_num " + str(fragment_num) +
               "\nreact_temp " + str(react_temp) + "\n\n")

for i in fragment_list:
    data_log.write("\n# " + i)
data_log.close()

# пишем скрипт для исполнения расчетов в DNAWorks. Подаем на вход сгенерированные файлы с последовательностями
# и условиями проведения реакции

bash_script = open("bash_next.sh", "w")
bash_script.write("#!/bin/bash\n")

for i in range(fragment_num):
    # запуск программы и подача в качестве аргументов файлы с фрагментами и параментрами реакции
    bash_script.write("./dnaworks " + project_name + "_" + str(i + 1) + ".inp\n")

bash_script.write("python3 Primer_Alg.py\n")  # тут произойдет запуск скрипта для анализа расчитанных олигов

# сохранение всех выводов в отдельную папку
bash_script.write("mkdir " + project_name + "\n")

# для переноса файлов вввода и вывода
for i in range(fragment_num):
    bash_script.write("cp " + project_name + "_" + str(i + 1) + ".inp " + project_name + "\n" +
                      "cp " + project_name + "_" + str(i + 1) + ".txt " + project_name + "\n")

# для перемещения файлов с расчетами термодинамики олигов
bash_script.write("cp " + project_name + "_hairpin.csv " + project_name + "\n" +
                  "cp " + project_name + "_high_temp.csv " + project_name + "\n" +
                  "cp " + project_name + "_SG_primers.fasta " + project_name + "\n" +
                  "cp " + project_name + "_sequence.fasta " + project_name + "\n")

bash_script.write("cp -r " + project_name + " ..\n")

# удаляю программу DNAWorks
bash_script.write("cd ..\n" +
    "rm -rf ./DNAWorks ./For_DNAWorks\n")

# упаковка всех выводов в архив
bash_script.write("tar -zcf " + project_name + ".tar.gz " + project_name +
                  "/\nrm -rf ./" + project_name)

bash_script.close()
