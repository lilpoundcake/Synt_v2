#!/usr/bin/python3

# вытаскивание данных из файлов выхода DNAWorks

import primer3
import numpy as np
import pandas as pd

primer_list = []  # записывает все последовательности праймеров

sequence = ''  # общая последовательность синтетического гена

start_primer = 0  # счетчики
start_seq = 0

data_log = open("data_log.txt", "r")

# из логфайла берем название проекта для названий расчитаных олигов
file_name = data_log.readline().split()[1]

for n in range(int(data_log.readline().split()[1])):
    # n - из логфайла берем число фрагментов для прохода по всем файлам выхода
    primers = open(file_name + '_' + str(n + 1) + '.txt', 'r')
    # итеративно берем последовательности олигов и фрагментов гена из файлов выхода - primers

    for i in primers:
        # пробегаемся по строкам файла выхода и ищем только определенные строки. строка - i
        if 'oligonucleotides need to be synthesized' in i:
            start_primer += 1  # включаем счетчик для записи последовательностей олигов
            number_of_primers = int(i.strip().split(' ')[0])  # смотрим число олигов для синтеза фрагмента

        elif 'The DNA sequence' in i:
            start_seq += 1  # включаем счетчик для записи последовательности синтетического гена

        elif start_seq == 1:
            for z in primers:  # NB - а нужен ли тут вообще цикл???
                # выписываем последовательность синтетического гена
                if "---" in z:
                    start_seq = 0  # обнуляем счетчик. нуклеотидную последовательность гена выписали
                    break
                else:
                    if len(z.strip().split(' ')) > 1:
                        sequence += z.strip().split(' ')[1]
                    else:
                        continue

        elif start_primer == 1:
            for _ in range(number_of_primers):  # NB - а тут??? вроде и так у нас сверху цикл
                # выписываем олиги для синтеза
                primer_list.append(primers.readline().strip().split(' ')[1])
            start_primer = 0
            break
    primers.close()


# создаем файл, в который будут выписаны олиги, в формате SnapGene
list_primer = open(file_name + '_SG_primers.fasta', 'w')
# записываем олиги
for i in range(len(primer_list)):
    list_primer.write('>' + file_name + '_%i' % i + '\n' + primer_list[i] + '\n')
list_primer.close()

# записываем нуклеотидную последовательность гена
project_seq = open(file_name + '_sequence.fasta', 'w')
project_seq.write('>' + file_name + '\n' + sequence)
project_seq.close()

# создаем словарь название: последовательность
primers_files = open(file_name + '_SG_primers.fasta', 'r')
primers = dict()

for i in primers_files:
    primers[i.strip('>').strip('\n')] = primers_files.readline().strip()


# генерируем датафрейм для гомодимеров
df1 = pd.DataFrame(index=primers.keys(), columns=('Seq', 'Tm', 'Hairpin'))
df1.index.name = 'Name'
react_temp = int(data_log.readline().split()[1])  # просто для удобства дальнейших тестов


for i in df1.index:
    # считаем термодинамику олигов
    df1.loc[i, 'Seq'] = primers[i].upper()
    df1.loc[i, 'Tm'] = primer3.calcTm(primers[i], mv_conc=25.0, dv_conc=1.5,
                                      dntp_conc=0.8, dna_conc=0.04)
    df1.loc[i, 'Hairpin'] = primer3.calcHairpinTm(primers[i], mv_conc=25.0,
                                                  dv_conc=1.5, dntp_conc=0.8, dna_conc=0.04, temp_c=25)
    df1.loc[i, 'Homodimer'] = primer3.calcHomodimerTm(primers[i], mv_conc=25.0,
                                                      dv_conc=1.5, dntp_conc=0.8, dna_conc=0.04, temp_c=25)
    df1.loc[i, 'GC content'] = int(100 * (primers[i].count('G') + primers[i].count('C')) / len(primers[i]))
    df1.loc[i, 'React_temp'] = react_temp

# генерируем датафрейм для анализа гетеродимеров
df2 = pd.DataFrame(index=primers.keys(), columns=primers.keys())

# тестируем температуры отжига для гетеродимеров
# NB - может создать отдельные матрицы для анализа олигов только внутри одной реакции?
for i in df1.index:
    for j in df1.index:
        if df2.loc[i, j] == df2.loc[j, i] and df2.loc[i, j] != 'NaN':
            continue
        else:
            df2.loc[i, j] = primer3.calcHeterodimerTm(primers[i], primers[j], mv_conc=25.0, dv_conc=1.5,
                                                      dntp_conc=0.8, dna_conc=0.04, temp_c=25)
            df2.loc[j, i] = df2.loc[i, j]

# выбираем те олиги - которые имеют высокую температуру плавления гетеродимеров
hetero_Tm = dict()
n = 0
# ключ - порядковый номер. значение - олиг1, олиг2, Tm
for i in df2.index:
    for j in df2.index:
        if df2.loc[i, j] > react_temp:
            hetero_Tm[n] = [i, j, df2.loc[i, j]]
            n += 1
# форматируем как датафрейм
df3 = pd.DataFrame(hetero_Tm, index=['Seq1', 'Seq2', 'Tm'])  # для анализа гетеродимеров с высокой температурой отжига
df3 = df3.transpose()

# анализируем GC состав гена. Окно анализа - 30 нуклеотидов
GC_30 = dict()

for i in range(len(sequence) - 30):
    GC_30[i] = int(100 * (sequence[i: i + 30].count('G') + sequence[i: i + 30].count('C')) / len(sequence[i: i + 30]))

GC_content = pd.DataFrame(GC_30, index=['GC_content'])
GC_content = GC_content.transpose()

mean_GC = int(100 * (sequence.count('G') + sequence.count('C')) / len(sequence))

for i in GC_content.index:
    GC_content.loc[i, 'Max'] = 100
    GC_content.loc[i, 'Min'] = 0
    GC_content.loc[i, 'Mean'] = mean_GC

# вывод GC содержания гена
print("\nMean GC: " + str(mean_GC))

GC_content.to_csv(file_name + '_GC.csv')  # импорт данных по GC составу синтетического гена

# NB - надо добавить код по расчету уровня оптимизации гена!

# df1 - термодинамика олигов
# df2 - гомодимеры
# df3 - гетеродимеры

# выводим проблемные олиги, в которых высокотемпературные шпильки и высокотемпературные гетеродимеры
if len(df1[df1.Hairpin > react_temp].loc[:, ['Hairpin', 'Seq']]) != 0:
    print('\nHairpins \n', df1[df1.Hairpin > react_temp].loc[:, ['Hairpin', 'Seq']])

if len(df3) != 0 :
    print('\n\nHeterodimers\n', df3)

# и создаем отдельный файл с олигами у которых высокотемпературные шпильки
high_hairpin = df1[df1.Hairpin > react_temp].loc[:, ['Hairpin', 'Seq']]
high_hairpin.to_csv(file_name + "_hairpin.csv")

# файл с олигами, формирующими высокотемпературные гетеродимеры
df3.to_csv(file_name + "_high_temp.csv")

primers_files.close()
