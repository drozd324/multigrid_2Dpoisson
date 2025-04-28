import re
import matplotlib.pyplot as plt

"""----------------------------------------------------------------------------------------------"""

filename = '../q3_part1.out'

l_values = []
time_taken_values = []
final_residual_values = []
num_iter_values = []

with open(filename, 'r') as file:
	content = file.read()

l_pattern = re.compile(r'l = (\d+)')  
time_pattern = re.compile(r'time taken = ([\d\.]+)')  
final_residual_pattern = re.compile(r'final residual \|\|r\|\| = ([\d\.]+)')  
num_iter_pattern = re.compile(r'total no. of coarse level solves = ([\d\.]+)')  

l_matches = l_pattern.findall(content)
time_matches = time_pattern.findall(content)
final_residual_matches = final_residual_pattern.findall(content)
num_iter_matches = num_iter_pattern.findall(content)

for i in range(len(l_matches)):
	l_values.append(int(l_matches[i]))  
	time_taken_values.append(float(time_matches[i]))  
	final_residual_values.append(float(final_residual_matches[i]))  
	num_iter_values.append(int(num_iter_matches[i]))

plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 1)
plt.plot(l_values, time_taken_values, marker='o', linestyle='-', color='b')
plt.xlabel('lmax')
plt.ylabel('time taken')

plt.subplot(1, 3, 2)
plt.plot(l_values, final_residual_values, marker='o', linestyle='-', color='r')
plt.xlabel('lmax')
plt.ylabel(r'$||r||_2$')

plt.subplot(1, 3, 3)
plt.plot(l_values, num_iter_values, marker='o', linestyle='-', color='g')
plt.xlabel('lmax')
plt.ylabel('number of coarse iterations')

plt.tight_layout()
plt.savefig("q3_part1.png", dpi=300)

"""----------------------------------------------------------------------------------------------"""

filename = '../q3_part2.out'

N_values = []
time_taken_values = []
final_residual_values = []
num_iter_values = []

with open(filename, 'r') as file:
	content = file.read()

N_pattern = re.compile(r'N = (\d+)')
time_pattern = re.compile(r'time taken = ([\d\.]+)')
final_residual_pattern = re.compile(r'final residual \|\|r\|\| = ([\d\.]+)')
num_iter_pattern = re.compile(r'total no. of coarse level solves = ([\d\.]+)')  

N_matches = N_pattern.findall(content)
time_matches = time_pattern.findall(content)
final_residual_matches = final_residual_pattern.findall(content)
num_iter_matches = num_iter_pattern.findall(content)

for i in range(len(N_matches)):
	N_values.append(int(N_matches[i]))  
	time_taken_values.append(float(time_matches[i]))  
	final_residual_values.append(float(final_residual_matches[i]))  
	num_iter_values.append(int(num_iter_matches[i]))

plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 1)
plt.plot(N_values, time_taken_values, marker='o', linestyle='-', color='b')
plt.xlabel('N')
plt.ylabel('Time Taken')

plt.subplot(1, 3, 2)
plt.plot(N_values, final_residual_values, marker='o', linestyle='-', color='r')
plt.xlabel('N')
plt.ylabel(r'$||r||_2$')

plt.subplot(1, 3, 3)
plt.plot(N_values, num_iter_values, marker='o', linestyle='-', color='g')
plt.xlabel('N')
plt.ylabel('number of coarse iterations')

plt.tight_layout()
plt.savefig("q3_part2.png", dpi=300)

plt.show()

