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
plt.plot(l_values, time_taken_values, marker='o', linestyle='-')
plt.xlabel('lmax')
plt.ylabel('time taken (seconds)')

plt.subplot(1, 3, 2)
plt.plot(l_values, final_residual_values, marker='o', linestyle='-')
plt.xlabel('lmax')
plt.ylabel(r'$||r||_2$')

plt.subplot(1, 3, 3)
plt.plot(l_values, num_iter_values, marker='o', linestyle='-')
plt.xlabel('lmax')
plt.ylabel('number of coarse iterations')

plt.tight_layout()
plt.savefig("q3_part1.png", dpi=300)

"""----------------------------------------------------------------------------------------------"""

filename2 = '../q3_part2_lmax2.out'
filename8 = '../q3_part2_lmax8.out'

N_values_2              = []
time_taken_values_2     = []
final_residual_values_2 = []
num_iter_values_2       = []

N_values_8              = []
time_taken_values_8     = []
final_residual_values_8 = []
num_iter_values_8       = []


with open(filename2, 'r') as file2:
	content2 = file2.read()

with open(filename8, 'r') as file8:
	content8 = file8.read()


N_pattern_2              = re.compile(r'N = (\d+)')
time_pattern_2           = re.compile(r'time taken = ([\d\.]+)')
final_residual_pattern_2 = re.compile(r'final residual \|\|r\|\| = ([\d\.]+)')
num_iter_pattern_2       = re.compile(r'total no. of coarse level solves = ([\d\.]+)')  

N_pattern_8              = re.compile(r'N = (\d+)')
time_pattern_8           = re.compile(r'time taken = ([\d\.]+)')
final_residual_pattern_8 = re.compile(r'final residual \|\|r\|\| = ([\d\.]+)')
num_iter_pattern_8       = re.compile(r'total no. of coarse level solves = ([\d\.]+)')  


N_matches_2              = N_pattern_2.findall(content2)
time_matches_2           = time_pattern_2.findall(content2)
final_residual_matches_2 = final_residual_pattern_2.findall(content2)
num_iter_matches_2       = num_iter_pattern_2.findall(content2)

N_matches_8              = N_pattern_8.findall(content8)
time_matches_8           = time_pattern_8.findall(content8)
final_residual_matches_8 = final_residual_pattern_8.findall(content8)
num_iter_matches_8       = num_iter_pattern_8.findall(content8)


for i in range(len(N_matches_2)):
	N_values_2.append(int(N_matches_2[i]))  
	time_taken_values_2.append(float(time_matches_2[i]))  
	final_residual_values_2.append(float(final_residual_matches_2[i]))  
	num_iter_values_2.append(int(num_iter_matches_2[i]))

for i in range(len(N_matches_8)):
	N_values_8.append(int(N_matches_8[i]))  
	time_taken_values_8.append(float(time_matches_8[i]))  
	final_residual_values_8.append(float(final_residual_matches_8[i]))  
	num_iter_values_8.append(int(num_iter_matches_8[i]))

plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 1)
plt.plot(N_values_2, time_taken_values_2, marker='o', linestyle='-', label="lmax=2")
plt.plot(N_values_8, time_taken_values_8, marker='o', linestyle='-', label="lmax=8")
plt.xlabel('N')
plt.ylabel('time taken (seconds)')
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(N_values_2, final_residual_values_2, marker='o', linestyle='-', label="lmax=2")
plt.plot(N_values_8, final_residual_values_8, marker='o', linestyle='-', label="lmax=8")
plt.xlabel('N')
plt.ylabel(r'$||r||_2$')
plt.legend()

plt.subplot(1, 3, 3)
plt.plot(N_values_2, num_iter_values_2, marker='o', linestyle='-', label="lmax=2")
plt.plot(N_values_8, num_iter_values_8, marker='o', linestyle='-', label="lmax=8")
plt.xlabel('N')
plt.ylabel('number of coarse iterations')
plt.legend()


plt.tight_layout()
plt.savefig("q3_part2.png", dpi=300)

plt.show()

