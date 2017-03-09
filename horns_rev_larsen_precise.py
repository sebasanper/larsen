__author__ = 'sebasanper'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
import wake_geometry as wake
from math import sqrt, log
import time

output = open('matrix_precise.dat', 'w')
output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
output2 = open('final_speed_precise.dat', 'w')
# output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
layout = open('horns_rev.dat', 'r')
windrose = open('horns_rev_windrose2.dat', 'r')
draw = open('draw_horns_rev_precise.dat', 'w')
draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
# draw2 = open('drawline.dat', 'w')
turb_data = open('turb_data_precise85ms.dat', 'w')
direction = open('direction_efficiency_precise.dat', 'w')
direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')
row = open('row_data_precise.dat', 'w')

def analysis():

    layout_x = []
    layout_y = []
    for line in layout:
        columns = line.split()
        layout_x.append(float(columns[0]))
        layout_y.append(float(columns[1]))

    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    layout.close()
    windrose.close()

    def Ct(U0):
        return 0.0001923077 * U0**4.0 + -0.0075407925 * U0**3.0 + 0.096462704 * U0**2.0 - 0.5012354312 * U0 + 1.7184749184

    # def power(U0):
    #     return -0.0071427272 * U0**5.0 + 0.5981106302 * U0**4.0 - 18.5218157059 * U0**3.0 + 251.0929636046 * U0**2.0 - 1257.8070377904 * U0 + 2043.2240149783

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 >= 4.0:
            return 19.7907842158 * U0 ** 2.0 - 74.9080669331 * U0 + 37.257967033  # From 4 to 11 m/s
    # def Cp(U0):
    #     return power(U0) / (0.5 * 1.225 * pi * r0**2.0 * U0**3.0)

    # for U0 in range(4, 20):
    nt = 80  # Number of turbines
    summation = 0.0
    p = [0.0 for x in range(nt)]
    ct = 0.0
    for wind in range(0, len(windrose_speed)):
        ct += 1
    # for wind in range(180-2, 180+3):
        # if wind in [100, 133, 271, 280, 313]:
        #     continue
        # U1 = windrose_speed[wind]  # Free stream wind speed
        # U0 = U1 * (70.0 / 10.0)**0.11 # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        U0 = 8.5
        k = 0.04  # Decay constant
        r0 = 40.0  # Turbine rotor radius
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        proportion = [[0.0 for x in range(nt)] for x in range(nt)]
        affected_matrix = [[] for x in range(nt)]

        for turbine in range(0, nt):
            flag = [False for x in range(nt)]
            for i in range(nt):
                if i != turbine:
                    proportion[i][turbine], flag[i] = wake.determine_if_in_wake(layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i], k, r0, angle3)
                elif i == turbine:
                    proportion[i][turbine] = 0.0
                    flag[i] = False
                if flag[i] is True:
                    affected_matrix[i].append(turbine)

        deficit = [0.0 for x in range(nt)]
        length = [0 for x in range(nt)]
        U = [0.0 for x in range(nt)]

        for i in range(nt):
            length[i] = len(affected_matrix[i])
            row.write('{0:f}\t{1:f}\t{2:d}\n'.format(layout_x[i], layout_y[i], len(affected_matrix[i])))
        row_vector = [[] for x in range(0, max(length) + 1)]

        for i in range(nt):
            for n in range(0, max(length) + 1):
                if length[i] == n:
                    row_vector[n].append(i)
        total_deficit = [0.0 for x in range(nt)]
        counter = 0
        for n in range(0, max(length) + 1):
        # for n in range(0, 5):
            if n == 0:
                for turb in row_vector[n]:
                    deficit[turb] = 0.0
                    U[turb] = U0
            else:
                for turb in row_vector[n]:
                    for i in affected_matrix[turb]:
                        # print U[i], turb, turb, i, i, counter
                        if U[i] == 0.0:
                            row_vector[n].append(turb)
                            counter += 1
                            continue
                        elif deficit_matrix[turb][i] == 0.0:
                            deficit_matrix[turb][i] = proportion[turb][i] * wake.wake_deficit(Ct(U[i]), k, wake.distance(layout_x[turb], layout_y[turb], layout_x[i], layout_y[i]), r0)
                            total_deficit[turb] += deficit_matrix[turb][i] ** 2.0
                    if counter >= 20000:
                        break
                    total_deficit[turb] = sqrt(total_deficit[turb])
                    U[turb] = U0 * (1.0 - total_deficit[turb])

        for turb in range(nt):
            for i in range(nt):
                output.write('{0:f}\t'.format(deficit_matrix[turb][i]))
            output.write('\n')
            output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(turb, layout_x[turb], layout_y[turb], total_deficit[turb], U[turb], int(angle), power(U[turb])))
        output2.write('\n')

        # Individual turbine efficiency
        # if angle == 0:
        # for g in range(nt):
        #     p[g] += power(U[g])/power(U[row_vector[0][0]])
        turb_data.write('{0:f}\t{1:f}\n'.format(angle, power(U[14])))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        for l in range(nt):
            profit += power(U[l])
        # print row_vector[0][0]
        efficiency = profit * 100.0 / (float(nt) * power(U[row_vector[0][0]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        print angle, efficiency
        direction.write('{0:f}\t{1:f}\n'.format(angle, efficiency))
        summation += efficiency_proportion[wind]
    print 'total farm efficiency is {0:f} %'.format(summation)

    # for g in range(nt):
    #     turb_data.write('{0:d}\t{1:f}\t{2:f}\t{3:f}\n'.format(g, layout_x[g], layout_y[g], p[g] / ct))
    # turb_data.write('\n')

    turb_data.close()
    output.close()
    output2.close()
    draw.close()
    direction.close()
    row.close()

if __name__ == '__main__':
    start_time = time.time()
    analysis()
print("--- %s seconds ---" % (time.time() - start_time))