# from __future__ import print_function
__author__ = 'sebasanper'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
from math import sqrt, log, pi
import time
import wake_geometry as wake
# 80 turbines

def analysis(v):
    output = open('matrix_larsen.dat', 'w')
    output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
    output2 = open('final_speed_larsen.dat', 'w')
    output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
    layout = open('../area_overlap/horns_rev.dat', 'r')  # Or horns_rev.dat for full layout
    windrose = open('horns_rev_windrose2.dat', 'r')
    draw = open('draw_horns_rev_larsen.dat', 'w')
    draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
    # draw2 = open('drawline.dat', 'w')
    turb_data = open('turb17_larsenviejo.dat', 'w')
    direction = open('direction_efficiency_larsen.dat', 'w')
    direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')
    nt = 80  # Or 80 for full layout
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
    summation = 0.0

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return 0.0003234808 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 **5.0 - 30.3162345004 * U0 **4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
        else:
            return 0.0

    r0 = 40.0  # Turbine rotor radius
    D = 2.0 * r0
    A = pi * r0 ** 2.0
    ct = 0.81
    deff = D * sqrt((1.0 + sqrt(1.0 - ct)) / (2.0 * sqrt(1.0 - ct)))
    H = 70.0  # Hub height
    ia = 0.08  # Ambient turbulence intensity according to vanlucanee. 8% on average
    rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
    r95 = 0.5 * (rnb + min(H, rnb))
    x0 = 9.5 * D / ((2.0 * r95 / deff) ** 3.0 - 1.0)
    c1 = (deff / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (ct * A * x0) ** (- 5.0 / 6.0)  # Prandtl mixing length

    for wind in range(0, len(windrose_speed)):
    # for wind in range(0, 1):
    #     U1 = windrose_speed[wind]  # Free stream wind speed
        # U0 = U1 * (70.0 / 10.0)**0.11 # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        U0 = 8.5
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        for turbine in range(nt):
            flag = [False for x in range(nt)]
            proportion = [0.0 for x in range(nt)]
            perpendicular_distance = [0.0 for x in range(nt)]
            parallel_distance = [0.0 for x in range(nt)]
            for i in range(nt):
                proportion[i], flag[i], perpendicular_distance[i], parallel_distance[i] = wake.determine_if_in_wake(layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i], A, c1, ct, angle3, r0)

            # Matrix with effect of each turbine <i = turbine> on every other turbine <j> of the farm
            for j in range(nt):
                if turbine != j and flag[j] is True:
                    # print (U0, ct, A, parallel_distance[j], perpendicular_distance[j], c1, x0, r95, rnb, deff)
                    wake_deficit_matrix[turbine][j] = proportion[j] * wake.wake_deficit(U0, ct, A, parallel_distance[j], perpendicular_distance[j], c1)
                elif turbine == j:
                    wake_deficit_matrix[turbine][j] = 0.0
                output.write('{0:f}\t'.format(wake_deficit_matrix[turbine][j]))
            output.write('\n')

        total_deficit = [0.0 for x in range(nt)]
        total_speed = [0.0 for x in range(nt)]
        for j in range(nt):
            for i in range(nt):
                total_deficit[j] += wake_deficit_matrix[i][j] ** 2.0
            total_deficit[j] = sqrt(total_deficit[j])
            total_speed[j] = U0 * (1.0 - total_deficit[j])
            output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(j, layout_x[j], layout_y[j], total_deficit[j], total_speed[j], int(angle), power(total_speed[j])))
        output2.write('\n')
        turb_data.write('{0:f}\t{1:f}\n'.format(angle, power(total_speed[14])))

        # Individual turbine efficiency
        # for g in range(nt):
        #      turb_data.write('{0:f}\n'.format(total_speed[g]))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        efficiency = 0.0
        for l in range(nt):
            profit += power(total_speed[l])
        efficiency = profit * 100.0 / (nt * power(U0))
        efficiency_proportion[wind] = efficiency / 100.0 * windrose_frequency[wind]
        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        direction.write('{0:f}\t{1:f}\n'.format(angle, efficiency))
        summation += efficiency_proportion[wind]
    print('total farm efficiency is {0:f} %'.format(summation))

    turb_data.close()
    output.close()
    output2.close()
    draw.close()
    direction.close()

if __name__ == '__main__':
    start_time = time.time()
    # for v in range(8, 50):
    #     try:
    #         # start_time2 = time.time()
    #         print float(v) / 2.0, analysis(float(v) / 2.0)
    #         # print("--- %s seconds ---" % (time.time() - start_time2))
    #     except ZeroDivisionError:
    #         continue
    print analysis(8.5)
    print("--- %s seconds ---" % (time.time() - start_time))