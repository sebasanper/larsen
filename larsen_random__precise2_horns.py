__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Larsen wake model with partial shadowing factor applied to horns rev.
import wake_geometry as wake
from math import sqrt, log, tan, cos, pi, floor, ceil
from numpy import deg2rad
import time
from os import makedirs, path, chdir
from random import random
powers = ['power7', 'power5', 'power3', 'powertable']
# powers = ['power7']
thrusts = ['Ct6', 'Ct4', 'Ct3', 'Cttable', 'Ctstep']
# thrusts = ['Ct6']

for powertype in powers:  # Loop over power curves
    for thrusttype in thrusts:  # Loop over thrust curves
        print powertype, thrusttype
        for withdata in [False]:  # Either measure execution time, or run once completely.
            for rose in ['360']:#, '30']:  # Either run with 30 degrees, or 360 degrees wind roses.

                newpath = path.join('larsen_results/', powertype, thrusttype, rose)

                def power(U0):
                    #power7, power5, power3, powertable, powerstep
                    if powertype == 'power7':
                        return power7(U0)
                    elif powertype == 'power5':
                        return power5(U0)
                    elif powertype == 'power3':
                        return power3(U0)
                    elif powertype == 'powertable':
                        return powertable(U0)
                    elif powertype == 'powerstep':
                        return powerstep(U0)

                def Ct(U0):
                    # Ct6, Ct4, Ct3, Cttable, Ctstep
                    if thrusttype == 'Ct6':
                        return Ct6(U0)
                    elif thrusttype == 'Ct4':
                        return Ct4(U0)
                    elif thrusttype == 'Ct3':
                        return Ct3(U0)
                    elif thrusttype == 'Cttable':
                        return Cttable(U0)
                    elif thrusttype == 'Ctstep':
                        return Ctstep(U0)

                chdir('/home/sebasanper/PycharmProjects/larsenEWTS2')
                layout = open('horns_rev.dat', 'r')
                if rose == '360':
                    windrose = open('horns_rev_windrose2.dat', 'r')
                elif rose == '30':
                    windrose = open('horns_rev_windrose.dat', 'r')

                if not path.exists(newpath):
                    makedirs(newpath)
                chdir(newpath)
                if withdata:
                    turb_data = open('random_turb14_power.dat', 'w', 1)
                    direction = open('random_direction_power.dat', 'w', 1)
                    eff = open('random_efficiency.dat', 'w', 1)
                if not withdata:
                    exe_time = open('random_exe_time.dat', 'w', 1)

                def Ct_table(U0):
                    v = U0
                    if v < 4: return 0.1
                    if v == 4: return 0.82
                    if v == 5: return 0.81
                    if v == 6: return 0.8
                    if v == 7: return 0.81
                    if v == 8: return 0.81
                    if v == 9: return 0.78
                    if v == 10: return 0.74
                    if v == 11: return 0.65
                    if v == 12: return 0.57
                    if v == 13: return 0.41
                    if v == 14: return 0.31
                    if v == 15: return 0.25
                    if v == 16: return 0.2
                    if v == 17: return 0.17
                    if v == 18: return 0.14
                    if v == 19: return 0.12
                    if v == 20: return 0.1
                    if v == 21: return 0.09
                    if v == 22: return 0.08
                    if v == 23: return 0.07
                    if v == 24: return 0.06
                    if v == 25: return 0.05

                def Cttable(U0):
                    return interpolate(floor(U0), Ct_table(floor(U0)), ceil(U0), Ct_table(ceil(U0)), U0)

                def Ctstep(U0):
                    if U0 < 4.0:
                        return 0.1
                    elif U0 <= 9.5:
                        return 0.80
                    elif U0 <= 25:
                        return 0.250625
                    else:
                        return 0.0

                def Ct6(U0):
                    if U0 < 4.0:
                        return 0.1
                    elif U0 <= 25.0:
                        return 7.3139922126945e-7 * U0 ** 6.0 - 6.68905596915255e-5 * U0 ** 5.0 + 2.3937885e-3 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
                    else:
                        return 0.0

                def Ct4(U0):
                    if U0 < 4.0:
                        return 0.1
                    elif U0 <= 25.0:
                        return - 3.10224672352816e-5 * U0 ** 4.0 + 0.0021367624 * U0 ** 3.0 - 0.0495873986 * U0 ** 2.0 + 0.3976324804 * U0 - 0.1608576035
                    else:
                        return 0.0

                def Ct3(U0):
                    if U0 < 4.0:
                        return 0.1
                    elif U0 <= 25.0:
                        return 3.374593e-4 * U0 ** 3.0 - 0.0136412226 * U0 ** 2.0 + 0.1118003309 * U0 + 0.5782039288
                    else:
                        return 0.0

                def power7(U0):
                    if U0 < 4.0:
                        return 0.0
                    elif U0 <= 25.0:
                        return 3.234808e-4 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 ** 5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
                    else:
                        return 0.0

                def power5(U0):
                    if U0 < 4.0:
                        return 0.0
                    elif U0 <= 25.0:
                        return - 0.0110778061 * U0 ** 5.0 + 0.8986075613 * U0 ** 4.0 - 27.2165513154 * U0 ** 3.0 + 368.8877606215 * U0 ** 2.0 - 1994.1905079276 * U0 + 3712.3986113386
                    else:
                        return 0.0

                def power3(U0):
                    if U0 < 4.0:
                        return 0.0
                    elif U0 <= 25.0:
                        return - 0.5308414162 * U0 ** 3.0 + 15.4948143381 * U0 ** 2.0 + 13.1508234816 * U0
                    else:
                        return 0.0

                def power_table(U0):
                    v = U0
                    if v < 4: return 0.0
                    if v == 4: return 66.3
                    if v == 5: return 152.0
                    if v == 6: return 280.0
                    if v == 7: return 457.0
                    if v == 8: return 690.0
                    if v == 9: return 978.0
                    if v == 10: return 1296.0
                    if v == 11: return 1598.0
                    if v == 12: return 1818.0
                    if v == 13: return 1935.0
                    if v == 14: return 1980.0
                    if v == 15: return 1995.0
                    if v == 16: return 1999.0
                    if v == 17: return 2000.0
                    if v == 18: return 2000.0
                    if v == 19: return 2000.0
                    if v == 20: return 2000.0
                    if v == 21: return 2000.0
                    if v == 22: return 2000.0
                    if v == 23: return 2000.0
                    if v == 24: return 2000.0
                    if v == 25: return 2000.0

                def powertable(U0):
                    return interpolate(floor(U0), power_table(floor(U0)), ceil(U0), power_table(ceil(U0)), U0)

                def interpolate(minx, miny, maxx, maxy, valx):
                    return miny + (maxy - miny) * ((valx - minx) / (maxx - minx))

                def powerstep(U0):
                    if U0 < 4.0:
                        return 0.0
                    elif U0 <= 13.0:
                        return 815.0333333
                    elif U0 <= 25:
                        return 2000.0
                    else:
                        return 0.0

                windrose_angle = []
                windrose_speed = []
                windrose_frequency = []
                for line in windrose:
                    columns = line.split()
                    windrose_angle.append(float(columns[0]))
                    windrose_speed.append(float(columns[1]))
                    windrose_frequency.append(float(columns[2]))

                def analysis():

                    layout_x = []
                    layout_y = []
                    nt = 80
                    for xx in range(nt):
                        l = random()
                        k = random()
                        layout_x.append(5457.0 * l)
                        if layout_x[-1] <= 412.0:
                            layout_y.append(k * 3907.0 + (1.0 - k) * (- 3907.0 / 412.0 * layout_x[-1] + 3907.0))
                        elif layout_x[-1] <= 5040.0:
                            layout_y.append(k * 3907.0)
                        else:
                            layout_y.append(k * (3907.0 / 417.0 * (- layout_x[-1] + 5457.0)))

                    # for U0 in range(4, 20):
                    nt = len(layout_y)  # Number of turbines
                    summation = 0.0
                    r0 = 40.0  # Turbine rotor radius
                    D = 2.0 * r0
                    A = pi * r0 ** 2.0
                    H = 70.0  # Hub height
                    ia = 0.08  # Ambient turbulence intensity according to vanluvanee. 8% on average

                    def deff(U0):
                        return D * sqrt((1.0 + sqrt(1.0 - Ct(U0))) / (2.0 * sqrt(1.0 - Ct(U0))))

                    rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
                    r95 = 0.5 * (rnb + min(H, rnb))

                    def x0(U0):
                        return 9.5 * D / ((2.0 * r95 / deff(U0)) ** 3.0 - 1.0)

                    def c1(U0):
                        return (deff(U0) / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (Ct(U0) * A * x0(U0)) ** (- 5.0 / 6.0)  # Prandtl mixing length

                    def distance_to_front(x, y, theta, r):
                        theta = deg2rad(theta)
                        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

                    for wind in range(0, len(windrose_angle)):
                    # for wind in range(90, 91):
                        # if wind in [100, 133, 271, 280, 313]:
                        #     continue
                        # U1 = windrose_speed[wind]  # Free stream wind speed
                        # U0 = U1 * (70.0 / 10.0) ** 0.11 # Power or log law for wind shear profile
                            # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
                        U0 = 8.5
                        angle = windrose_angle[wind]
                        angle3 = angle + 180.0
                        deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
                        distance = [[0.0 for x in range(2)] for x in range(nt)]

                        U = [U0 for x in range(nt)]
                        total_deficit = [0.0 for x in range(nt)]

                        for tur in range(nt):
                            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle, 100000000.0), tur]
                        distance.sort()

                        for turbine in range(nt):
                            for num in range(turbine):
                                total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
                            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
                            U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
                            flag = [False for x in range(nt)]
                            proportion = [0.0 for x in range(nt)]
                            perpendicular_distance = [0.0 for x in range(nt)]
                            parallel_distance = [0.0 for x in range(nt)]
                            for i in range(turbine + 1, nt):
                                proportion[distance[i][1]], flag[distance[i][1]], perpendicular_distance[distance[i][1]], parallel_distance[distance[i][1]] = wake.determine_if_in_wake(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]], A, c1(U[distance[turbine][1]]), Ct(U[distance[turbine][1]]), angle3, r0, x0(U[distance[turbine][1]]))
                                if parallel_distance[distance[i][1]] > 0.0: ## Add if proportion is 0, skip operation and deficit_matrix == 0.
                                    deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[i][1]] * wake.wake_deficit(U[distance[turbine][1]], Ct(U[distance[turbine][1]]), A, parallel_distance[distance[i][1]] + x0(U[distance[turbine][1]]), perpendicular_distance[distance[i][1]], c1(U[distance[turbine][1]]))
                                else:
                                    deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

                        # for turb in range(nt):
                        #     for i in range(nt):
                        #         output.write('{0:f}\t'.format(deficit_matrix[turb][i]))
                        #     output.write('\n')
                        #     output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(turb, layout_x[turb], layout_y[turb], total_deficit[turb], U[turb], int(angle), power(U[turb])))
                        # output2.write('\n')

                        # for n in range(nt):
                        #     aver[n] += power(U[n]) / 360.0
                        # ---------------------------------------TODO UNCOMMENT -----------------------------------------
                        if withdata:
                            turb_data.write('{0:f} {1:f}\n'.format(angle, power(U[14])))

                        # Farm efficiency
                        profit = 0.0
                        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
                        for l in range(nt):
                            profit += power(U[l])
                        efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
                        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
                    # ---------------------------------------TODO UNCOMMENT ------------------------------------------
                        if withdata:
                            direction.write('{0:f} {1:f}\n'.format(angle, profit))
                        summation += efficiency_proportion[wind]
                    return summation
                    # for n in range(nt):
                    #     turb_data.write('{0:f}\n'.format(aver[n]))

                    # turb_data.close()
                    # output.close()
                    # output2.close()
                    # draw.close()
                    # direction.close()
                    # row.close()
                    # data.close()

                if not withdata:
                    for h in range(20):
                        start_time = time.time()
                        analysis()
                        exe_time.write('{0:f}\n'.format(time.time() - start_time))
                    exe_time.close()
                else:
                    eff.write('{0:f}'.format(analysis()))
                    eff.close()
