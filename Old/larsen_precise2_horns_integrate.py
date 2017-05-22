__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Larsen wake model with partial shadowing factor applied to horns rev.
from math import sqrt, log, tan, cos, pi, sin
from numpy import deg2rad
import time
from integrate_simpson import simpson_integrate

output = open('matrix_larsen.dat', 'w')
output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
output2 = open('final_speed_larsen.dat', 'w')
# output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
layout = open('horns_rev.dat', 'r')
windrose = open('horns_rev_windrose2.dat', 'r')
draw = open('draw_horns_rev_larsen.dat', 'w')
draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
# draw2 = open('drawline.dat', 'w')
turb_data = open('x0turb17_larsen.dat', 'w')
direction = open('x0direction_efficiency_larsen.dat', 'w')
direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')
row = open('row_data_larsen.dat', 'w')
data = open('rows_larsen.dat', 'w')

def analysis():

    def wake_radius(c1, ct, A, x):
        return ((35.0 / 2.0 / pi) ** (1.0 / 5.0)) * ((3.0 * c1 ** 2.0) ** (1.0 / 5.0)) * ((ct * A * x) ** (1.0 / 3.0))

    def wake_speed(U0, ct, A, x, r, c1):
        return U0 * (1.0 - ((ct * A * x ** (- 2.0)) ** (1.0 / 3.0)) / 9.0 * (r ** (3.0 / 2.0) * (3.0 * c1 ** 2.0 * ct * A * x) ** (- 1.0 / 2.0) - (35.0 / 2.0 / pi) ** (3.0 / 10.0) * (3.0 * c1 ** 2.0) ** (- 1.0 / 5.0)) ** 2.0)

    def wake_deficit(U0, ct, A, x, r, c1):
        return sqrt((1.0 / A) * simpson_integrate(integrand, 0.0, 40.0, 5, 0.0, 2 * pi, 10)) / U0

    def integrand(dr, th):
        Rw = wake_radius(c1(U[distance[turbine][1]]), Ct(U[distance[turbine][1]]), A, parallel_distance[distance[i][1]] + x0(U[distance[turbine][1]]))
        if sqrt((perpendicular_distance[distance[i][1]] + dr * cos(th)) ** 2.0 + (dr * sin(th)) ** 2.0) <= Rw:
            return dr * (U0 - wake_speed(U[distance[turbine][1]], Ct(U[distance[turbine][1]]), A, parallel_distance[distance[i][1]] + x0(U[distance[turbine][1]]) + dr * cos(th), dr * sin(th), c1(U[distance[turbine][1]]))) ** 2.0
        else:
            return 0.0

    def distance(t1_x, t1_y, t2_x, t2_y):
        return abs(sqrt((t1_x-t2_x) ** 2.0 + (t1_y - t2_y) ** 2.0))

    def crosswind_distance(w_d, t1_x, t1_y, t2_x, t2_y):
        return abs(t1_x * (- sin(w_d)) + t1_y * cos(w_d) - t2_x * (- sin(w_d)) - t2_y * cos(w_d))

    def determine_if_in_wake(xt, yt, xw, yw, A, c1, ct, alpha, r0, x0):  # According to Larsen Model only
        # Eq. of centreline is Y = tan (d) (X - Xt) + Yt
        # Distance from point to line
        alpha = deg2rad(alpha + 180)
        distance_to_centre = abs(- tan(alpha) * xw + yw + tan(alpha) * xt - yt) / sqrt(1.0 + tan(alpha) ** 2.0)
            # print distance_to_centre
        # Coordinates of the intersection between closest path from turbine in wake to centreline.
        X_int = (xw + tan(alpha) * yw + tan(alpha) * (tan(alpha) * xt - yt)) / (tan(alpha) ** 2.0 + 1.0)
        Y_int = (- tan(alpha) * (- xw - tan(alpha) * yw) - tan(alpha) * xt + yt) / (tan(alpha) ** 2.0 + 1.0)
        # Distance from intersection point to turbine
        distance_to_turbine = sqrt((X_int - xt) ** 2.0 + (Y_int - yt) ** 2.0)
        # Radius of wake at that distance
        radius = wake_radius(c1, ct, A, distance_to_turbine + x0)
        # print radius
        if (xw - xt) * cos(alpha) + (yw - yt) * sin(alpha) <= 0.0:
            if abs(radius) >= abs(distance_to_centre):
                if abs(radius) >= abs(distance_to_centre) + r0:
                    value = True
                    return value, distance_to_centre, distance_to_turbine
                elif abs(radius) < abs(distance_to_centre) + r0:
                    value = True
                    return value, distance_to_centre, distance_to_turbine
            elif abs(radius) < abs(distance_to_centre):
                if abs(radius) <= abs(distance_to_centre) - r0:
                    value = False
                    return value, distance_to_centre, distance_to_turbine
                elif abs(radius) > abs(distance_to_centre) - r0:
                    value = True
                    return value, distance_to_centre, distance_to_turbine
        else:
            return 0.0, False, distance_to_centre, distance_to_turbine

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

    def power(U0):
        if U0 < 4.0:
            return 0.0
        elif U0 <= 25.0:
            return 0.0003234808 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 ** 5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
        else:
            return 0.0

    # for U0 in range(4, 20):
    nt = 80  # Number of turbines
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

    aver = [0.0 for x in range(nt)]
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
                flag[distance[i][1]], perpendicular_distance[distance[i][1]], parallel_distance[distance[i][1]] = determine_if_in_wake(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]], A, c1(U[distance[turbine][1]]), Ct(U[distance[turbine][1]]), angle3, r0, x0(U[distance[turbine][1]]))
                if parallel_distance[distance[i][1]] > 0.0 and flag[distance[i][1]] == True: ## Add if proportion is 0, skip operation and deficit_matrix == 0.
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = wake_deficit(U[distance[turbine][1]], Ct(U[distance[turbine][1]]), A, parallel_distance[distance[i][1]] + x0(U[distance[turbine][1]]), perpendicular_distance[distance[i][1]], c1(U[distance[turbine][1]]))
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

        turb_data.write('{0:f}\n'.format(power(U[14])))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        for l in range(nt):
            profit += power(U[l])
        efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        direction.write('{0:f}\n'.format(profit))
        summation += efficiency_proportion[wind]
    print 'total farm efficiency is {0:f} %'.format(summation)
    # for n in range(nt):
    #     turb_data.write('{0:f}\n'.format(aver[n]))

    turb_data.close()
    output.close()
    output2.close()
    draw.close()
    direction.close()
    row.close()
    data.close()

if __name__ == '__main__':
    start_time = time.time()
    analysis()
    print("--- %s seconds ---" % (time.time() - start_time))