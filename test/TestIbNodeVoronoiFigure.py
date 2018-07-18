import os

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


bg_gray = '0.75'
bg_line_width = 0.5
font_size = 8

mycol_gr = '#006532'  # Green
mycol_or = '#F39200'  # Orange
mycol_bl = '#2C2E83'  # Blue
mycol_lb = '#0072bd'  # Light blue

# Set LaTeX font rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='computer modern roman')
plt.rc('figure', figsize=[6, 5])  # inches
plt.rc('axes', linewidth=bg_line_width, edgecolor=bg_gray, axisbelow=True)
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.rc('xtick.major', size=0, pad=4)
plt.rc('ytick.major', size=0, pad=4)


data_dir = os.path.join(os.environ['CHASTE_TEST_OUTPUT'], 'VertexIbComp', 'IbNodeVoronoiFigure')

node_csv = os.path.join(data_dir, 'node_locations.csv')

vertex_csvs = []
for f_name in os.listdir(data_dir):
    if f_name.startswith('vertex_locations_'):
        vertex_csvs.append(os.path.join(data_dir, f_name))

output = os.path.join(data_dir, 'VoronoiDiagram.pdf')


node_pos = np.genfromtxt(node_csv, delimiter=',')

vertex_pos_list = []
for f_name in vertex_csvs:
    vertex_pos_list.append(np.genfromtxt(f_name, delimiter=','))


plt.scatter(node_pos[:, 0], node_pos[:, 1])

for cell in vertex_pos_list:
    if len(np.shape(cell)) == 1:
        plt.plot(cell[0], cell[1], markersize=0.3)
    else:
        plt.plot(cell[:, 0], cell[:, 1], marker='o', markersize=0.3, linewidth=0.1)

axes = plt.gca()
axes.axis('equal')

plt.savefig(output)



def rf_samples():
    """
    Generate the 3x3 figure grid for random field sample variance
    """
    rf_sample_files = []

    for f in os.listdir(rf_directory):
        if f.endswith('.rfs'):
            rf_sample_files.append(f)

    for file_num, rf_file in enumerate(sorted(rf_sample_files)):
        length_scale, proportion_of_trace = [float(x) for x in rf_file.replace('.rfs', '').split('-')]

        full_file_path = os.path.join(rf_directory, rf_file)

        grf = np.genfromtxt(full_file_path, delimiter='\n')

        # Check to see if nan values were produced
        grf_nan = grf[np.isnan(grf)]
        if grf_nan.size > 0:
            print(rf_file + " contains nan values")
            continue

        # plt.figure(1, figsize=(6, 6), dpi=72)
        plt.subplot(3, 3, 1+file_num)

        mu, sigma = 0, 1

        plt.axis([-3*sigma, 3*sigma, 0, 0.6])

        if file_num < 3:
            plt.title(r'$t = %s$' % proportion_of_trace, size='xx-large')
        if file_num % 3 == 0:
            plt.ylabel(r'$l = %s$' % length_scale, size='xx-large')

        plt.grid(True)
        n, bins, patches = plt.hist(grf, 50, normed=1, facecolor=mycol_gr, alpha=0.75)

        # Alter axes
        ax = plt.gca()
        ax.grid(which='major', linestyle='dashed', color='0.9')
        ax.yaxis.tick_right()
        if not (file_num + 1) % 3 == 0:
            ax.set_yticklabels([])
        if file_num < 6:
            ax.set_xticklabels([])

        # add a 'best fit' line
        y1 = mlab.normpdf(bins, mu, sigma)
        l1 = plt.plot(bins, y1, linestyle='dashed', color=mycol_bl, linewidth=2)

    plt.savefig(os.path.join(rf_directory, 'RandomFieldSamples.pdf'), bbox_inches='tight', pad_inches=0)
    plt.close()



