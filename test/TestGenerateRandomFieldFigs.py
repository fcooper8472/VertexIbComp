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


rf_directory = os.path.join(os.environ['CHASTE_TEST_OUTPUT'], 'VertexIbComp', 'RandomFields')


def main():
    rf_samples()
    # rf_viz()
    # rf_eigs()


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


def rf_viz():
    """
    Generate the 3x3 figure grid for random field visualisations
    """
    rf_viz_files = []

    for f in os.listdir(rf_directory):
        if f.endswith('.rfv'):
            rf_viz_files.append(f)

    for file_num, rf_file in enumerate(sorted(rf_viz_files)):
        length_scale, num = rf_file.replace('.rfv', '').split('-')
        length_scale = float(length_scale)
        num = int(num)

        full_file_path = os.path.join(rf_directory, rf_file)

        grf = np.genfromtxt(full_file_path, delimiter=',')
        points = grf[:, :2]  # first 2 cols are the (x,y) pairs
        values = grf[:, 2]   # third col is the z-value

        # Check to see if nan values were produced
        grf_nan = values[np.isnan(values)]
        if grf_nan.size > 0:
            print(rf_file + " contains nan values")
            continue

        # Tile the points, 3x3, to properly represent the periodic field
        x_periodic = np.array(points[:, 0])
        y_periodic = np.array(points[:, 1])
        values_periodic = np.array(values)

        offsets = [[-20, 20], [0, 20], [20, 20], [-20, 0], [20, 0], [-20, -20], [0, -20], [20, -20]]
        # Tile the points, 3x3, to properly represent the periodic field
        for offset in offsets:
            for i in range(points.shape[0]):
                x_periodic = np.append(x_periodic, x_periodic[i] + offset[0])
                y_periodic = np.append(y_periodic, y_periodic[i] + offset[1])
            values_periodic = np.append(values_periodic, values)

        x_grid = np.linspace(0, 20, 50)
        y_grid = np.linspace(0, 20, 50)
        grid_data = mlab.griddata(x_periodic, y_periodic, values_periodic, x_grid, y_grid, interp='linear')

        plt.subplot(3, 3, 1+file_num)
        plt.imshow(grid_data, interpolation='gaussian', extent=(0, 20, 0, 20), origin='lower')
        # plt.scatter(points[:, 0], points[:, 1], 0.1, marker='.', color='black')

        plt.xlim(0, 20)
        plt.ylim(0, 20)

        if file_num < 3:
            plt.title(r'$%s$' % num, size='xx-large')
        if file_num % 3 == 0:
            plt.ylabel(r'$l = %s$' % length_scale, size='xx-large')

        # Alter axes
        ax = plt.gca()
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    plt.savefig(os.path.join(rf_directory, 'RandomFieldViz.pdf'), bbox_inches='tight', pad_inches=0)
    plt.close()


def rf_eigs():
    """
    Generate the 3x1 figure grid for number of eigenvalues vs trace proportion
    """
    rf_eigs_files = []

    for f in os.listdir(rf_directory):
        if f.endswith('.rfe'):
            rf_eigs_files.append(f)

    list_of_trace_props = []
    list_of_eval_props = []
    list_of_length_scales = []

    for file_num, rf_file in enumerate(sorted(rf_eigs_files)):
        length_scale = float(rf_file.replace('.rfe', ''))

        full_file_path = os.path.join(rf_directory, rf_file)

        data = np.genfromtxt(full_file_path, delimiter=',')

        list_of_trace_props.append(data[:, 0])
        list_of_eval_props.append(data[:, 1])
        list_of_length_scales.append(length_scale)

    highlights = [0, 4, 8]
    labels = ['(a)', '(b)', '(c)']
    for i, highlight in enumerate(highlights):
        plt.figure(1, figsize=(6, 1.8))
        plt.subplot(1, 3, 1+i)
        plt.title(r'$l = %s$' % list_of_length_scales[highlight], size=10)
        plt.xlabel(r'%s' % labels[i], size=10)

        for line_num in range(len(list_of_trace_props)):
            plt.plot(list_of_trace_props[line_num], list_of_eval_props[line_num], color=bg_gray, linewidth=2.0)
        plt.plot(list_of_trace_props[highlight], list_of_eval_props[highlight], color=mycol_gr, linewidth=3.0)

        ax = plt.gca()
        x_ticks = [0.6, 0.7, 0.8, 0.9, 1.0]
        y_ticks = [0.0, 0.2, 0.4, 0.6, 0.8]

        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)

        if i == 1:
            ax.set_yticklabels([])
        if i == 2:
            ax.yaxis.tick_right()

        ax.grid(which='major', linestyle='dashed', color='0.9')

    # plt.subplots_adjust(top=0.8)
    plt.savefig(os.path.join(rf_directory, 'RandomFieldEvals.pdf'), bbox_inches='tight', pad_inches=0)
    plt.close()


if __name__ == "__main__":
    main()
