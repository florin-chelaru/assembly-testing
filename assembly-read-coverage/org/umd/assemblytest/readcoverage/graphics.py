'''
Created on Nov 11, 2013

@author: kzampog
'''

import math
import matplotlib.pyplot as plt

def plot_results(cov_dict, over_dict, under_dict, width, height, filename):
    # m rows, n cols
    n = int(math.ceil(math.sqrt((1.0 * len(cov_dict) * width) / height)))
    n = min(n, len(cov_dict))
    m = int(math.ceil((1.0 * len(cov_dict)) / n))
    global_max = max([max(cov_dict[cid]) for cid in cov_dict])
    plot_max = 1.2 * global_max
    plt.figure(figsize=(width, height))
    cnt = 1
    for cid in cov_dict:
        plt.subplot(m, n, cnt)
        x = range(len(cov_dict[cid]))
        y1 = cov_dict[cid]
        y2 = 0.5 * plot_max * over_dict[cid]
        y3 = plot_max * (1 - 0.5 * under_dict[cid])
        plt.plot(x, y1, x, y2, x, y3)
        plt.ylim(0, plot_max)
        plt.title('Contig ' + cid)
        cnt += 1
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.10, hspace=None)
    if filename != '':
        plt.savefig(filename)
    else:
        plt.show()
    return

'''
if __name__ == '__main__':
    cov_dict = {}
    cov_dict['a'] = range(10)
    cov_dict['b'] = range(10)[::-1]
    res_dict = {}
    res_dict['a'] = range(10)[::-1]
    res_dict['b'] = range(10)

    print cov_dict
    print res_dict
    plot_results(cov_dict, res_dict, 8, 5, '')
'''