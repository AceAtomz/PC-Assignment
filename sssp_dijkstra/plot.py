import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    filepath = './sssp_dijkstra/sssp.out'
    names = ["Serial", "OpenMP T4", "MPI T4"]
    x = list(range(0, 8, 1))  # this is the range of input sizes tested
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            y = list(map(float, line.split(',')))
            plt.plot(x, y, label=names[i])

    plt.xlabel('Input graphs')
    plt.ylabel('Time taken (ms)')
    plt.title('Dijkstras SSSP comparisons')
    plt.legend()
    plt.savefig('./sssp_dijkstra/sssp_comparison.png')
    plt.show()