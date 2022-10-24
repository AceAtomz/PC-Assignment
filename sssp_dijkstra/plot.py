import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    filepath = './sssp_dijkstra/sssp.out'
    names = ["Serial", "OpenMP T4", "MPI T4"]
    x = [6,8,16,64,128,256,384,512]  # this is the range of input sizes tested
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            y = list(map(float, line.split(',')))
            plt.plot(x, y, label=names[i])

    plt.xlabel('Number of Vertices')
    plt.ylabel('Time taken (ms)')
    plt.title('Serial implementation of Dijkstras SSSP algorithm')
    plt.legend()
    plt.savefig('./support_docs/latex_template/pics/serial_sssp.png')
    plt.show()