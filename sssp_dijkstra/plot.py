import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    filepath = './sssp_dijkstra/sssp.out'
    names = ["Serial", "MPI T2", "MPI T3", "MPI T4", "MPI T8"]
    x = [6,8,16,64,128,256,384,512]  # this is the range of input sizes tested
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            y = list(map(float, line.split(',')))
            plt.plot(x, y, label=names[i])

    plt.xlabel('Number of Vertices')
    plt.ylabel('Time taken (ms)')
    plt.title('MPI implementation of Dijkstras SSSP algorithm')
    plt.legend()
    plt.savefig('./support_docs/latex_template/pics/sssp_mpi.png')
    plt.show()