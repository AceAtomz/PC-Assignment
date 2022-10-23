import matplotlib.pyplot as plt

if __name__ == '__main__':

    filepath = './scan.out'
    names = ["Serial", "OpenMP" , "MPI"]  
    x = list(range(1000, 20001, 1000))  # this is the range of input sizes tested
    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            y = list(map(int, line.split(',')))
            plt.plot(x, y, label=names[i])

    plt.xlabel('Input size')
    plt.ylabel('Time taken (ns)')
    plt.title('Comparison of Prefix Sum Procedures')
    plt.legend()
    plt.savefig('scan_comparison.pdf')
    plt.show()