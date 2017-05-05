import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import click

@click.command()
@click.argument("output")


def main(output):
    print("Output file: ", output)
    df = pd.read_table(output, sep='\t')
    print(df.head(n=5))

    fig = plt.figure(figsize=(16, 8))

    fig.suptitle('Unscented Kalman Filter', fontsize=18, fontweight='bold')
    ax1 = fig.add_subplot(121)
    fig.subplots_adjust(top=0.85)

    ax1.scatter(df["px_ground_truth"], df["py_ground_truth"],alpha=0.4, label = "ground truth")
    ax1.scatter(df["px_state"], df["py_state"],alpha=0.7, marker='x', label = "UKF")
    ax1.set_title('UKF', fontsize=18, fontweight='bold')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend(loc='upper left')

    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(top=0.85)
    ax2.plot(df['NIS'])
    ax2.plot((0, 1400), (7, 7), 'k-', label = '7')
    ax2.legend(loc='upper right')
    ax2.axis([0, 510 , 0, 30])
    ax2.set_title('NIS', fontsize=18, fontweight='bold')
    ax2.set_xlabel('time')
    ax2.set_ylabel('chi2')
    fig.tight_layout()
    plt.savefig('result.png')
    plt.show()

if __name__ == '__main__':
    main()