
import csv
from pathlib import Path


def main():
    elegans_path = Path("data/stress_proteins/elegans.txt")
    fly_path = Path("data/stress_proteins/fly.txt")
    with open(fly_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            print(line[0].split(":")[-1])


if __name__ == "__main__":
    main()