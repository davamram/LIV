import argparse

def extract_second_column(file_path, col):
    second_column_numbers = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) == 5:
                second_column_numbers.append(columns[col])
    return ", ".join(second_column_numbers[-16:])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Y values")
    parser.add_argument("file_name", help="Name of the input yoda file")
    args = parser.parse_args()

    file_path="/home/amram/Documents/LorentzPhotons/Rivet/LIV/Plots/"+args.file_name+"/TEST_ANALYSIS/d01-x01-y01.dat"
    result = extract_second_column(file_path, 2)
    print("values")
    print(result)
    print("err :")
    result = extract_second_column(file_path, 3)
    print(result)