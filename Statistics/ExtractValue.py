import argparse

def extract_second_column(file_path):
    second_column_numbers = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) == 5:  # Vérifie si la ligne a au moins 2 colonnes
                second_column_numbers.append(columns[2])
    return ", ".join(second_column_numbers)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Y values")
    parser.add_argument("file_name", help="Name of the input yoda file")
    args = parser.parse_args()

    file_path="/home/amram/Documents/LorentzPhotons/Rivet/LIV/Plots/"+args.file_name+"/TEST_ANALYSIS/d01-x01-y01.dat"
    result = extract_second_column(file_path)
    print(result)