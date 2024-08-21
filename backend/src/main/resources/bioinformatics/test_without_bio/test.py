import sys
import json

import os



## reference part ( arg1 == 1 )
# input ( argument of python script )
#   - reference id ( string )
# result of bio function
#   - meta data ( json )
# output ( to backend controller )
#   - meta data ( json 형식의 string )
def return_reference_metadata(reference_id):  

    if reference_id == "NC_045512" :
        file_path = './backend/src/main/resources/bioinformatics/test_without_bio/metadata.json'

        with open(file_path, 'r') as file:
            result_json = json.load(file)

        print(json.dumps(result_json, indent=4))
    else :
        print("reference id error")
        

def return_analyze_result(fasta_string):


    # 출력 전달
    file_path = '/Users/delione/Desktop/VD_final/VirusDecode/backend/src/main/resources/bioinformatics/test_without_bio/analyze.json'
    with open(file_path, 'r') as file:
        result_json = json.load(file)
    print(json.dumps(result_json, indent=4))




def main():
    option = int(sys.argv[1])

    # metadata ( inputSeq page )
    if option == 1 : 
        reference_id = sys.argv[2]
        return_reference_metadata(reference_id)
    if option == 2 :
        fasta_string= sys.argv[2]
        return_analyze_result(fasta_string)
    



if __name__ == "__main__":
    main()