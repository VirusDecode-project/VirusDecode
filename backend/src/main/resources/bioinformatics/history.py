import os
import sys
import shutil
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
option = int(sys.argv[1])

# create history
if option == 1:
    historyName = sys.argv[2]
    new_dir = os.path.join(current_dir, "history/" + historyName)

    os.makedirs(new_dir, exist_ok=True)
    data_dir = os.path.join(current_dir, "data")

    # data 디렉토리의 모든 파일을 새로 생성한 디렉토리로 복사
    for filename in os.listdir(data_dir):
        file_path = os.path.join(data_dir, filename)
        if os.path.isfile(file_path):  # 파일인지 확인
            shutil.copy(file_path, new_dir)

# delete history
elif option == 2:
    historyName = sys.argv[2]
    new_dir = os.path.join(current_dir, "history/" + historyName)

    # history 디렉토리 삭제
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)

# rename history
elif option == 3:
    historyName = sys.argv[2]
    newHistoryName = sys.argv[3]
    new_dir = os.path.join(current_dir, "history/" + historyName)
    new_dir2 = os.path.join(current_dir, "history/" + newHistoryName)

    # history 디렉토리 이름 변경
    if os.path.exists(new_dir):
        os.rename(new_dir, new_dir2)

# get history
elif option == 4:
    historyName = sys.argv[2]
    source_dir = os.path.join(current_dir, "history", historyName)
    data_dir = os.path.join(current_dir, "data")

    # source_dir 디렉토리가 존재하는지 확인
    if os.path.exists(source_dir):
        if os.path.exists(data_dir):
            # data_dir의 모든 파일 삭제
            for filename in os.listdir(data_dir):
                file_path = os.path.join(data_dir, filename)
                if os.path.isfile(file_path):
                    os.remove(file_path)

        # source_dir의 모든 파일을 data_dir로 복사
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            if os.path.isfile(file_path):
                shutil.copy(file_path, data_dir)
# get history list
elif option == 5:
    data_dir = os.path.join(current_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    history_list = os.listdir(os.path.join(current_dir, "history"))
    # JSON으로 변환
    history_json = json.dumps(history_list, indent=4)
    
    # JSON 파일로 저장
    json_file_path = os.path.join(current_dir, "history_list.json")
    with open(json_file_path, 'w', encoding='utf-8') as json_file:
        json_file.write(history_json)
    
#     print(f"History list saved to {json_file_path}")