# VirusDecode

VirusDecode is a project that analyzes virus sequence data and interprets the characteristics of the virus using various bioinformatics tools. This project is implemented using a Java Spring backend, a React frontend, and Python scripts.

## Development Environment

### Operating System
- Linux

### Backend
- Java 21 (JDK 21)
- Spring Boot

### Frontend
- React
- npm

### Programming Language
- Python 3.11.9
- Python 2.7 (required to run LinearDesign)

### Required Packages
- Biopython 1.83

### Alignment Tool
- MUSCLE v3.8.1551 by Robert C. Edgar

## Installation Instructions

### 1. Update system and install necessary packages
```sh
sudo apt update
sudo apt install python2
sudo apt install python3
sudo apt install openjdk-21-jdk
curl -fsSL https://deb.nodesource.com/setup_lts.x | sudo -E bash -
sudo apt install nodejs
pip install biopython
pip install requests
```
### 2. Install necessary alignment tool
```sh
mkdir muscle
cd muscle
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar -xvzf muscle_src_3.8.1551.tar.gz
make
sudo cp muscle /usr/local/bin/
cd ..
rm -rf muscle
```

### 3. Clone and build project
```sh
git clone https://github.com/VirusDecode-project/VirusDecode.git

# Build LinearDesign
cd VirusDecode/backend/src/main/resources/bioinformatics/LinearDesign
make
cd ../../../../../..
```

## 4. Execution
```sh
# Run the frontend
cd frontend
npm install
npm start
```

```sh
# Run the backend
cd backend
chmod +x gradlew
./gradlew bootRun
```

