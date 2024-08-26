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






## Table of Contents
- [Project Overview](#project-overview)
- [Development Environment](#development-environment)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)




## Installation

To set up the development environment for VirusDecode, follow these steps:

1. **Alignment Tool**
    - Download and install [MUSCLE v3.8.1551](https://drive5.com/muscle/):
      ```bash
      mkdir muscle
      cd muscle
      wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
      tar -xvzf muscle_src_3.8.1551.tar.gz
      make
      sudo cp muscle /usr/local/bin/
      cd ..
      rm -rf muscle
      ```

2. **Python Environment**
    - Install Python 3.11.9 and Python 2.7.
    - Set up a virtual environment and install required packages:
      ```bash
      python3.11 -m venv venv
      source venv/bin/activate
      pip install biopython==1.83
      ```
    - If using LinearDesign, ensure Python 2.7 is correctly configured.

3. **Clone the Repository**
    - Clone the project
      ```bash
      git clone https://github.com/VirusDecode-project/VirusDecode.git
      cd VirusDecode
      ```

4. **Backend Setup**
    - Ensure Java 21 (JDK 21) is installed.
    - Install [Spring Boot](https://spring.io/guides/gs/spring-boot/) if necessary.
    - Navigate to the backend directory and run the Spring Boot application:
      ```bash
      cd backend
      cd src/main/resources/bioinformatics/LinearDesign
      make
      cd ../../../../
      ```

5. **Frontend Setup**
    - Install [Node.js and npm](https://nodejs.org/).
    - Navigate to the frontend directory and install dependencies:
      ```bash
      cd frontend
      npm install
      ```





## Usage

To run the VirusDecode application:

1. **Start the Backend Server**
    ```bash
    cd backend
    chmod +x gradlew
    ./gradlew bootRun
    ```

2. **Start the Frontend Development Server**
    ```bash
    cd frontend
    npm start
    ```

3. **Analyzing Virus Sequences**
    - Use the web interface to upload and analyze virus sequences. Follow the on-screen instructions for different types of analysis.
