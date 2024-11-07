# VirusDecode
**VirusDecode** is an open-source project designed to streamline virus sequence analysis by integrating various bioinformatics tools, supporting rapid mRNA vaccine development. With a Java Spring backend, a React frontend, and Biopython for bioinformatics, VirusDecode addresses traditional inefficiencies in genome analysis, enabling quick responses to viral mutations and fostering scientific collaboration.

Please note that two additional open-source tools, MUSCLE for sequence alignment and LinearDesign for mRNA structure prediction, are required but not included in this repository. To set up these tools, follow the provided installation guidelines.


The project structure is organized as follows:
```bash
.
├── frontend
├── backend
└── bioinformatics
    └── analysis
        ├── metadata.py
        ├── alignment.py   # Uses MUSCLE
        ├── mRNA_design.py # Uses LinearDesign
        └── viewer_3d.py
```
The bioinformatics directory contains the core components of the bioinformatics solution, with each script utilizing tools like MUSCLE and LinearDesign via predefined paths. These tools are accessed through the analysis scripts, which interface with the web API to perform sequence analysis and return results to users.

# For Clients #
For clients who interested in using the VirusDecode platform, a server is available at the following link:

[https://virusdecode.com](https://virusdecode.com)


# For Developers #
The following sections are intended for developers who wish to set up the VirusDecode project locally. This guide provides instructions for running the web project and includes steps for using Docker Compose.

## 1. Development Environment
  - OS: Linux (Ubuntu 22.04 LTS)
  - Backend: Java Spring (JDK 21)
  - Frontend: React (Node.js v18.20.4, npm v10.7.0)
  - Database: SQLite (For simplicity in the development environment; the production server uses MySQL)

## 2. Installation
To set up the development environment for VirusDecode, follow these steps:

1. **Python Environment**
    - Install Python 2.7, Python 3.11 and pip
      ```bash
      sudo apt install python2
      sudo apt install python3.11
      sudo apt install python3-pip
      ```
    - Set up a virtual environment and install required packages:
      ```bash
      pip install biopython==1.83
      pip install requests==2.32.3
      ```
    - For using LinearDesign, ensure Python 2.7 is correctly configured.


2. **Alignment Tool**
    - Download and install [MUSCLE v3.8.1551](https://drive5.com/muscle/):
      ```bash
      mkdir muscle
      cd muscle
      wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
      tar -xvzf muscle_src_3.8.1551.tar.gz
      make
      sudo cp muscle /usr/local/bin/
      cd ..
      rm -r muscle
      ```

3. **Clone LinearDesign**
    - VirusDecode uses the LinearDesign open-source software for specific bioinformatics analysis. You must clone LinearDesign separately as per its license terms:
      ```bash
      git clone https://github.com/LinearDesignSoftware/LinearDesign.git
      cd LinearDesign
      make
      cd ..
      ```
    - Note: The LinearDesign code is free for academic, non-profit, and research use. Redistribution of the code with or without modification is not permitted without explicit written permission from the lead corresponding author. If you intend to use this software for commercial purposes, please contact the lead corresponding author for licensing.
    - Using LinearDesign with VirusDecode

      In our project, we used LinearDesign with a specific option, --lambda 3, for optimal performance. The command we used is as follows:
        
        "./lineardesign --lambda 3"
        




## 3. Usage
To run the VirusDecode application:
1. **Clone our Repository**
    - Clone the project
      ```bash
      git clone https://github.com/VirusDecode-project/VirusDecode.git
      cd VirusDecode
      ```
    - Integration with Project Directory: After cloning the repository, it should appear in the project directory structure as follows:
      ``` bash
      .
      ├── LinearDesign
      └── VirusDecode
          ├── frontend
          ├── backend
          ├── bioinformatics
          │   ...
          
      ```

2. **Start the Backend Server**

    - Build the Backend
      - Ensure Java 21 (JDK 21) is installed.
        ```bash
        sudo apt install openjdk-21-jdk 
        ```
      - Navigate to the backend directory and build the Spring Boot application:
        ```bash
        cd backend
        chmod +x gradlew
        ./gradlew bootRun
        ```

3. **Start the Frontend Development Server**

    - Build the Frontend
      - Install [Node.js and npm](https://nodejs.org/) (nodejs: v18.20.4, npm: v10.7.0)
        ```bash
        curl -fsSL https://deb.nodesource.com/setup_18.x | sudo bash -
        sudo apt-get install -y nodejs
        ```
      - Navigate to the frontend directory and install dependencies:
        ```bash
        cd frontend
        npm install
        npm start

4. **Analyzing Virus Sequences**
   - Use the web interface to upload and analyze virus sequences. Follow the on-screen instructions for different types of analysis.



## Docker Setup for Easy Usage
Follow these steps to get started with Docker for VirusDecode:
## 1. Initial Setup (Run Once)
1. **Install Docker**
    - Follow the official Docker installation guide: [Docker Installation](https://docs.docker.com/get-docker/)

2. **Clone the VirusDecode repository and Create the Docker containers**
    ```bash
    git clone https://github.com/VirusDecode-project/VirusDecode.git
    cd VirusDecode
    docker compose build
    ```
## 2. Running the Application
  1. **Run both Backend and Frontend simultaneously**
      ```bash
      docker compose up
      ```
      - Open your web browser and navigate to http://localhost:3000 to access the frontend.
  2. **Run only the Backend on port 8080**
      ```bash
      docker compose up backend
      ```
  3. **Run only the Frontend on port 3000**
      ```bash
      docker compose up frontend
      ```
  - Note: Access the Application
    - Backend will be accessible at: http://localhost:8080
    - Frontend will be accessible at: http://localhost:3000
  4. **After using docker, close docker compose**
      ```bash
      docker compose down
      ```
## User Guide
![guideImg1](https://github.com/user-attachments/assets/09463626-df15-419c-b52a-c581827761c2)
![guideImg2](https://github.com/user-attachments/assets/9dd9a96d-9f9b-45af-a6db-81c406bc3484)
![guideImg3](https://github.com/user-attachments/assets/bfea20d1-1b44-43d5-b2e3-250f62553a72)
![guideImg4](https://github.com/user-attachments/assets/b9abe89c-a27e-41ee-bbb1-2cf2c9722201)


