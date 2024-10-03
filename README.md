# VirusDecode
VirusDecode is an open-source project that integrates various bioinformatics tools to streamline virus sequence analysis and support rapid mRNA vaccine development. Built with a Java Spring backend, React frontend, and Python scripts, it addresses the inefficiencies of traditional genome analysis, enabling quick responses to viral mutations and fostering collaboration in the scientific community.

# Client Access #
For users interested in utilizing the VirusDecode platform, a server has been set up for your convenience. You can access it at the following link:

www.virusdecode.site

This server allows you to explore the features of VirusDecode without needing to set up a local development environment.


## 1. Development Environment
  - Linux (Ubuntu 22.04 LTS)

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
      │   ...
      └── VirusDecode
          ├── backend
          │   ...
          ├── frontend
          │   ...
          └── virusdecode.py
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
        ./gradlew build
        ```

    - Run the Backend
      - After building, run the Spring Boot application:
        ```bash
        java -jar build/libs/virusdecode.jar
        ```

3. **Start the Frontend Development Server**

    - Build the Frontend
      - Install [Node.js and npm](https://nodejs.org/) (nodejs: v18.20.4, npm: v10.7.0)
        ```bash
        curl -fsSL https://deb.nodesource.com/setup_18.x | sudo bash -
        sudo apt-get install -y nodejs
        sudo npm install -g serve
        ```
      - Navigate to the frontend directory and install dependencies:
        ```bash
        cd frontend
        npm install
        npm run build
        ```

    - Run the Frontend
      - Serve the built frontend using the `serve` package:
        ```bash
        serve -s build
        ```

4. **Analyzing Virus Sequences**
   - Use the web interface to upload and analyze virus sequences. Follow the on-screen instructions for different types of analysis.


<br/>
<br/>

## Docker Setup for Easy Usage
Follow these steps to get started with Docker for VirusDecode:
## 1. Initial Setup (Run Once)
1. **Install Docker**
    - Follow the official Docker installation guide: [Docker Installation](https://docs.docker.com/get-docker/)

2. **Clone the VirusDecode repository and Create the Docker containers**
    ```bash
    git clone https://github.com/VirusDecode-project/VirusDecode.git
    cd VirusDecode
    docker-compose build
    ```
## 2. Running the Application
  1. **Run both Backend and Frontend simultaneously**
      ```bash
      docker-compose up
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
