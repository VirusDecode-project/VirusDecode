# VirusDecode
VirusDecode is an open-source project that integrates various bioinformatics tools to streamline virus sequence analysis and support rapid mRNA vaccine development. Built with a Java Spring backend, React frontend, and Python scripts, it addresses the inefficiencies of traditional genome analysis, enabling quick responses to viral mutations and fostering collaboration in the scientific community.

#### Development Environment
  - Linux

#### Installation

To set up the development environment for VirusDecode, follow these steps:

1. **Python Environment**
    - Install Python 2.7. and Python 3.11
    - Set up a virtual environment and install required packages:
      ```bash
      sudo apt install python2
      sudo apt install python3
      pip install biopython==1.83
      pip install requests
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
      rm -rf muscle
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

#### Usage
To run the VirusDecode application:
1. **Clone our Repository**
    - Clone the project
      ```bash
      git clone https://github.com/VirusDecode-project/VirusDecode.git
      cd VirusDecode
      ```


2. **Start the Backend Server**
    - Ensure Java 21 (JDK 21) is installed.
    - Navigate to the backend directory and run the Spring Boot application:
      ```bash
      cd backend
      chmod +x gradlew
      ./gradlew bootRun
      ```

3. **Start the Frontend Development Server**
    - Install [Node.js and npm](https://nodejs.org/).
    - Navigate to the frontend directory and install dependencies:
      ```bash
      cd frontend
      npm install
      npm start
      ```

3. **Analyzing Virus Sequences**
    - Use the web interface to upload and analyze virus sequences. Follow the on-screen instructions for different types of analysis.
