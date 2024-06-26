//
//  main.cpp
//  ate
//
//  Created by Koen van Brero on 12-05-14.
//  Copyright (c) 2014 Koen van Brero. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <thread>
#include <mutex>
#include <math.h>
#include <queue>
#include <utility> 
#include <array>

#define SIZE 2048
// #include <semaphore.h>

const int BUFLEN = 20000;   // de lengte van de queue
const int BLOCK_SIZE = 2048; // de grote van blocken
double ba1, ba2, bb0, bb1, bb2;
double ta1, ta2, tb0, tb1, tb2;

// Bereken de bass coefficients
void bassCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2)
{
    double frequency = 330;
    double qFactor = 0.5;
    double gain = intensity;
    double sampleRate = 44100;

    double pi = 4.0 * atan(1);
    double a = pow(10.0, gain / 40);
    double w0 = 2 * pi * frequency / sampleRate;
    double alpha = sin(w0) / (2.0 * qFactor);
    double a0 = (a + 1) + (a - 1) * cos(w0) + 2.0 * sqrt(a) * alpha;

    *a1 = -(-2.0 * ((a - 1) + (a + 1) * cos(w0))) / a0;
    *a2 = -((a + 1) + (a - 1) * cos(w0) - 2.0 * sqrt(a) * alpha) / a0;

    *b0 = (a * ((a + 1) - (a - 1) * cos(w0) + 2.0 * sqrt(a) * alpha)) / a0;
    *b1 = (2 * a * ((a - 1) - (a + 1) * cos(w0))) / a0;
    *b2 = (a * ((a + 1) - (a - 1) * cos(w0) - 2.0 * sqrt(a) * alpha)) / a0;
    std::cout << "bb0: " << *b0 << " bb1: " << *b1 << " bb2: " << *b2 << " ba1: " << *a1 << " ba2: " << *a2 << std::endl;

}

// Bereken de treble coefficients
void trebleCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2)
{
    double frequency = 3300;
    double qFactor = 0.5;
    double gain = intensity;
    double sampleRate = 44100;

    double pi = 4.0 * atan(1);
    double a = pow(10.0, gain / 40);
    double w0 = 2 * pi * frequency / sampleRate;
    double alpha = sin(w0) / (2.0 * qFactor);
    double a0 = (a + 1) - (a - 1) * cos(w0) + 2.0 * sqrt(a) * alpha;

    *a1 = -(2.0 * ((a - 1) - (a + 1) * cos(w0))) / a0;
    *a2 = -((a + 1) - (a - 1) * cos(w0) - 2.0 * sqrt(a) * alpha) / a0;

    *b0 = (a * ((a + 1) + (a - 1) * cos(w0) + 2.0 * sqrt(a) * alpha)) / a0;
    *b1 = (-2.0 * a * ((a - 1) + (a + 1) * cos(w0))) / a0;
    *b2 = (a * ((a + 1) + (a - 1) * cos(w0) - 2.0 * sqrt(a) * alpha)) / a0;

    std::cout << "tb0: " << *b0 << " tb1: " << *b1 << " tb2: " << *b2 << " ta1: " << *a1 << "ta2: " << *a2 << std::endl;
}

// Een blok met data
class Block {
private:
    std::array<signed short, SIZE / 2> data; // 1024 samples
    int phase; // 0 = created, 1 = filled, 2 = bassed, 3 = bassed+trebled, 4 = written
    int index;

public:
    Block(int newIndex) : phase(0), index(newIndex) {
        data.fill(0);
    }

    void setData(const std::array<signed short, SIZE / 2> newData) {
        data = newData;
    }

    void setPhase(int newPhase) {
        phase = newPhase;
    }

    int getPhase() const {
        return phase;
    }

    int getIndex() const {
        return index;
    }

    std::array<signed short, SIZE / 2> getData() {
        return data;
    }
};


// De ronde wachtrij
// Die moet je dus thread-safe maken!
class Queue {
private:
    Block* buffer[BUFLEN];
    int getpos, putpos, count;
    std::mutex mtx;
    std::condition_variable cond_not_empty;
    std::condition_variable cond_not_full;
    bool done;

public:
    Queue() : getpos(0), putpos(0), count(0), done(false) { 
        for (int i = 0; i < BUFLEN; i++) buffer[i] = nullptr;           // Vul de buffer met nullptrs
    }

    void put(Block* block) {
        std::unique_lock<std::mutex> lock(mtx);                         // Vergrendel de mutex
        cond_not_full.wait(lock, [this] { return count < BUFLEN; });    // Wacht tot de buffer niet vol is
        buffer[putpos] = block;                                         
        putpos = (putpos + 1) % BUFLEN;							       
        count++;													
        cond_not_empty.notify_one();								    // Signaleer dat de buffer niet leeg is
    }

    Block* get() {
        std::unique_lock<std::mutex> lock(mtx);                         // Vergrendel de mutex
        cond_not_empty.wait(lock, [this] { return count > 0 || done; });// Wacht tot de buffer niet leeg is of done is
        if(count == 0 && done) {                                       
            return nullptr;
        }
        Block* block = buffer[getpos];                                  
        getpos = (getpos + 1) % BUFLEN;								    
        count--;
        cond_not_full.notify_one();									    // Signaleer dat de buffer niet vol is
        return block;
    }

    void setDone() {
        std::unique_lock<std::mutex> lock(mtx);                         // Vergrendel de mutex
        done = true;
        cond_not_empty.notify_all();									// Signaleer dat de buffer niet leeg is
    }

    bool isDone() {
    std::unique_lock<std::mutex> lock(mtx);
    return done;
    }

};

// Fijne globale variabele: de wachtrij
Queue queue;

// De worker thread om de blocks te bewerken
void Worker() {
    std::array<signed short, SIZE / 2> y = { 0 };                         
    signed short ctb0 = static_cast<signed short>(ta1), ctb1 = static_cast<signed short>(ta2),                              // Coefficients voor de treble
        ctb2 = static_cast<signed short>(tb0), cta1 = static_cast<signed short>(tb1), cta2 = static_cast<signed short>(tb2),
        cbb0 = static_cast<signed short>(ba1), cbb1 = static_cast<signed short>(ba2),							            // Coefficients voor de bass
        cbb2 = static_cast<signed short>(bb0), cba1 = static_cast<signed short>(bb1), cba2 = static_cast<signed short>(bb2);

    std::cout << "ctb0: " << ctb0 << " ctb1: " << ctb1 << " ctb2: " << ctb2 << " cta1: " << cta1 << " cta2: " << cta2 << std::endl;
    std::cout << "cbb0: " << cbb0 << " cbb1: " << cbb1 << " cbb2: " << cbb2 << " cba1: " << cba1 << " cba2: " << cba2 << std::endl;
    
    while (true) {
        Block* block = queue.get();
        if (!block) break;

        if (block->getPhase() == 1) { // Als block gevuld is
            y.fill(0);
            // Apply bass equalizer
            std::array<signed short, SIZE / 2> x = block->getData(); // haal de data uit het block
            for (int n = 0; n < y.size(); n++) {
                // Bewerk de data
                if (n >= 2) { 
                    y[n] = cbb0 * x[n] + cbb1 * x[n - 1] + cbb2 * x[n - 2] + cba1 * y[n - 1] + cba2 * y[n - 2];
                }
                else if (n == 1) {
                    y[n] = cbb0 * x[n] + cbb1 * x[n - 1] + cba1 * y[n - 1];
                }
                else {
                    y[n] = cbb0 * x[n];
                }
            }
            block->setData(y);  // Zet de data terug in het block
            block->setPhase(2); // Zet de fase naar bassed
            std::cout << "Processed BASS, block with size " << block->getData().size() << "\n";
            queue.put(block);   // Zet het block terug in de queue
        }
        else if (block->getPhase() == 2) {
            // Apply treble equalizer
            std::array<signed short, SIZE / 2> x = block->getData(); // haal de data uit het block
            for (int n = 0; n < y.size(); n++) {
                // Bewerk de data
                if (n >= 2) {
                    y[n] = ctb0 * x[n] + ctb1 * x[n - 1] + ctb2 * x[n - 2] + cta1 * y[n - 1] + cta2 * y[n - 2];
                }
                else if (n == 1) {
                    y[n] = ctb0 * x[n] + ctb1 * x[n - 1] + cta1 * y[n - 1];
                }
                else {
                    y[n] = ctb0 * x[n];
                }
            }
            block->setData(y);  // Zet de data terug in het block
            block->setPhase(3); // Zet de fase naar bassed+trebled
            std::cout << "Processed TREBLE, block with size " << block->getData().size() << "\n";
            queue.put(block); // Zet het block terug in de queue
        }
        else // Als block in fase 3 is en om te verkomen dat het block uit de queue wordt gehaald en niet wordt teruggezet
        {
            			queue.put(block);
        }
    }
}

// Lees de PCM data
void readPCM(const std::string& inputFile) {
    std::ifstream file(inputFile, std::ios::binary); // Open het bestand
    if (!file.is_open()) {
        std::cerr << "Unable to open input file!" << std::endl;
        return;
    }

    int index = 0;
    while (file.good()) {
        std::array<signed short, SIZE / 2> buffer;
        file.read(reinterpret_cast<char*>(buffer.data()), BLOCK_SIZE); // Lees de data in de buffer 
        if (file.gcount() == 0) break; // Als er niks is gelezen, stop

        Block* block = new Block(index++); // Maak een nieuw block aan
        block->setData(buffer);
        block->setPhase(1);
        std::cout << "Read block with size " << block->getData().size() << "\n";
        queue.put(block); // Zet het block in de queue
    }
    file.close();
}

// Schrijf de PCM data
void writePCM(const std::string& outputFile) {
    std::ofstream file(outputFile, std::ios::binary); // Open het bestand
    if (!file.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        return;
    }

    std::priority_queue<std::pair<int, Block*>, std::vector<std::pair<int, Block*>>, std::greater<>> blockQueue; // Maak een priority queue aan

    // Zet de blocks in de priority queue
    while (true) {
        Block* block = queue.get();
        std::cout << "Got block from queue\n";
        if (!block && queue.isDone()) break; // Als er geen block is en de queue is klaar, stop
        if (!block) continue;
        std::cout << "Block is not empty\n";

        blockQueue.push(std::make_pair(block->getIndex(), block)); // Zet het block in de queue 
    }

    // Schrijf de blocks naar het bestand
    while (!blockQueue.empty()) {
        Block* block = blockQueue.top().second; // Haal het block uit de queue
        blockQueue.pop();

        // Als het block in bassed+trebled fase is
        if (block->getPhase() == 3) {
            std::cout << "Block in phase 3\n";
            std::array<signed short, SIZE / 2> data = block->getData();
            std::cout << "Writing block with size " << data.size() << " and index " << block->getIndex() << "\n";
            file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(signed short)); // Schrijf de data naar het bestand
        }
        delete block;
    }
    file.close();
}

// Main functie
int main(int argc, const char* argv[])
{
    // Check of er genoeg argumenten zijn
    if (argc != 6)
    {
        std::cout << argv[0] << std::endl;
        std::cout << "Error: Incorrect number of parameters !" << argc << std::endl;
        std::cout << "Usage: ate -p:<number of threads> -b:<bass intensity> -t:<trebble intensity> <input file> <output file>" << std::endl;
        exit(1);
    }

    

    // Parse de argumenten
    int numThreads = std::stoi(argv[1] + 3);
    int bassIntensity = std::stoi(argv[2] + 3);
    int trebleIntensity = std::stoi(argv[3] + 3);
    std::string inputFileStr = argv[4];
    std::string outputFileStr = argv[5];

    // Bereken de coefficients
    bassCoefficients(bassIntensity, &bb0, &bb1, &bb2, &ba1, &ba2);
    trebleCoefficients(trebleIntensity, &tb0, &tb1, &tb2, &ta1, &ta2);

    std::cout << "Number of threads: " << numThreads << std::endl;
    std::cout << "Bass intensity: " << bassIntensity << std::endl;
    std::cout << "Treble intensity: " << trebleIntensity << std::endl;
    std::cout << "Input file: " << inputFileStr << std::endl;
    std::cout << "Output file: " << outputFileStr << std::endl;

    // Maak de threads aan
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(Worker);
    }
    std::cout << "Reached checkpoint 1...\n"; // Thread aanmaken is gelukt
    readPCM(inputFileStr);
    std::cout << "Reached checkpoint 2...\n"; // PCM lezen is gelukt
    for (int i = 0; i < numThreads; ++i) {
        queue.put(nullptr); // Signal workers to stop
    }
    std::cout << "Reached checkpoint 3...\n"; // Queue is klaar 
    for (auto& thread : threads) {
       thread.join();
    }
    queue.setDone();
    std::cout << "Reached checkpoint 4...\n"; // Threads zijn klaar 
    writePCM(outputFileStr);

    std::cout << "Reached checkpoint 5...\n"; // PCM schrijven is gelukt

    return 0;
}

