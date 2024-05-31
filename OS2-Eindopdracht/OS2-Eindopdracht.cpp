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

// #include <semaphore.h>

const int BUFLEN = 20;   // de lengte van de queue
const int BLOCK_SIZE = 2048;
double ba1, ba2, bb0, bb1, bb2;
double ta1, ta2, tb0, tb1, tb2;


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
}

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
}

class Block {
private:
    std::vector<signed short> data;
    int phase; // 0 = created, 1 = filled, 2 = bassed, 3 = bassed+trebled, 4 = written
    int index;

public:
    Block(int newIndex) : data(BLOCK_SIZE), phase(0), index(newIndex) {}

    void setData(const std::vector<signed short>& newData) {
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

    std::vector<signed short>& getData() {
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
        for (int i = 0; i < BUFLEN; i++) buffer[i] = nullptr;
    }

    void put(Block* block) {
        std::unique_lock<std::mutex> lock(mtx);
        cond_not_full.wait(lock, [this] { return count < BUFLEN; });
        buffer[putpos] = block;
        putpos = (putpos + 1) % BUFLEN;
        count++;
        cond_not_empty.notify_one();
    }

    Block* get() {
        std::unique_lock<std::mutex> lock(mtx);
        cond_not_empty.wait(lock, [this] { return count > 0 || done; });
        if(count == 0 && done) {
            return nullptr;
        }
        Block* block = buffer[getpos];
        getpos = (getpos + 1) % BUFLEN;
        count--;
        cond_not_full.notify_one();
        return block;
    }

    void setDone() {
        std::unique_lock<std::mutex> lock(mtx);
        done = true;
        cond_not_empty.notify_all();
    }

    bool isDone() {
    std::unique_lock<std::mutex> lock(mtx);
    return done;
    }

};

// Fijne globale variabele: de wachtrij
Queue queue;

void Worker() {
    double x[3] = {0, 0, 0};
    double y[3] = {0, 0, 0};
    while (true) {
        Block* block = queue.get();
        if (!block) break;

        if (block->getPhase() == 1) {
            // Apply bass equalizer
            std::vector<signed short>& data = block->getData();
            for (size_t i = 2; i < data.size(); ++i) {
                x[0] = x[1];
                x[1] = x[2];
                x[2] = data[i];
                y[0] = y[1];
                y[1] = y[2];
                y[2] = (bb0/2 * x[0] + bb1/2 * x[1] + bb2/2 * x[2] - ba1 * y[0] - ba2 * y[1]) / 2;

                data[i] = y[2];
            }
            block->setPhase(2);
            std::cout << "Processed block with size " << block->getData().size() << "\n";
            queue.put(block);
        }
        else if (block->getPhase() == 2) {
            // Apply treble equalizer
            std::vector<signed short>& data = block->getData();
            for (size_t i = 2; i < data.size(); ++i) {
                x[0] = x[1];
                x[1] = x[2];
                x[2] = data[i];
                y[0] = y[1];
                y[1] = y[2];
                y[2] = (tb0/2 * x[0] + tb1/2 * x[1] + tb2/2 * x[2] - ta1 * y[0] - ta2 * y[1]) / 2;

                data[i] = y[2];
            }
            block->setPhase(3);
            std::cout << "Processed block with size " << block->getData().size() << "\n";
            queue.put(block);
        }
    }
}

void readPCM(const std::string& inputFile) {
    std::ifstream file(inputFile, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Unable to open input file!" << std::endl;
        return;
    }

    int index = 0;
    while (file.good()) {
        std::vector<signed short> buffer(BLOCK_SIZE / sizeof(signed short));
        file.read(reinterpret_cast<char*>(buffer.data()), BLOCK_SIZE);
        if (file.gcount() == 0) break;

        Block* block = new Block(index++);
        block->setData(buffer);
        block->setPhase(1);
        std::cout << "Read block with size " << block->getData().size() << "\n";
        queue.put(block);
    }
    file.close();
}

void writePCM(const std::string& outputFile) {
    std::ofstream file(outputFile, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Unable to open output file!" << std::endl;
        return;
    }

    while (true) {
        Block* block = queue.get(); 
        std::cout << "Got block from queue\n";
        if (!block && queue.isDone()) break; // Only break if block is nullptr and done is true
        if (!block) continue;
        std::cout << "Block is not empty\n";


        if (block->getPhase() == 3) {
            std::cout << "Block in phase 3\n";
            const std::vector<signed short>& data = block->getData();

            std::cout << "Writing block with size " << data.size() << "\n";
            file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(signed short));
        }
        delete block;
    }
    file.close();
}

int main(int argc, const char* argv[])
{

    if (argc != 6)
    {
        std::cout << argv[0] << std::endl;
        std::cout << "Error: Incorrect number of parameters !" << argc << std::endl;
        std::cout << "Usage: ate -p:<number of threads> -b:<bass intensity> -t:<trebble intensity> <input file> <output file>" << std::endl;
        exit(1);
    }

    

    bassCoefficients(3, &bb0, &bb1, &bb2, &ba1, &ba2);
    trebleCoefficients(-4, &tb0, &tb1, &tb2, &ta1, &ta2);

    int numThreads = std::stoi(argv[1] + 3);
    int bassIntensity = std::stoi(argv[2] + 3);
    int trebleIntensity = std::stoi(argv[3] + 3);
    std::string inputFileStr = argv[4];
    std::string outputFileStr = argv[5];


    std::cout << "Number of threads: " << numThreads << std::endl;
    std::cout << "Bass intensity: " << bassIntensity << std::endl;
    std::cout << "Treble intensity: " << trebleIntensity << std::endl;
    std::cout << "Input file: " << inputFileStr << std::endl;
    std::cout << "Output file: " << outputFileStr << std::endl;
    

    //FILE* inputFile;
    //FILE* outputFile;

    //inputFile = fopen(inputFileStr.c_str(), "rb");
    //outputFile = fopen(outputFileStr.c_str(), "wb");

    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(Worker);
    }
    std::cout << "Reached checkpoint 1...\n";
    readPCM(inputFileStr);
    std::cout << "Reached checkpoint 2...\n";
    for (int i = 0; i < numThreads; ++i) {
        queue.put(nullptr); // Signal workers to stop
    }
    std::cout << "Reached checkpoint 3...\n";
    for (auto& thread : threads) {
        thread.join();
    }
    queue.setDone();
    std::cout << "Reached checkpoint 4...\n";
    writePCM(outputFileStr);

    std::cout << "Reached checkpoint 5...\n";

    return 0;
}

