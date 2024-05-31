#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>
#include <condition_variable>

const int BUFLEN = 20;
const int BLOCK_SIZE = 2048;

void bassCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2);
void trebleCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2);

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

class Queue {
private:
    Block* buffer[BUFLEN];
    int getpos, putpos, count;
    std::mutex mtx;
    std::condition_variable cond_not_empty;
    std::condition_variable cond_not_full;

public:
    Queue() : getpos(0), putpos(0), count(0) {
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
        cond_not_empty.wait(lock, [this] { return count > 0; });
        Block* block = buffer[getpos];
        getpos = (getpos + 1) % BUFLEN;
        count--;
        cond_not_full.notify_one();
        return block;
    }
};

Queue queue;

double ba1, ba2, bb0, bb1, bb2;
double ta1, ta2, tb0, tb1, tb2;

void worker() {
    while (true) {
        Block* block = queue.get();
        if (!block) break;

        if (block->getPhase() == 1) {
            // Apply bass equalizer
            std::vector<signed short>& data = block->getData();
            for (size_t i = 2; i < data.size(); ++i) {
                // Example processing logic
                // data[i] = ...; // Apply biquad filter for bass
            }
            block->setPhase(2);
            queue.put(block);
        } else if (block->getPhase() == 2) {
            // Apply treble equalizer
            std::vector<signed short>& data = block->getData();
            for (size_t i = 2; i < data.size(); ++i) {
                // Example processing logic
                // data[i] = ...; // Apply biquad filter for treble
            }
            block->setPhase(3);
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
        if (!block || block->getPhase() != 3) break;

        const std::vector<signed short>& data = block->getData();
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(signed short));
        delete block;
    }
    file.close();
}

int main(int argc, const char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: ate -p:<number of threads> -b:<bass intensity> -t:<treble intensity> <input file> <output file>" << std::endl;
        return 1;
    }

    int numThreads = std::stoi(argv[1] + 3);
    int bassIntensity = std::stoi(argv[2] + 3);
    int trebleIntensity = std::stoi(argv[3] + 3);
    std::string inputFile = argv[4];
    std::string outputFile = argv[5];

    bassCoefficients(bassIntensity, &bb0, &bb1, &bb2, &ba1, &ba2);
    trebleCoefficients(trebleIntensity, &tb0, &tb1, &tb2, &ta1, &ta2);

    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(worker);
    }

    readPCM(inputFile);
    for (int i = 0; i < numThreads; ++i) {
        queue.put(nullptr); // Signal workers to stop
    }

    for (auto& thread : threads) {
        thread.join();
    }

    writePCM(outputFile);

    return 0;
}

void bassCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2) {
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

void trebleCoefficients(int intensity, double* b0, double* b1, double* b2, double* a1, double* a2) {
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
